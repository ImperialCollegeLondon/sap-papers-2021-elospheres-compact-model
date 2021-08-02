function[] = mvdr_20210601_compact(...
    in_wav_file_path ...          % path to received signal
    ,out_wav_file_path ...        % path where enhanced signal(s) should be written
    ,ht_file_path ...             % ignored - path to head tracker signal
    ,listener_characteristics ... % ignored - path to listener characteristics profile
    ,in_params ...                % structure with config information for test bench
    ,oracle_data ...              % structure containing information and/or paths to information which would not be available in real-world application but may be useful in development
    ,saved_data_dir ...           % path to folder where any supplementary/intermediate results should be written
    ,temp_data_dir ...            % path to folder where temporary files should be written - normally empty but may be useful for debugging
    )
% This is an updated version of 
% mvdr_20190128_01_binaural_steered_oracle_circiso_sensor_pw 
% changes:
% - new minimisation (Mike) put in
% - adapted for eigenmike
% - uses parfor loops for speed
% - allows a known noise source direction to be specified and does the
% necessary calculations
% - allows any combination of free and fixed directional components (max. 1
% free)
% fixed the issues saving the direction

% mvdr_20200811_01_eig_parameterised_fastsearch implements an mvdr
% beamformer where the PSD/covariance matrix is modelled as a combination of
% circularly isotropic noise, spatially white noise and a planewave interferer 
% from an estimated direction.
% The look direction is steered towards the known source direction in response to
% array rotation. Unlike other beamformers, this one works directly in frequency
% domain ignoring circular convolution issues for simplicity
%
% Usage:
%  mvdr_20190128_01_binaural_steered_oracle_circiso_sensor_pw(...
%      in_wav_file_path, out_wav_file_path, ht_file_path, listener_characteristics, ...
%      in_params, oracle_data, saved_data_dir, temp_data_dir)
%
% Outputs:
%  None
%
% Inputs:
%  in_wav_file_path:
%      received signal
%  out_wav_file_path:
%      path where enhanced signal(s) should be written
%  ht_file_path:
%      head tracker signal
%  listener_characteristics:
%      [not yet implemented]
%  in_params:
%      struct
%  oracle_data:
%      sruct containing information which would not be available in a
%      real-world application but may be useful in development
%  saved_data_dir:
%      folder for storing intermediate data for exchange between modules or
%      supplementary results which should not be included in metrics.
%      May be empty if test bench dictates that this data should not be retained
%      after the experiment, in which case intermediate data should be
%      saved to temporary storage (see tempname).
%  temp_data_dir:
%      folder for storing data which is not required once processing is
%      complete. Normally empty, in which case use standard temporary directory
%      (see tempname) but useful for debugging.
%
% Alastair Moore, January 2019
% Rebecca Vos, May 2020
% Sina Hafezi, June 2021

fprintf('start time: %s\n', datestr(now,'HH:MM:SS'));
%% validate the input
% not really required in this simple function but will use it as a template

% in_wav_file_path
if ~exist(in_wav_file_path,'file'), error('Couldn''t find %s',in_wav_file_path), end

% final_data_dir
final_data_dir = fileparts(out_wav_file_path);
check_output_dir_exists(final_data_dir)

% ht_file_path
if ~exist(ht_file_path,'file'), error('Couldn''t find %s',ht_file_path), end

% listener_characteristics (haven't specified this yet)
% if not present, assume default listener?

% oracle_data (struct)
% check the required field(s) are present
required_oracle_data = {...
    'target_position_csv_path'...            % original doa inclination
    ,'ema_fcn'...
    };
for ireq = 1:length(required_oracle_data)
    if ~isfield(oracle_data,required_oracle_data{ireq})
        error('oracle_data is missing %s field',required_oracle_data{ireq})
    end
end

% saved_data_dir (string)
if nargin<7 || isempty(saved_data_dir)
    error('saved_data_dir should not be empty - want to retain the data')
end
check_output_dir_exists(saved_data_dir);

% temp_data_dir (string)
% if nargin<8 || isempty(temp_data_dir)
%     temp_data_dir = tempname;
% end
% check_output_dir_exists(temp_data_dir);

%% defualt parameters
params.fs = [];
params.c = soundspeed();
params.frame_duration = 0.016;
params.frame_overlap_frac = 0.5;
params.max_condition_number = 1000;
params.mvdr_fcn = @fcn_20170222_01_mvdr_weights;
params.mvdr_debug_level = 0;%

params.do_export_filters_as_wav_per_frame = 0;
params.filename_prefix = 'mvdr_moving_target';


% specific to noise estimation
params.noise_only_init_time = 0.1; % seconds
params.time_constant = 0.05; %seconds
params.nFree = 1;
params.nFixed = 0;
params.circ_az = deg2rad((0:5:355).' );         % assuming circular diffuse field
params.max_planewave_to_circiso_pow_ratio = 4;

allowed_values.nFree = [0, 1];
allowed_values.nFixed = 0;      % this version does not allow fixed direction components


%% update with input parameters
if ~isempty(in_params)
    params = override_valid_fields(params,in_params,allowed_values);
end

% safety check
if params.nFixed ~=0
    error('Fixed components not supported')
end

%% read in data and do pre-processing
fprintf('Loading wav file...\n')
[y, fs_in] = readwav(in_wav_file_path,'g'); %[nSamples,nChans]
if isempty(params.fs)
    params.fs = fs_in;
else
    y = resample(y,params.fs,fs_in);
end

% optional audio trimming (if longer than 20 s)
%%{
if size(y,1)>round(20*params.fs)
    y=y(1:round(20*params.fs),:);
end
%}

% nReqOutSamples = size(y,1);
nChans = size(y,2);

% we now now the sample rate. Do the stft to determine the frame times
fs = params.fs;
stft_params = convert_in_params_to_stft_params(params);

fprintf('Doing STFT...\n')
[Y,x_tail_anal,pm] = stft_v2('fwd', y, ...
    stft_params.win_anal, ...
    stft_params.sig_frame_inc,...
    stft_params.nfft,...
    stft_params.fs);
[nFreq,nChan,nFrames] = size(Y);


%% instantiate the microphone array
% - from this we obtain the impulse responses for different directions
ema = oracle_data.ema_fcn(); % instatiate the array class
ema.prepareData(fs);

%% read in rotation/position data and get values at each frame time
fprintf('Interpolating array rotation...\n')
% array rotation
% - stored as roll pitch yaw but we want it as quaternions
[rpy,fs_rpy] = load_ht_rpy_as_wav(ht_file_path);
t_rpy = (0:size(rpy,1)-1).' ./fs_rpy;
qr_in = v_roteu2qr('xyz',rpy.'); % [4 nSamples]: [w x y z]
qr_w = griddedInterpolant(t_rpy,qr_in(1,:),'linear');
qr_x = griddedInterpolant(t_rpy,qr_in(2,:),'linear');
qr_y = griddedInterpolant(t_rpy,qr_in(3,:),'linear');
qr_z = griddedInterpolant(t_rpy,qr_in(4,:),'linear');
qr = zeros(4,nFrames);
qr(1,:) = qr_w(pm.t);
qr(2,:) = qr_x(pm.t);
qr(3,:) = qr_y(pm.t);
qr(4,:) = qr_z(pm.t);
clear rpy fs_rpy qr_in qr_x qr_y qr_z

% target direction
fprintf('Interpolating target position...\n')
% - obtained from cartesian position
target_pos = readtable(oracle_data.target_position_csv_path,...
    'ReadVariableNames',0);
if size(target_pos,2)~=4
    error('csv file %s does not have 4 columns',...
        oracle_data.target_position_csv_path)
end
target_pos.Properties.VariableNames = {'t','x','y','z'};
target_x = griddedInterpolant(target_pos.t,target_pos.x,'linear');
target_y = griddedInterpolant(target_pos.t,target_pos.y,'linear');
target_z = griddedInterpolant(target_pos.t,target_pos.z,'linear');
target_doa_vec = zeros(nFrames,3);
target_doa_vec(:,1) = target_x(pm.t);
target_doa_vec(:,2) = target_y(pm.t);
target_doa_vec(:,3) = target_z(pm.t);
clear target_pos target_x target_y target_z



%% precompute model covariance bases
fprintf('Computing basis covariance matrices...\n')
model_params.fs = params.fs;
model_params.stft_params = stft_params;
model_params.ema = ema;


% circularly/cylindrical isotropic
%{
model_params.doa_circ.az = params.circ_az;
model_params.doa_circ.inc = deg2rad( 90 * ones(size(model_params.doa_circ.az)));

model = fcn_20210215_ahm_02_get_compact_model_bases(model_params);
%}

% spherical isotropic
%%{
aa=params.circ_az;
bb=deg2rad([80:5:100].');
[AA,BB]=meshgrid(aa,bb);
all_look_doas=reshape(cat(2,AA',BB'),[],2);

[xxx yyy zzz]=mysph2cart(all_look_doas(:,1),all_look_doas(:,2),ones(size(all_look_doas,1),1));
look_dirs_cart=[xxx yyy zzz];

% get moving target direction to remove nearby estimated doas
ang_spacing_threshold=10; % deg around which to avoid search in optimisation (preventing potential mis-steering)
dirs_to_remove=cell(nFrames,1);
for iframe = 1:nFrames
ema.setPoseQuaternion(qr(:,iframe));
[target_az, target_inc,~] = mycart2sph(target_doa_vec(iframe,:));
[target_az,target_inc] = ema.planeWaveDoaWrtArray(target_az,target_inc);
[xxx yyy zzz]=mysph2cart(target_az,target_inc,1);
steer_cart=[xxx yyy zzz];
look_dirs_ang=cellfun(@(xxr) rad2deg(atan2(norm(cross(xxr,steer_cart)),dot(xxr,steer_cart))),table2cell(table(look_dirs_cart)));
dirs_to_remove{iframe}=find(look_dirs_ang<=ang_spacing_threshold);
end
model_params.doas = all_look_doas;
model_params.iso_grid.azi_step=10;
model_params.iso_grid.inc_step=10;
model = fcn_20210506_get_compact_model_bases(model_params);
model.Rcirc = model.Rdoas;
%}


%% estimate model values
fprintf('Signal covariance...\n')
% convert bases to REAL vectorised versions with right dimensions 
r_circ = permute(cdmpcovr(model.Rcirc),[1 3 2]);     %[nChan^2, nDOA, nFreq]
r_iso = permute(cdmpcovr(model.Riso),[1 3 2]);      %[nChan^2, 1, nFreq]
r_white = permute(cdmpcovr(model.Rwhite),[1 3 2]);  %[nChan^2, 1, nFreq]
fixed_components = cat(2, r_iso, r_white);  % [nChan^2, 2, nFreq]

% smoothed covariance in real vectorised form
tinc = (pm.t(2)-pm.t(1));               % time increment between frames
ax = exp(-tinc./params.time_constant);  % smoothing factor
r_y = cdmpcovs(permute(Y,[3,1,2]),ax);  % [nChan^2,nFrame,nFreq]

% estimate the model parameters
fprintf('Estimate model parameters...\n')
% [rmodel,pow_component,imin,err] = cdmpfitras(...)
%           rmodel: reconstruced covariance matrix in real vectorised form
%                   this isn't so helpful for MVDR so easirer to
%                   reconstruct it outselves from pow_component
%    pow_component: Component weights: directional component first followed
%                   by fixed components [nb,nt,nf]
%             imin: index of directional component [nt,nf]
switch params.nFree
    case 0
        [rmodel,pow_component,imin,err]=cdmpfitras(r_y, fixed_components);
    case 1
        [rmodel,pow_component,imin,err]=cdmpfitras(r_y, fixed_components, ...
            r_circ);
    otherwise
        error('How did we get here? params.nFree was validated to be 0 or 1!')
end

% (thresholding) excluding target-nearby directions adjacent by setting their power component to zero
% pow_component(1,:,:)); size: 3 x nFrames x nFreq
% imin size: nFrames x nFreq
for iframe = 1:nFrames
    id_to_remove=find(ismember(imin(iframe,:),dirs_to_remove{iframe}));
    pow_component(1,iframe,id_to_remove)=0;
end

% for later clarity rearrange to be consistent with our normal arrangement
pow.iso = permute(pow_component(2,:,:), [3 2 1]);
pow.white = permute(pow_component(3,:,:), [3 2 1]);
pow.pw = permute(pow_component(1,:,:), [3 2 1]);
doa_index = imin.';

Rmodel = permute(cdmpcovc(rmodel),[1 2 4 3]); % [nChan nChan nFreq nFrame]

%% main processing
% preallocate variables to be populated
w = zeros(nFreq,nChan,nFrames); % beamformer weights for output channel
Z = zeros(nFreq,1,nFrames);      % output signals in stft domain

reverseStr = '';
for iframe = 1:nFrames
    msg = sprintf('Processing %2.2f s of %2.2f', pm.t(iframe), pm.t(end));
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));

    % update steering vector accounting for array rotation and source
    % movement
    ema.setPoseQuaternion(qr(:,iframe));
    [target_az, target_inc,~] = mycart2sph(target_doa_vec(iframe,:));
    
    d_t = fcn_20210223_01_get_transposed_rtf(ema,target_az,target_inc,stft_params.nfft);
    
    for ifreq = 1:nFreq
        % reconstruct the modeled covariance matrix
        R = pow.iso(ifreq,iframe) * model.Riso(:,:,ifreq) ...
            + pow.white(ifreq,iframe) * model.Rwhite(:,:,ifreq) ...
            + pow.pw(ifreq,iframe) * model.Rcirc(:,:,ifreq, doa_index(ifreq,iframe));
        
        if ~isequal(Rmodel(:,:,ifreq,iframe),R)
            diff_mat = Rmodel(:,:,ifreq,iframe) - R;
            if diff_mat>1e-18
                fprintf('freq %2.2f Hz, frame %2.2f s, max(abs(diff)) %2.2e\n',...
                    pm.f(ifreq), pm.t(iframe), max(abs(diff_mat(:))))
            end
        end
        
        
        % impose diagonal loading if necessary
        R = fcn_20180928_02_cap_condition_number_by_diagonal_loading(...
            R,params.max_condition_number);
     
        % compute beamformer weights
        w(ifreq,:,iframe) = conj(params.mvdr_fcn(R, d_t(:,ifreq), params.mvdr_debug_level).'); %
    end
end

% apply directly in frequency domain then obtain time domain signals
Z(:,1,:) = sum(w .* Y,2);
out = stft_v2('inv',Z,pm);

%% write out
fprintf('Saving results...\n')
elobes_writewav_wrapper(out,fs,out_wav_file_path);

if 1
    % check you're not overwriting data - if so append a _ to the name
    while exist(fullfile(saved_data_dir,sprintf('%s_filters.mat',params.filename_prefix)),'file')
        params.filename_prefix = strcat(params.filename_prefix,'_');
    end

    % now save everything
    save(fullfile(saved_data_dir,sprintf('%s_filters.mat',params.filename_prefix)),'w','pm','-v7.3');

    %{
    save(fullfile(saved_data_dir,sprintf('%s_model_estimates.mat',params.filename_prefix)),...
        'model','pow','doa_index','model_params','-v7.3');    
    %}
end
fprintf('end time: %s\n', datestr(now,'HH:MM:SS'));