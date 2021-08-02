function[] = mvdr_20210601_isotropic(...
    in_wav_file_path ...          % path to received signal
    ,out_wav_file_path ...        % path where enhanced signal(s) should be written
    ,ht_file_path ...             % ignored - path to head tracker signal
    ,listener_characteristics ... % ignored - path to listener characteristics profile
    ,in_params ...                % structure with config information for test bench
    ,oracle_data ...              % structure containing information and/or paths to information which would not be available in real-world application but may be useful in development
    ,saved_data_dir ...           % path to folder where any supplementary/intermediate results should be written
    ,temp_data_dir)            % path to folder where temporary files should be written - normally empty but may be useful for debugging

% mvdr_20210223_01_moving_target_isotropic_field implements an mvdr
% beamformer where the PSD/covariance matrix is modelled as a 
% circularly isotropic noise.
% The look direction is steered towards the known source direction in response to
% array rotation. Unlike other beamformers, this one works directly in frequency
% domain ignoring circular convolution issues for simplicity
%
% Usage:
%  mvdr_20210223_01_moving_target_isotropic_field(...
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
% Alastair Moore, February 2021
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
params.circ_az = deg2rad((0:5:355).' );         % assuming circular diffuse field

%% update with input parameters
if ~isempty(in_params)
    params = override_valid_fields(params,in_params);
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


%% precompute model covariance bases
fprintf('Computing diffuse covariance matrices...\n')
model_params.fs = params.fs;
model_params.stft_params = stft_params;
model_params.ema = ema;


% cylindrical isotropic
%{
model_params.doa_circ.az = params.circ_az;
model_params.doa_circ.inc = deg2rad( 90 * ones(size(model_params.doa_circ.az)));
model = fcn_20210215_ahm_02_get_compact_model_bases(model_params);
%}


% spherical isotropic (azi-inc grid w/ quadrature weighting)
%%{
model_params.iso_grid.azi_step=10;
model_params.iso_grid.inc_step=10;
model_params.doas=[0 0]; % not used in isotorpic
model = fcn_20210506_get_compact_model_bases(model_params);
%}

%% diffuse-field assumption
Rest = repmat(model.Riso,[1, 1, 1, nFrames]);

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
    % original
    ema.setPoseQuaternion(qr(:,iframe));
    [target_az, target_inc,~] = mycart2sph(target_doa_vec(iframe,:));
    
    d_t = fcn_20210223_01_get_transposed_rtf(ema,target_az,target_inc,stft_params.nfft);

    for ifreq = 1:nFreq    
        % impose diagonal loading if necessary
        R = fcn_20180928_02_cap_condition_number_by_diagonal_loading(...
            Rest(:,:,ifreq,iframe), params.max_condition_number);
     
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
       
end
fprintf('end time: %s\n', datestr(now,'HH:MM:SS'));