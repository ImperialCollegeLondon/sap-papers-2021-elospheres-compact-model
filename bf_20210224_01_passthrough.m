function[] = bf_20210224_01_passthrough(...
    in_wav_file_path ...          % path to received signal
    ,out_wav_file_path ...        % path where enhanced signal(s) should be written
    ,ht_file_path ...             % ignored - path to head tracker signal
    ,listener_characteristics ... % ignored - path to listener characteristics profile
    ,in_params ...                % structure with config information for test bench
    ,oracle_data ...              % structure containing information and/or paths to information which would not be available in real-world application but may be useful in development
    ,saved_data_dir ...           % path to folder where any supplementary/intermediate results should be written
    ,temp_data_dir ...            % path to folder where temporary files should be written - normally empty but may be useful for debugging
    )
% just copies reference channel to output. resamples if necessary
%
% Usage:
%  bf_20210224_01_passthrough(...
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
% Alastair Moore, February 2021 2019

%% validate the input
% not really required in this simple function but will use it as a template

% in_wav_file_path
if ~exist(in_wav_file_path,'file'), error('Couldn''t find %s',in_wav_file_path), end

% final_data_dir
final_data_dir = fileparts(out_wav_file_path);
check_output_dir_exists(final_data_dir)


%% defualt parameters
params.fs = [];

% update with input parameters
if ~isempty(in_params)
    params = override_valid_fields(params,in_params);
end

%% read in data and do pre-processing
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

% instantiate the microphone array
% - only actually need the ref channel so don't need to call prepareData
ema = oracle_data.ema_fcn();

% copy reference channel to output
out = y(:,ema.refChan);

% write out
elobes_writewav_wrapper(out,params.fs,out_wav_file_path);
