%% setup & settings
setup;
% i/o root paths
input_root = './demo/';
output_root = './outputs/';

%ema_fcn = @RigidSphere_20210224_01_em32_locata_refchan_13; % for simulated RIR based on locata em32 class
ema_fcn = @RigidSphere_20171114_01_ideal_em32_refchan_13; % for measured RIR

% main params (max condition number & tau)
in_params.time_constant = 0.05;
in_params.max_condition_number = 1e2;

% subfolder - indicating the trial case
task_dir = 'case 1';
target_src ='src_1'; % target source

%% initialisation
fprintf('case: %s \ntarget: %s \n',task_dir,target_src);
in_wav_file_path = fullfile(input_root,task_dir,'mixed.wav');
out_dir = fullfile(output_root,task_dir,target_src);
if ~exist(out_dir) mkdir(out_dir); end
ht_file_path = fullfile(input_root,task_dir,'array_rot.wav');
listener_characteristics = [];

oracle_data.ema_fcn = ema_fcn;
oracle_data.target_position_csv_path = fullfile(input_root,task_dir,[target_src '.csv']);

saved_data_dir = out_dir;
temp_data_dir = out_dir;

%% Passthrough
%%{
out_wav_file_path = fullfile(out_dir,'enhanced_passthrough.wav');

bf_20210224_01_passthrough(...
    in_wav_file_path ...          % path to received signal
    ,out_wav_file_path ...        % path where enhanced signal(s) should be written
    ,ht_file_path ...             % ignored - path to head tracker signal
    ,listener_characteristics ... % ignored - path to listener characteristics profile
    ,in_params ...                % structure with config information for test bench
    ,oracle_data ...              % structure containing information and/or paths to information which would not be available in real-world application but may be useful in development
    ,saved_data_dir ...           % path to folder where any supplementary/intermediate results should be written
    ,temp_data_dir)
%}
%% Compact MVDR
%%{
prefix=['compact_' 'conNum_' num2str(in_params.max_condition_number) '_tau_' num2str(in_params.time_constant)];
out_wav_file_path = fullfile(out_dir,['enhanced_' prefix '.wav']);
in_params.filename_prefix = prefix;
mvdr_20210601_compact(...
    in_wav_file_path ...          % path to received signal
    ,out_wav_file_path ...        % path where enhanced signal(s) should be written
    ,ht_file_path ...             % ignored - path to head tracker signal
    ,listener_characteristics ... % ignored - path to listener characteristics profile
    ,in_params ...                % structure with config information for test bench
    ,oracle_data ...              % structure containing information and/or paths to information which would not be available in real-world application but may be useful in development
    ,saved_data_dir ...           % path to folder where any supplementary/intermediate results should be written
    ,temp_data_dir)
%}
%% Isotropic MVDR 
%%{
prefix=['isotropic_' 'conNum_' num2str(in_params.max_condition_number) '_tau_' num2str(in_params.time_constant)];
out_wav_file_path = fullfile(out_dir,['enhanced_' prefix '.wav']);
in_params.filename_prefix = prefix;
mvdr_20210601_isotropic(...
    in_wav_file_path ...          % path to received signal
    ,out_wav_file_path ...        % path where enhanced signal(s) should be written
    ,ht_file_path ...             % ignored - path to head tracker signal
    ,listener_characteristics ... % ignored - path to listener characteristics profile
    ,in_params ...                % structure with config information for test bench
    ,oracle_data ...              % structure containing information and/or paths to information which would not be available in real-world application but may be useful in development
    ,saved_data_dir ...           % path to folder where any supplementary/intermediate results should be written
    ,temp_data_dir)
%}
%% MPDR
%%{
prefix=['mpdr_' 'conNum_' num2str(in_params.max_condition_number) '_tau_' num2str(in_params.time_constant)];
out_wav_file_path = fullfile(out_dir,['enhanced_' prefix '.wav']);
in_params.filename_prefix = prefix;
mvdr_20210601_mpdr(...
    in_wav_file_path ...          % path to received signal
    ,out_wav_file_path ...        % path where enhanced signal(s) should be written
    ,ht_file_path ...             % ignored - path to head tracker signal
    ,listener_characteristics ... % ignored - path to listener characteristics profile
    ,in_params ...                % structure with config information for test bench
    ,oracle_data ...              % structure containing information and/or paths to information which would not be available in real-world application but may be useful in development
    ,saved_data_dir ...           % path to folder where any supplementary/intermediate results should be written
    ,temp_data_dir)
%}
%% Oracle MVDR
%%{
prefix=['oracle_' 'conNum_' num2str(in_params.max_condition_number) '_tau_' num2str(in_params.time_constant)];
out_wav_file_path = fullfile(out_dir,['enhanced_' prefix '.wav']);
in_params.filename_prefix = prefix;
mvdr_20210601_oracle(...
    in_wav_file_path ...          % path to received signal
    ,out_wav_file_path ...        % path where enhanced signal(s) should be written
    ,ht_file_path ...             % ignored - path to head tracker signal
    ,listener_characteristics ... % ignored - path to listener characteristics profile
    ,in_params ...                % structure with config information for test bench
    ,oracle_data ...              % structure containing information and/or paths to information which would not be available in real-world application but may be useful in development
    ,saved_data_dir ...           % path to folder where any supplementary/intermediate results should be written
    ,temp_data_dir)
%}