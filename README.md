# sap-papers-2021-elospheres-compact-model
SAP Group - Papers - 2021 - Compact Model

*** INPUTS ***
1. have a folder per test case scenario
2. make sure the scenario contains the following mandatory files:
  2a. src_xx.csv      containing the cartesian position of a single source over time. The content is a Nx4 matrix, where each row is a time instant and columns are time (second), x, y , z position in m. The identifier 'xx' in the file name should be an integer number associated with source id. e.g. src_1, src_2 etc.
  2b. src_xx.wav      the anechoic source signal at the microphone array. The content is a Mx32 matrix, where each column is time series signal for a microphone in the array. The identifier 'xx' in the file name should be an integer number associated with source id. e.g. src_1, src_2 etc.
  2c. mixed.wav       the microphone signal. The content is a Mx32 matrix, where each column is time series signal for a microphone in the array.
  2d. array_rot.wav   the orientation of the array over time. The content is a Lx3 matrix, where columns are respectively roll, pitch and yaw normalised to [-1 1] range (-1: -pi or 180, 0: 0, +1: +pi or 180). 
3. folder demo includes an example test case materials. It is based on METU SPARG measured Room Impulse Response for em2 Eigenmike with the following setup:
  3a. Room: 8.3×6.5×2.9 m with T60=1.12 s
  3b. Array: em32 Eigenmike located at ( 4.65 , 3.25 , 1.5 )m
  3c. Sources: two stationary sources at 0 and 27 deg azimuth, both with 90 deg (horizontal plane) inclination. approx. 1m source-array distance. C50=5.4, 4.5 dB respectively.
  3d. Babble: 8 sources with close-to-uniform circular distribution on the horizontal plane with azimuths=18~342 deg and C50=2.3~3.3 dB. 
  3e. signal-to-sensor noise=40 dB. signal-to-interference=-10 dB. signal-to-babble-noise=10 dB.
  
***  SETUP ***
1. make sure setup.m is updated with the correct path to the external modules/library (available via Imperial College GitHub)
2. make sure run_all_methods is updated with the correct root path to the input parent folder. 
3. make sure 'task_dir' parameter is the subfolder name for the case of interest in the input's parent folder and 'target_src' parameter is the name of the target source file (e.g. src_1) as specified in the input folder. 

***  RUN ***
1. (optional) you can change the main settings 'in_params.time_constant' (temporal reaction time) and 'in_params.max_condition_number' (regularisation parameter for covariance matrix)
2. run 'run_all_methods'

***  OUTPUTS ***
The output of the beamformers at the reference microphone are stored in wav files. 
The beamformer weights in the STFT domain along with extra info are stored in mat files.
