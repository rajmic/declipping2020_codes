%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%                                %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%      DECLIPPING MAIN FILE      %%%%%%%%%%%%%%%%%%%%%%%     
%%%%%%%%%%%%%%%%%%%%              C-OMP             %%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%                                %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% The implementation of the Constrained OMP Declipper was kindly provided
% by Valentin Emiya.
%
% Main file to fit the common interface of the decliping toolbox created by 
% Pavel Záviška, Brno University of Technology, 2020

% using CVX toolbox

addpath('Tools')
addpath('Sounds')
addpath(genpath('Methods/Constrained OMP'))
addpath(genpath('Evaluation_algorithms'))

close all
clear variables


%% input file settings   
% load audio-file
audio_file = 'a08_violin';    %  'a08_violin'
                              %  'a16_clarinet'
                              %  'a18_bassoon'
                              %  'a25_harp'	
                              %  'a35_glockenspiel'
                              %  'a41_celesta'
                              %  'a42_accordion'
                              %  'a58_guitar_sarasate'
                              %  'a60_piano_schubert'
                              %  'a66_wind_ensemble_stravinsky'
                                
fprintf(['Loading audio ''', audio_file , '.wav''\n']);
[x, fs] = audioread(['Sounds/' audio_file '.wav']);


%% settings

% input SDR of the clipped signal
inputSDR = 7;     % set the input SDR value

fprintf('Setting up the frame parameters\n')

param.N = 1024; % frame length
param.inpaintFrame = @inpaintFrame_consOMP_Gabor; % solver function
param.OMPerr = 0.001;
param.sparsityDegree = 64;
param.D_fun = @Gabor_Dictionary; % Dictionary (function handle)
param.OLA_frameOverlapFactor = 4;
param.redundancyFactor = 2; % Dictionary redundancy
param.wd = @wRect; % Weighting window for dictionary atoms
param.wa = @wRect; % Analysis window
param.OLA_ws = @wSine; % Synthesis window
param.SKIP_CLEAN_FRAMES = true; % do not process frames where there is no missing samples
param.MULTITHREAD_FRAME_PROCESSING = false; % not implemented yet
param.VERBOSE = true;

%% Process input signal

fprintf('Generating clipped signal\n');
[y, masks, clipping_threshold, trueSDR, percentage] = clip_sdr(x, inputSDR);
fprintf('Clipping threshold %4.3f, true SDR value is %4.2f dB and %4.2f%% samples has been clipped \n', clipping_threshold, trueSDR, percentage)

problemData.x = y;
problemData.IMiss = ~masks.Mr;


%% Declip using cOMP

fprintf('Starting the Constrained OMP algorithm.\n')
tic
[x_est, ~] = inpaintSignal_IndependentProcessingOfFrames(problemData, param);
time = toc;


%% Time & SDR evaluation

% time
fprintf('Result obtained in %4.3f seconds.\n', time);

% crop the result
L = min(length(x), length(x_est));
x = x(1:L);
y = y(1:L);
x_est = x_est(1:L);

% SDR
sdr_clip = sdr(x, y);
sdr_rec = sdr(x, x_est);
fprintf('SDR of the clipped signal is %4.3f dB.\n', sdr_clip);
fprintf('SDR of the reconstructed signal is %4.3f dB.\n', sdr_rec);
fprintf('SDR improvement is %4.3f dB.\n', sdr_rec-sdr_clip);
