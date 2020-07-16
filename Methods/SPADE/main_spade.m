%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%                                %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%      DECLIPPING MAIN FILE      %%%%%%%%%%%%%%%%%%%%%%%     
%%%%%%%%%%%%%%%%%%%%              SPADE             %%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%                                %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Pavel Záviška, Brno University of Technology, 2020
%
% using toolbox LTFAT
ltfatstart

addpath('Tools')
addpath('Sounds')
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
[data, fs] = audioread(['Sounds/' audio_file '.wav']);

% signal length
param.Ls = length(data);


%% settings
fprintf('Setting up the frame parameters\n')

% input SDR of the clipped signal
param.inputSDR = 7;     % set the input SDR value

% window parameters
param.w = 8192;       % window length
param.a = param.w/4;  % window shift
param.wtype = 'hann'; % options available on: http://ltfat.github.io/doc/sigproc/firwin.html

% DFT parameters
param.F = frame('dft');
param.F.redundancy = 2;  % non-native, our own additional parameter
param.F.frana = @(insig)dft([insig; zeros(length(insig)*(param.F.redundancy-1),1)]);
param.F.frsyn = @(insig)postpad(idft(insig),length(insig)/param.F.redundancy);

% paramsolver parameters
paramsolver.s = 1;   % increment of k
paramsolver.r = 2;   % every r-th iteration increment k by s   
paramsolver.epsilon = 0.1;  % stopping criterion of termination function
paramsolver.maxit = ceil(floor(param.w*param.F.redundancy/2+1)*paramsolver.r/paramsolver.s); % maximum number of iterations

paramsolver.store_sdr = 1; 
paramsolver.store_obj = 1;

paramsolver.verbose = true; % 0 - no listing, 1 - inform about processed blocks

% algorithm
param.algorithm = 'aspade'; % algorithm to compute declipping, options: 'aspade', 'sspade_new'
fprintf('Setting up algorithm: %s\n', param.algorithm)


%% clipping 
fprintf('Generating clipped signal\n')
[data_clipped, param.masks, param.theta, trueSDR, percentage] = clip_sdr(data, param.inputSDR); % clipping of original signal
fprintf('Clipping threshold %4.3f, true SDR value is %4.2f dB and %4.2f%% samples has been clipped \n', param.theta, trueSDR, percentage)


%% Main algorithm
fprintf('Starting the %s algorithm\n', param.algorithm)

tic; 

[data_rec, sdr_iter, obj_iter] = spade_segmentation(data_clipped, param, paramsolver, data);

time = toc;


%% Time & SDR evaluation

% time
fprintf('Result obtained in %4.3f seconds.\n', time);

% SDR
sdr_clip = sdr(data, data_clipped);
sdr_rec = sdr(data, data_rec);
fprintf('SDR of the clipped signal is %4.3f dB.\n', sdr_clip);
fprintf('SDR of the reconstructed signal is %4.3f dB.\n', sdr_rec);
fprintf('SDR improvement is %4.3f dB.\n', sdr_rec-sdr_clip);

