%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%                                %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%      DECLIPPING MAIN FILE      %%%%%%%%%%%%%%%%%%%%%%%     
%%%%%%%%%%%%%%%%%%%%           (P(W))CSL1           %%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%                                %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Implementation of the declipping methods according to B. Defraene et al: 
% Declipping of audio signals using perceptual compressed sensing
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

% sampling frequency
param.fs = fs;


%% settings
fprintf('Setting up the frame parameters\n')

% input SDR of the clipped signal
param.inputSDR = 20;     % set the input SDR value

% window parameters
param.w = 8192;       % window length
param.a = param.w/4;  % window shift
param.wtype = 'hann'; % options available on: http://ltfat.github.io/doc/sigproc/firwin.html

% DFT parameters
param.F = frame('dft');
param.F.redundancy = 2;  % redundancy of the DFT transform
param.F.frana = @(insig)dft([insig; zeros(length(insig)*(param.F.redundancy-1),1)]);
param.F.frsyn = @(insig)postpad(idft(insig),length(insig)/param.F.redundancy);

% algorithm
param.algorithm = 'csl1'; % algorithm to compute declipping, options: 'csl1', 'pcsl1', 'pwcsl1'
fprintf('Setting up algorithm: %s\n', param.algorithm)

% paramsolver parameters
paramsolver = settings_csl1(param.algorithm);

paramsolver.store_dsdr = 1;
paramsolver.store_obj = 1;

paramsolver.verbose = true;

%% clipping 
fprintf('Generating clipped signal\n')
[data_clipped, param.masks, param.theta, trueSDR, percentage] = clip_sdr(data, param.inputSDR); % clipping of original signal
fprintf('Clipping threshold %4.3f, true SDR value is %4.2f dB and %4.2f%% samples has been clipped \n', param.theta, trueSDR, percentage)


%% Main algorithm
fprintf('Starting the %s algorithm\n', param.algorithm)

tic; 

[data_rec, sdr_iter, obj_iter] = csl1_segmentation(data_clipped, param, paramsolver, data);

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

