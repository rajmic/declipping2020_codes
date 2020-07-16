%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%                                %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%      DECLIPPING MAIN FILE      %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%             Janssen            %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%                                %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Ondøej Mokrý and Pavel Záviška, Brno University of Technology, 2020
%
% using toolbox LTFAT
ltfatstart

addpath('Tools')
addpath('Sounds')
addpath('Evaluation_algorithms')

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
fprintf('Setting up the parameters\n')

% input SDR of the clipped signal
param.inputSDR = 7;     % set the input SDR value

% window parameters
param.w = 8192;       % window length
param.a = param.w/4;  % window shift
param.wtype = 'hann'; % options available on: http://ltfat.github.io/doc/sigproc/firwin.html

% Janssen param
paramsolver.NIt = 3;    % number of inner iterations
paramsolver.p = 512;    % order of the model

paramsolver.verbose = true; % display parameter

%% clipping 
fprintf('Generating clipped signal\n')
[data_clipped, param.masks, param.theta, trueSDR, percentage] = clip_sdr(data, param.inputSDR); % clipping of original signal
fprintf('Clipping threshold %4.3f, true SDR value is %4.2f dB and %4.2f%% samples has been clipped \n', param.theta, trueSDR, percentage)

% mask of the reliable samples
param.mask = param.masks.Mr;

        
%% Main algorithm    

tic
data_rec = janssen_2(data_clipped, param, paramsolver);
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
