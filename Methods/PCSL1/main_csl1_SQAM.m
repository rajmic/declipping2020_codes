%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%                                %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%      DECLIPPING MAIN FILE      %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%           (P(W))CSL1           %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%          (all sounds)          %%%%%%%%%%%%%%%%%%%%%%%
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

%% Database
sounds = {'a08_violin', 'a16_clarinet', 'a18_bassoon', 'a25_harp', 'a35_glockenspiel', 'a41_celesta', ...
    'a42_accordion', 'a58_guitar_sarasate', 'a60_piano_schubert', 'a66_wind_ensemble_stravinsky'};

algorithms = {'csl1', 'pcsl1', 'pwcsl1'};
input_SDRs = [1, 3, 5, 7, 10, 15, 20];

% load audio samples
load('Sounds\Sound_database.mat')

% initialization of matrices for SDR and dSDR values, and computational time
SDR = NaN(length(sounds), length(input_SDRs));
dSDR = NaN(length(sounds), length(input_SDRs));
TIME = NaN(length(sounds), length(input_SDRs));

% initialization of counter
cnt = 0;

for sound = 1:length(sounds)
    for clip_idx = 1:length(input_SDRs)
        for algorithm = 1
            
            cnt = cnt+1;
            %% input file settings
            % load audio-file
            eval(['data = ', sounds{sound}, ';']);
            
            % signal length
            param.Ls = length(data);
            
            % sampling frequency
            param.fs = 44100; % Hz

            
            %% settings
            param.inputSDR = input_SDRs(clip_idx);    % set (symetric) clipping threshold
            
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
            param.algorithm = algorithms{algorithm}; % algorithm to compute declipping, options: 'aspade', 'sspade', 'sspade_new'

            % paramsolver parameters
            paramsolver = settings_csl1(param.algorithm);
            
            paramsolver.store_dsdr = 0; 
            paramsolver.store_obj = 0;
            
            paramsolver.verbose = false;
                                              
            %% clipping
            [data_clipped, param.masks, param.theta, trueSDR, percentage] = clip_sdr(data, param.inputSDR); % clipping of original signal
            
            
            %% Main algorithm
            tic;
            
            [data_rec, sdr_iter, obj_iter] = csl1_segmentation(data_clipped, param, paramsolver, data);
            
            time = toc;
            
            
            %% Time & SDR evaluation
            
            % Time
            TIME(sound, clip_idx) = time;

            % Save restored file
            if input_SDRs(clip_idx) < 10
                eval([sounds{sound} '_rec_' algorithms{algorithm} '_0' num2str(input_SDRs(clip_idx)) ' = data_rec;']);
            else
                eval([sounds{sound} '_rec_' algorithms{algorithm} '_' num2str(input_SDRs(clip_idx)) ' = data_rec;']);
            end

            % SDR
            sdr_clip = sdr(data, data_clipped);
            sdr_rec = sdr(data, data_rec);

            SDR(sound, clip_idx) = sdr_rec;
            dSDR(sound, clip_idx) = sdr_rec - sdr_clip;


            disp(['Done: ' num2str(cnt), '/70'  ]);
        end
        %save('Results/CSL1_Sounds.mat');
    end
end

