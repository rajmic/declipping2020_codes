%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%                                %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%      DECLIPPING MAIN FILE      %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%             Janssen            %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%           (all sounds)         %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%                                %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Ondøej Mokrý and Pavel Záviška, Brno University of Technology, 2020
%
% using toolbox LTFAT
ltfatstart

close all
clear variables 

%% Order of the model
p = 512;

%% Database
sounds = {'a08_violin', 'a16_clarinet', 'a18_bassoon', 'a25_harp', 'a35_glockenspiel', 'a41_celesta', ...
    'a42_accordion', 'a58_guitar_sarasate', 'a60_piano_schubert', 'a66_wind_ensemble_stravinsky'};

input_SDRs = [1, 3, 5, 7, 10, 15, 20];

% load audio samples
load('Sounds\Sound_database.mat')

% initialization of matrices for SDR and dSDR values, and computational time
TIME = NaN(length(sounds), length(input_SDRs));

% initialization of counter
cnt = 0;


for sound = 1:length(sounds)
    for clip_idx = 1:length(input_SDRs)
        
        cnt = cnt+1;
        
        %% input file 
        % load audio-file              
        eval(['data = ', sounds{sound}, ';']);
    
        % signal length
        param.Ls = length(data);

        
        %% settings
        param.inputSDR = input_SDRs(clip_idx);    % set (symetric) clipping threshold

        % window parameters
        param.w = 8192;       % window length
        param.a = param.w/4;  % window shift
        param.wtype = 'hann'; % options available on: http://ltfat.github.io/doc/sigproc/firwin.html
        
        % Janssen param
        paramsolver.NIt = 3; %[ 1 2 3 5 8 10 ];
        paramsolver.p = p;
        paramsolver.verbose = false;
        
        %% clipping
        [data_clipped, param.masks, param.theta, ~, ~] = clip_sdr(data, param.inputSDR); % clipping of original signal
        
        % mask of the reliable samples
        param.mask = param.masks.Mr;
        
        
        %% Main algorithm    
        
        tic
        data_rec_all = janssen_2(data_clipped, param, paramsolver);
        time = toc;
        
        %% Time & SDR evaluation
        
        TIME(sound, clip_idx) = time;
        
        nitcount = 0;
        for nit = paramsolver.NIt
        
            nitcount = nitcount + 1;
            data_rec = data_rec_all(:,nitcount);
            
            if input_SDRs(clip_idx) < 10
                eval([sounds{sound} '_rec_' num2str(p) '_' num2str(nit) '_0' num2str(input_SDRs(clip_idx)) ' = data_rec;']);
            else
                eval([sounds{sound} '_rec_' num2str(p) '_' num2str(nit) '_' num2str(input_SDRs(clip_idx)) ' = data_rec;']);
            end

            % SDR
            sdr_clip = sdr(data, data_clipped);
            sdr_rec  = sdr(data, data_rec);          
            
            eval(['SDR_janssen_'  num2str(p) '_' num2str(nit) '(sound, clip_idx) = sdr_rec;']);
            eval(['dSDR_janssen_' num2str(p) '_' num2str(nit) '(sound, clip_idx) = sdr_rec - sdr_clip;']);

            %save('Janssen_512_Sounds.mat');
        
        end
        
        disp(['Done: ' num2str(cnt), '/70'  ]);
    end
end