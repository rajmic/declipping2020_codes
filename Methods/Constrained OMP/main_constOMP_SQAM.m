%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%                                %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%      DECLIPPING MAIN FILE      %%%%%%%%%%%%%%%%%%%%%%%     
%%%%%%%%%%%%%%%%%%%%              C-OMP             %%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%           (all sounds)         %%%%%%%%%%%%%%%%%%%%%%%
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

% Database
sounds = {'a08_violin', 'a16_clarinet', 'a18_bassoon', 'a25_harp', 'a35_glockenspiel', 'a41_celesta', ...
    'a42_accordion', 'a58_guitar_sarasate', 'a60_piano_schubert', 'a66_wind_ensemble_stravinsky'};

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

        cnt = cnt+1;
        %% Input file
        % load audio-file              
        eval(['x = ', sounds{sound}, ';']); 
        
        %% settings

        % input SDR of the clipped signal
        inputSDR = input_SDRs(clip_idx);     % set the input SDR value
        
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
        param.VERBOSE = false;
        
        %% clipping
        [y, masks, clipping_threshold, ~, ~] = clip_sdr(x, inputSDR); % clipping of original signal
        
        problemData.x = y;
        problemData.IMiss = ~masks.Mr;

        %% Declip using cOMP
        tic
        [x_est, ~] = inpaintSignal_IndependentProcessingOfFrames(problemData, param);
        time = toc;
        
        %% Time & SDR evaluation

        % Time
        TIME(sound, clip_idx) = time;

        % crop the result
        L = min(length(x), length(x_est));
        x = x(1:L);
        y = y(1:L);
        x_est = x_est(1:L);
        
        % Save restored file
        if input_SDRs(clip_idx) < 10
            eval([sounds{sound} '_rec_consOMP_0' num2str(input_SDRs(clip_idx)) ' = x_est;']);
        else
            eval([sounds{sound} '_rec_consOMP_' num2str(input_SDRs(clip_idx)) ' = x_est;']);
        end

        % SDR
        sdr_clip = sdr(x, y);
        sdr_rec = sdr(x, x_est);

        SDR(sound, clip_idx) = sdr_rec;
        dSDR(sound, clip_idx) = sdr_rec - sdr_clip;


        disp(['Done: ' num2str(cnt), '/70'  ]);
        
        %save('Results/consOMP_Sounds.mat')
    end
end
