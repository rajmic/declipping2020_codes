%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%                                %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%      DECLIPPING MAIN FILE      %%%%%%%%%%%%%%%%%%%%%%%     
%%%%%%%%%%%%%%%%%%%% CONSISTENT DICTIONARY LEARNING %%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%          (all sounds)          %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%                                %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% The implementation of the Constrained OMP Declipper was kindly provided
% by Lucas Rencker.
%
% Main file to fit the common interface of the decliping toolbox created by 
% Pavel Záviška, Brno University of Technology, 2020


addpath('Tools')
addpath('Sounds')
addpath(genpath('Methods/Consistent Dictionary Learning'))
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

% Parameters
param.N = 1024;

param.hop = 0.25*param.N;
param.redundancyFactor = 2; % Dictionary redundancy
param.M = param.N * param.redundancyFactor; % number of atoms
param.target_error = 0;
param.wa = @wRect;
param.ws = param.wa;

D_DCT = DCT_Dictionary(param);

% initialization of counter
cnt = 0;



for sound = 1:length(sounds)
    for clip_idx = 1:length(input_SDRs)
    
        cnt = cnt+1;
        %% Input file
        % load audio-file              
        eval(['x = ', sounds{sound}, ';']);

        %% Process input singal
        
        % crop signal
        L = length(frames2signal(signal2frames(x,param),param));
        x = x(1:L);
        
        % clip signal
        inputSDR = input_SDRs(clip_idx);    % set (symetric) clipping threshold
        [y, masks, clipping_threshold, trueSDR, percentage] = clip_sdr(x, inputSDR);

        y = y(1:L);

        % decompose signal into frames
        Y = signal2frames(y,param);

        Nframes = size(Y,2);

        reliable_samples_mat = binary_vec2mat(masks.Mr,param);

        distortion_param.name = 'clipping';
        distortion_param.reliable_samples_mat = reliable_samples_mat;
        distortion_param.positive_clipped = ~reliable_samples_mat & (Y>=0); 
        distortion_param.negative_clipped = ~reliable_samples_mat & (Y<=0); 

        
        %% Iterative Hard Thresholding (IHT) with adaptive sparsity for initialization:

        paramIHT.Nit = 2000; % max number of iterations
        paramIHT.rate = 10; % increases sparsity level every 10 iterations
        paramIHT.target_cost = sum(sum(reliable_samples_mat)) * 1e-6; 
        paramIHT.A_init = zeros(param.M,Nframes);
        paramIHT.save_log = 0;

        %t = cputime;
        tic
        [A_IHT,alg_log_IHT] = consIHT(Y,D_DCT,paramIHT,distortion_param);

        K_IHT = floor(nnz(A_IHT)/Nframes);  % keep average sparsity level in memory


        %% Dictionary learning:

        % DL parameters:
        paramDL.Nit = 20;
        paramDL.A_init = A_IHT;  % initialize with coefficients from IHT
        paramDL.D_init = D_DCT;
        paramDL.warm_start = 1;
        paramDL.loud = 0; % print results

        % sparse coding parameters:
        paramSC.alg = @consIHT;
        paramSC.Nit = 20;
        paramSC.rate = 0;
        paramSC.K = K_IHT;
        paramSC.accuracy = 0; 
        paramSC.save_log = 0;

        % dict update parameters:
        paramDictUpdate.Nit = 20;
        paramDictUpdate.accelerate = 1;
        paramDictUpdate.save_log = 0;

        [D,A, alg_log] = consDL(Y,paramDL,paramSC,paramDictUpdate, distortion_param);
        time = toc;

        % Reconstruct signal:
        X_est = D*A;       
        x_est = frames2signal(X_est,param);


        %% Time & SDR evaluation

        % Time
        TIME(sound, clip_idx) = time;

        % Save restored file
        if input_SDRs(clip_idx) < 10
            eval([sounds{sound} '_rec_DL_0' num2str(input_SDRs(clip_idx)) ' = x_est;']);
        else
            eval([sounds{sound} '_rec_DL_' num2str(input_SDRs(clip_idx)) ' = x_est;']);
        end

        % SDR
        sdr_clip = sdr(x, y);
        sdr_rec = sdr(x, x_est);

        SDR(sound, clip_idx) = sdr_rec;
        dSDR(sound, clip_idx) = sdr_rec - sdr_clip;


        disp(['Done: ' num2str(cnt), '/70'  ]);

        %save('Results/DL_Sounds.mat');
    end
end

