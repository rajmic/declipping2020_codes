%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%                                %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%      DECLIPPING MAIN FILE      %%%%%%%%%%%%%%%%%%%%%%%     
%%%%%%%%%%%%%%%%%%%% CONSISTENT DICTIONARY LEARNING %%%%%%%%%%%%%%%%%%%%%%% 
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

%% Input file settings   
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


%% Parameters

inputSDR = 7;       % set the input SDR value

fprintf('Setting up the frame parameters\n')

param.N = 1024;     % window length

param.hop = 0.25*param.N;   % window shift
param.redundancyFactor = 2; % Dictionary redundancy
param.M = param.N * param.redundancyFactor; % number of atoms
param.target_error = 0;
param.wa = @wRect;
param.ws = param.wa;

N = param.N;
M = param.M;

D_DCT = DCT_Dictionary(param);


%% Process input signal

% crop signal
L = length(frames2signal(signal2frames(x,param),param));
x = x(1:L);

% clip signal
fprintf('Generating clipped signal\n');
[y, masks, clipping_threshold, trueSDR, percentage] = clip_sdr(x, inputSDR);
fprintf('Clipping threshold %4.3f, true SDR value is %4.2f dB and %4.2f%% samples has been clipped \n', clipping_threshold, trueSDR, percentage)


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
paramIHT.A_init = zeros(M,Nframes);
paramIHT.save_log = 0;

%t = cputime;
fprintf('Starting the Consistent Iterative Hard Thresholding alorithm.\n')
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

fprintf('Starting the Consistent Dictionary Learning algorithm.\n')
[D,A, alg_log] = consDL(Y,paramDL,paramSC,paramDictUpdate, distortion_param);
time = toc;

% Reconstruct signal:
X_est = D*A;       
x_est = frames2signal(X_est,param);


%% Time & SDR evaluation

% time
fprintf('Result obtained in %4.3f seconds.\n', time);

% SDR
sdr_clip = sdr(x, y);
sdr_rec = sdr(x, x_est);
fprintf('SDR of the clipped signal is %4.3f dB.\n', sdr_clip);
fprintf('SDR of the reconstructed signal is %4.3f dB.\n', sdr_rec);
fprintf('SDR improvement is %4.3f dB.\n', sdr_rec-sdr_clip);
