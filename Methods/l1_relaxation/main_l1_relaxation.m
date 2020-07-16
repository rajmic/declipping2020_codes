%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%                                %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%      DECLIPPING MAIN FILE      %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%          l1-relaxation         %%%%%%%%%%%%%%%%%%%%%%%
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
fprintf('Setting up parameters\n')

% input SDR of the clipped signal
param.inputSDR = 7;     % set the input SDR value

% DGT parameters
param.wtype = 'hann';  % window type
param.w = 8192;        % window length 
param.a = param.w / 4; % window shift
param.M = 2*8192;        % number of frequency channels     

% Algorithm
param.algorithm = 'DR'; % algorithm to compute declipping, options: Douglas-Rachford (synthesis model): 'DR', Chamboll-Pock (analysis model): 'CP'
fprintf('Setting up algorithm: %s\n', param.algorithm)

param.replaceReliable = 0;  % after finishing all iterations, replace reliable samples in declipped signal with original signal

% reweighting (Rl1CC)
param.reweighting = 0; % enables reweighting of the algorithms (Rl1CC algorithms)

% weighting
param.weighting = 0; % enables weighting of the coefficients during the restoration
param.weighting_type = 'par'; % options: 'ath', 'mpeg', 'par'
param.weights_assignment = 1; % options: 1, 2, 3
param.psychoacoustic_model_source = 'clipped'; % options: 'clipped', 'reconstructed', 'original'

% paramsolver parameters
paramsolver = settings_l1_relaxation(param.algorithm, param.reweighting);
paramsolver.verbose = 1;      % display parameter
paramsolver.comp_dsdr = 1;    % compute and store dSDR during iterations
paramsolver.comp_obj = 1;     % compute and store objective function values during iterations


%% clipping
fprintf('Generating clipped signal\n')
[data_clipped, param.masks, param.theta, trueSDR, percentage] = clip_sdr(data, param.inputSDR); % clipping of original signal
fprintf('Clipping threshold %4.3f, true SDR value is %4.2f dB and %4.2f%% samples has been clipped \n', param.theta, trueSDR, percentage)

%% construction of frame
fprintf('Creating the frame\n')

param.F = frame('dgtreal', {param.wtype, param.w}, param.a, param.M); % creating requested frame

param.F = frametight(param.F); % creating Parseval tight frame
param.F = frameaccel(param.F, param.Ls);  % precomputation for a fixed signal length


% construction of the weights
if param.weighting
    fprintf('Constructing the weights\n')
    switch param.weighting_type
        case {'ATH', 'ath'}
            param.weights = weights_ath(0, fs/2, floor(param.M/2)+1, param.F, param.Ls, param.weights_assignment);
        case {'MPEG', 'mpeg'}
            switch param.psychoacoustic_model_source
                case {'clipped', 'clip'}
                    model_source = data_clipped;
                case {'original', 'orig'}
                    model_source = data;
                case {'reconstructed', 'rec'}
                    if exist('data_rec', 'var')
                        model_source = data_rec;
                    else
                        error('The recostructed signal is not available. Please run the algorithm without weighting first to obtain the reconstructed signal.')
                    end
                otherwise
                    error('Invalid source for the psychoacoustic model!')
            end
            param.weights = weights_mpeg(param.w, param.M, param.a, param.F, model_source, fs, param.weights_assignment);
        case {'PAR', 'par'}
            param.weights = weights_triag(0, fs/2, floor(param.M/2)+1, param.F, param.Ls, fs/2, 2);
    end
else
    param.weights = 1;
end


%% Proximal algorithm
fprintf('Starting the proximal algorithm\n')

tic;
switch param.algorithm
    case {'Douglas-Rachford', 'DR'}
        if param.reweighting
            [data_rec, dSDR_process, objective_process] = reweighted_douglas_rachford(data_clipped, param, paramsolver, data);
        else
            [data_rec, dSDR_process, objective_process] = douglas_rachford(data_clipped, param, paramsolver, data);
        end
    case {'Chambolle-Pock', 'CP'}
        if param.reweighting
            [data_rec, dSDR_process, objective_process] = reweighted_chambolle_pock(data_clipped, param, paramsolver, data);                        
        else
            [data_rec, dSDR_process, objective_process] = chambolle_pock(data_clipped, param, paramsolver, data);
        end
    otherwise
        error('Invalid algorithm is set!');
end
time = toc;

% replace reliable samples with original from clipped signal and compute error on reliable samples
if param.replaceReliable
    reliableDifference = norm(data_clipped(param.masks.Mr) - data_rec(param.masks.Mr));
    data_rec(param.masks.Mr) = data_clipped(param.masks.Mr);
    fprintf('l2-norm of difference on reliable samples after proximal algorithm is %4.3f \n', reliableDifference);
end


%% Time & SDR evaluation

% time
fprintf('Result obtained in %4.3f seconds.\n', time);

% SDR
sdr_clip = sdr(data, data_clipped);
sdr_rec = sdr(data, data_rec);
fprintf('SDR of the clipped signal is %4.3f dB.\n', sdr_clip);
fprintf('SDR of the reconstructed signal is %4.3f dB.\n', sdr_rec);
fprintf('SDR improvement is %4.3f dB.\n', sdr_rec-sdr_clip);





