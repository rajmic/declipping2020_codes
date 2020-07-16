%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%                                %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%      DECLIPPING MAIN FILE      %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%         Social Sparsity        %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%                                %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% The implementation of the Social Sparsity Declipper was kindly provided
% by Matthieu Kowalski.
% The codes were then slightly edited to fit the common interface of the
% declipping toolbox by
% Pavel Záviška, Brno University of Technology, 2020
%
% using toolbox LTFAT
ltfatstart;

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
Ls = length(data);


%% settings
fprintf('Setting up the frame parameters\n')

% input SDR of the clipped signal
inputSDR = 7;    % set the input SDR value

% DGT parameters
wtype = 'hann';  % window type
w = 8192;        % window length
a = w / 4;       % window shift
M = 2*8192;      % number of frequency channels

% set shrinkage operator
shrinkage = 'EW'; % 'L', 'WGL', 'EW', 'PEW';
fprintf('Setting up shrinkage operator: %s\n', shrinkage)

verbose = true;

      
%% clipping
fprintf('Generating clipped signal\n')
[data_clipped, masks, theta, trueSDR, percentage] = clip_sdr(data, inputSDR); % clipping of original signal
fprintf('Clipping threshold %4.3f, true SDR value is %4.2f dB and %4.2f%% samples has been clipped \n', theta, trueSDR, percentage)


data_r = zeros(Ls,1);
data_r(masks.Mr) = data_clipped(masks.Mr);

data_c = zeros(Ls,1);
data_c(masks.Mh) = theta;
data_c(masks.Ml) = -theta;
        

%% construction of the frame
g = gabwin({'tight', {'hann', w}}, a, M);

G_clip = dgtreal(data_clipped, g, a, M);
[nf,nt] = size(G_clip);

gab.analysis = @(x) dgtreal(x,g,a,M);
gab.synthesis = @(x) idgtreal(x,g,a,M,Ls);
        
%% Setting of the Shrinkage operator

switch shrinkage
    case 'L'
        shrink = struct('type', 'l', 'lambda', 1e-4, 'mu', 1, 'glabel', 'time','neigh', 1, 'center', [1 1], 'expo', 1,'orth',0);
    case 'WGL'
        shrink = struct('type', 'l', 'lambda', 1e-4, 'mu', 1, 'glabel', 'time','neigh', ones(3,7), 'center', [2,4], 'expo', 1,'orth',0);
    case 'EW'
        shrink = struct('type', 'l', 'lambda', 1e-4, 'mu', 1, 'glabel', 'time','neigh', 1, 'center', [1 1], 'expo', 2,'orth',0);
    case 'PEW'
        shrink = struct('type', 'l', 'lambda', 1e-4, 'mu', 1, 'glabel', 'time','neigh', ones(3,7), 'center', [2,4], 'expo', 2,'orth',0);
    otherwise
        error('Invalid shrinkage operator!')
end

Gdeclip = zeros(nf,nt);
GdeclipZ = Gdeclip;
        
%% ISTA
fprintf('Starting the ISTA algoritm\n')

tic

number_lambdas = 20;
inner_iterations = 500;
delta = 1e-4;
iter_cnt = 0;

dSDR_process = NaN(number_lambdas*inner_iterations, 1);
data_rec = data_clipped;

for lambda=logspace(-1,-4,number_lambdas)
    shrink.lambda = lambda;
    for k=1:inner_iterations

        data_rec_old = data_rec;
        iter_cnt = iter_cnt + 1;

        % forward step
        GdeclipOLD = Gdeclip;

        xdeclipZ = gab.synthesis(GdeclipZ);

        r1 = zeros(Ls,1);
        r1(masks.Mr) = data_r(masks.Mr) - xdeclipZ(masks.Mr);

        grad1 = -gab.analysis(r1);

        r2 = zeros(Ls,1);
        r2(masks.Mh) = data_c(masks.Mh) - xdeclipZ(masks.Mh);
        r2(masks.Ml) = data_c(masks.Ml) - xdeclipZ(masks.Ml);

        r2(abs(xdeclipZ)>theta) = 0;
        grad2 = -gab.analysis(r2);


        GdeclipZ = GdeclipZ - grad1 - grad2;


        % thresholding step
        Gdeclip = gen_thresh(GdeclipZ, shrink);


        % relaxation step
        GdeclipZ = Gdeclip + (k-1)/(k+5) * (Gdeclip - GdeclipOLD);



        data_rec = gab.synthesis(Gdeclip);
        if norm(data_rec_old - data_rec) < delta
            break
        end


        if verbose
            sdr_rec = sdr(data, data_rec);
            fprintf('  lambda = %f -- it = %d -- SDR = %f -- dSDR = %f \n',lambda,k,sdr_rec,sdr_rec-inputSDR);
        end

    end
end

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
