function [ReconstSignal1 ReconstSignal2] = inpaintSignal_IndependentProcessingOfFrames(problemData,param)
%
%
% Usage:
%
%
% Inputs:
%          - 
%          - 
%          - 
%          - 
%          - 
%          - 
%          - 
%          - 
%
% Outputs:
%          - 
%          - 
%          - 
%          - 
%
% Note that the CVX library is needed.
%
% -------------------
%
% Audio Inpainting toolbox
% Date: June 28, 2011
% By Valentin Emiya, Amir Adler, Maria Jafari
% This code is distributed under the terms of the GNU Public License version 3 (http://www.gnu.org/licenses/gpl.txt).
% ========================================================
% Perform Audio De-clipping with overlapping blocks
% Synthesis Approach, union of overcomplete DCT dictionary
% Date: 14 Apr. 2010
% Inputs:
%          - x: Clipped signal
%          - ClipMask: Indices of clipped samples
%          - Optional parameters [and default values]:
%             - param.N: frame length [256]
%             - param.frameOverlapFactor: overlap factor between frames [2]
%             - param.wa: weighting analysis window [@sineWin]
%             - param.ws: weighting synthesis window [@sineWin]
%             - param.OMPerr: error threshold to stop OMP iterations [0.001]
%             - param.sparsityDegree: max number of non-zero components to
%             stop OMP iterations [param.N/4];
%             - other fields: see the documentation of UDCT_Dictionary
%
% Outputs:
%          ReconstSignal1 - reconstructed signal (all samples generated
%          from the synthesis model)
%          ReconstSignal2 - reconstructed signal (only clipped samples are generated
%          from the synthesis model)
%
% By Valentin Emiya - SMALL Project, 2010
% 
% ========================================================

% Check parameters
defaultParam.N = 256;
defaultParam.OLA_frameOverlapFactor = 4;
defaultParam.wa = @wSine;
defaultParam.OLA_ws = @wSine;
defaultParam.OLA_par_waitingTime_mainProcess = 0.2;
defaultParam.OLA_par_waitingTime_thread = 0.2;
defaultParam.OLA_par_frameBlockSize = 1;
defaultParam.TCPIP_port = 3000;
defaultParam.COM_DISP = false;
defaultParam.STATE_DISP = false;

if nargin<2
    param = defaultParam;
else
    names = fieldnames(defaultParam);
    for k=1:length(names)
        if ~isfield(param,names{k}) || isempty(param.(names{k}))
            param.(names{k}) = defaultParam.(names{k});
        end
    end
end

x = problemData.x;
ClipMask = find(problemData.IMiss);

% According to this flag, switch between a parallel multithread processing
% and a singlethread processing. This can fasten the computations but does
% not affect the results.
if param.MULTITHREAD_FRAME_PROCESSING
   [ReconstSignal1 ReconstSignal2] = multithreadProcessing(x,ClipMask,param);
else
   [ReconstSignal1 ReconstSignal2] = singlethreadProcessing(x,ClipMask,param);
end

return;

function [ReconstSignal1 ReconstSignal2] = singlethreadProcessing(x,ClipMask,param)
% ========================================================
% Overlap-and-add processing of a signal with missing samples.
% Decomposition into overlapping frames, processing of each
% frame independently and OLA reconstruction.
% Date: 01 Jun. 2010
% Inputs:
%          - x: Clipped signal
%          - ClipMask: Indices of clipped samples
%          - Optional parameters [and default values]:
%             - param.N: frame length [256]
%             - param.inpaintFrame: function handle for inpainting a frame
%             [@inpaintFrame_OMP]
%             - param.OLA_frameOverlapFactor: overlap factor between frames [2]
%             - param.wa: weighting analysis window [@sineWin]
%             - param.OLA_ws: weighting synthesis window [@sineWin]
%             - param.SKIP_CLEAN_FRAMES: flag to skip frames with no
%             missing samples [true]
%             - other fields: see the documentation the inpainting method
%
% Outputs:
%          ReconstSignal1 - reconstructed signal (all samples generated
%          from the synthesis model)
%          ReconstSignal2 - reconstructed signal (only clipped samples are generated
%          from the synthesis model)
%
% By Amir Adler, Maria Jafari, Valentin Emiya - SMALL Project, 2010
% 
% ========================================================

% Check parameters
defaultParam.N = 256;
defaultParam.inpaintFrame = @inpaintFrame_OMP;
defaultParam.OLA_frameOverlapFactor = 2;
defaultParam.wa = @wSine;
defaultParam.ws = @wSine;
defaultParam.SKIP_CLEAN_FRAMES = true;

if nargin<3
    param = defaultParam;
else
    names = fieldnames(defaultParam);
    for k=1:length(names)
        if ~isfield(param,names{k}) || isempty(param.(names{k}))
            param.(names{k}) = defaultParam.(names{k});
        end
    end
end
% if ~isfield(param,'D')
%    param.D = param.D_fun(param);
% end

bb=param.N; % block size

% modify signal length to accommodate integer number of blocks
L=floor(length(x)/bb)*bb;
x=x(1:L);
ClipMask(ClipMask>L) = [];

% Extracting the signal blocks with 50% overlap
Ibegin = (1:bb/param.OLA_frameOverlapFactor:length(x)-bb);
if Ibegin(end)~=L-bb+1
    Ibegin(end+1) = L-bb+1;
end
Iblk = ones(bb,1)*Ibegin+(0:bb-1).'*ones(size(Ibegin));
wa = param.wa(bb); % analysis window
xFrames=diag(wa)*x(Iblk);

% Generating the block mask
Mask=ones(size(x));
Mask(ClipMask)=0;
blkMask=Mask(Iblk);

% Declipping the Patches
[n,P]=size(xFrames);

Reconst = zeros(n,P);
for k=1:1:P,

    if param.SKIP_CLEAN_FRAMES && all(blkMask(:,k))
        Reconst(:,k) = xFrames(:,k);
        continue
    end
    frameProblemData.x = xFrames(:,k);
    frameProblemData.IMiss = ~blkMask(:,k);
    
    Reconst(:,k) = ...
        param.inpaintFrame(frameProblemData,param);
    
    if param.VERBOSE
        fprintf('  Processed signal blocks: %d/%d \n', k, P)
    end
    
end;

% Overlap and add

% The completly reconstructed signal
ReconstSignal1 = zeros(size(x));
ws = param.OLA_ws(bb); % synthesis window
wNorm = zeros(size(ReconstSignal1));
for k=1:size(Iblk,2)
    ReconstSignal1(Iblk(:,k)) = ReconstSignal1(Iblk(:,k)) + Reconst(:,k).*ws(:);
    wNorm(Iblk(:,k)) = wNorm(Iblk(:,k)) + ws(:).*wa(:);
end
ReconstSignal1 = ReconstSignal1./wNorm;

% Only replace the clipped samples with the reconstructed ones
ReconstSignal2=x;
ReconstSignal2(ClipMask)=ReconstSignal1(ClipMask);

return;

function [ReconstSignal1 ReconstSignal2] = multithreadProcessing(x,ClipMask,param)
% Send parameters to the threads
% initParamFilename = [param.OLA_par_threadDir 'par_param.mat'];
% save(initParamFilename,'param');

bb=param.N; % block size

% modify signal length to accommodate integer number of blocks
L=floor(length(x)/bb)*bb;
x=x(1:L);
ClipMask(ClipMask>L) = [];

% Extracting the signal blocks with 50% overlap
Ibegin = (1:round(bb/param.OLA_frameOverlapFactor):length(x)-bb);
if Ibegin(end)~=L-bb+1
    Ibegin(end+1) = L-bb+1;
end
Iblk = ones(bb,1)*Ibegin+(0:bb-1).'*ones(size(Ibegin));
wa = param.wa(bb); % analysis window
xFrames=diag(wa)*x(Iblk);

% Generating the block mask
Mask=ones(size(x));
Mask(ClipMask)=0;
blkMask=Mask(Iblk);

% Declipping the Patches
[n,P]=size(xFrames);
Reconst = NaN(n,P);
% initializedThreads = [];

% find the block of frames to process
k_lists = {};
kTrame = 1;
while kTrame<=P
    k_list = zeros(param.OLA_par_frameBlockSize,1);
    ind = 0;
    while ind<param.OLA_par_frameBlockSize && kTrame<=P
        if param.SKIP_CLEAN_FRAMES && all(blkMask(:,kTrame))
            kTrame=kTrame+1;
            continue
        end
        ind = ind+1;
        k_list(ind) = kTrame;
        kTrame=kTrame+1;
    end
    if ind==0
        break;
    end
    k_lists{end+1} = k_list(1:ind);
end
k_list_all = cell2mat(k_lists');

% Create a server
serverSocket = createServer(param);
% Definition of the client states
stateDef;

param.COM_DISP = false;

kBlock=1;
initializedClientIDs = [];
while any(isnan(Reconst(1,k_list_all)))
    
    % Wait for a new client
    currentClient = waitClient(serverSocket);
    
    % Receive client state
    clientIDState = readData(currentClient,param.COM_DISP);
    clientID = clientIDState(1);
    clientState = clientIDState(2);
    switch clientState
        case INIT
            if param.STATE_DISP
                fprintf('INIT %d\n',clientID);
            end
            if 0
               sendData(currentClient,initParamFilename,param.COM_DISP);
            else
               sendData(currentClient,param,param.COM_DISP);
            end
            initializedClientIDs(end+1) = clientID;
        case FREE
            if ~ismember(clientID,initializedClientIDs)
                sendData(currentClient,INIT_ORDER,param.COM_DISP); % INIT
            elseif kBlock<=length(k_lists)
                k_list = k_lists{kBlock};
                if param.STATE_DISP
                    fprintf('TO PROCESS %d:',clientID);
                    arrayfun(@(x)fprintf(' %d',x),k_list);
                    fprintf('\n');
                end
                sendData(currentClient,k_list,param.COM_DISP);
                sendData(currentClient,xFrames(:,k_list),param.COM_DISP);
                sendData(currentClient,find(blkMask(:,k_list)),param.COM_DISP);
                kBlock = kBlock+1;
            else
                if param.STATE_DISP
                    fprintf('WAIT %d\n',clientID);
                end
                sendData(currentClient,WAIT_ORDER,param.COM_DISP); % no data to process
            end
        case PROCESSED
            processed_k_list = readData(currentClient,param.COM_DISP);
            y = readData(currentClient,param.COM_DISP);
            y = reshape(y,[],length(processed_k_list));
            if param.STATE_DISP
                fprintf('PROCESSED %d:',clientID);
                arrayfun(@(x)fprintf(' %d',x),processed_k_list);
                fprintf('\n');
            end
            if ~isempty(processed_k_list)
                Reconst(:,processed_k_list) = y;
            end
        otherwise
            error('switch:UndefinedClientState','Undefined client state');
    end
    
    closeClientConnection(currentClient);
end;

% Close the server
closeServer(serverSocket);


% Overlap and add

% The completly reconstructed signal
ReconstSignal1 = zeros(size(x));
ws = param.OLA_ws(bb); % synthesis window
wNorm = zeros(size(ReconstSignal1));
for k=1:size(Iblk,2)
    ReconstSignal1(Iblk(:,k)) = ReconstSignal1(Iblk(:,k)) + Reconst(:,k).*ws(:);
    wNorm(Iblk(:,k)) = wNorm(Iblk(:,k)) + ws(:).*wa(:);
end
ReconstSignal1 = ReconstSignal1./wNorm;

% Only replace the clipped samples with the reconstructed ones
ReconstSignal2=x;
ReconstSignal2(ClipMask)=ReconstSignal1(ClipMask);

killClients(param);

return;

% ========================================================
% ========================================================



% ========================================================
% ========================================================
