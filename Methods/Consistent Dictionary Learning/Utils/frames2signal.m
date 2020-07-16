function y = frames2signal(yFrames, param)

% Reconstruct signal from frames using overlap and add
% 
% Input: -yFrames: overlapping time frames
%        -param.N: frame length
%         param.hop: hop size
%         param.wa: analysis window
%         param.ws: synthesis window
% 
% Output: y: Reconstructed signal
% 
% 

N = param.N;
hop = param.hop;
Nframes = size(yFrames,2);
ws = param.wa(N);
wa = param.wa(N);

L = hop*(Nframes-1)+N; % size of signal to reconstruct

istart = 1:hop:(L-N+1); % start index of each frame
iframes = bsxfun(@plus,repmat(istart,N,1),(0:N-1)'); % index map of each frame

% Initialize signal and reconstructed norm:
y = zeros(L,1);
wNorm = zeros(L,1);

for frame = 1:Nframes
    y(iframes(:,frame)) = y(iframes(:,frame)) + yFrames(:,frame).*ws(:);
    wNorm(iframes(:,frame)) = wNorm(iframes(:,frame)) + ws(:).*wa(:);
end

wNorm(wNorm==0) = eps; %avoid dividing by zero
y = y./wNorm;


end