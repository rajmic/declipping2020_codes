function yFrames = signal2frames(y, param)

% Divide signal into overlapping time-frames
% 
% Input: -y: signal
%        -param.N: frame length
%         param.hop: hop size
%         param.wa: Analysis window
% 
% Output: -yFrames: matrix containing each frame
% 

N = param.N;
L = length(y);
hop = param.hop;
wa = param.wa(N);

istart = 1:hop:(L-N+1); % start index of each frame
iframes = bsxfun(@plus,repmat(istart,N,1),(0:N-1)'); % index map of each frame

yFrames = diag(wa) * y(iframes); % Overlapping time frames

end
