function s = hard_thresholding(a, k)
% HARD_THRESHOLDING performs hard-thresholding of the DFT coefficients
% taking into account the comples conjugate coefficients.
%
% Input parameters
%       a      vector of the DFT coefficients
%       k      number of coefficients to be selected 
%
% Pavel Záviška, Brno University of Technology, 2020



odd = mod(length(a),2);

% taking only half of the spectrum + dc coefficient
a = a(1:floor(length(a)/2)+1);
a(1) = a(1)/2;

% sorting the coefficients
[~, ind] = sort(abs(a), 'descend');
s = zeros(length(a),1);
if k < length(a)
    s(ind(1:k)) = a(ind(1:k));
else
    s = a;
    warning('Variable k is greater than the length of the DFT coefficients. The coefficients will not be thresholded.')
end

% compute conjugates of selected coefficients
s(1) = s(1)*2;

if odd
    s_conj = conj(flip(s(2:end)));
else
    s_conj = conj(flip(s(2:end-1)));
end

s = [s; s_conj];

end 