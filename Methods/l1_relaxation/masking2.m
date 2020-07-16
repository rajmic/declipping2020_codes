function [SF, L] = masking2(x, fs, model)

% INPUTS:
% x        input signal (row signal expected!)
% fs       sampling frequency
% model    spreading function: '2slope'    = 2-slope triangle function
%                              'mpeg1'     = MPEG Psychoacoustic model 1
%                              'mpeg2'     = MPEG Psychoacoustic model 2
%                              'schroeder' = modified Schroeder function
%
% OUTPUT:
% SF       masked threshold
%
% Jiøí Schimmel and Pavel Záviška, Brno University of Technology, 2020


if (nargin<3), model='mpeg1'; end

if size(x,1) > size(x,2), x=x.'; end % ensuring row input vector x

%% FFT
N = length(x);
fa = linspace(0, fs/2, floor(N/2)+1);
X = fft((hann(N).').*x);

%% z mapping
z = 13*atan(0.76.*fa/1000)+3.5*atan((fa/7500).^2);

%% SPL computation
L = 96 + 10.*log10( 4./(N.^2) * 0.25*(abs(X(1:length(fa)))).^2 );
% Threshold in quiet
TiQ = 3.64*(fa/1000).^(-0.8) - 6.5*exp(-0.6*(fa/1000-3.3).^2) + 10.^(-3)*(fa/1000).^4;
Lm = L;
SF = 10.^(TiQ/10);

Lm(Lm<=TiQ) = -10; 

%% Masker decimation
chk = 0;
nbr = 0;

while (isempty(chk)~=1)

    % find maximum
    Lmax = max(Lm);
    Lmax_index = find(L==Lmax);

    if isempty(Lmax_index)
        break;
    end

    if(length(Lmax_index)>1), Lmax_index = Lmax_index(1,1); end
    f_Lmax_index = (Lmax_index-1) * fs / N;

    zc = z(Lmax_index);
    dz = (z-zc);

    if strcmp(model, 'mpeg1')
        % MPEG model 1
        Ft = -17.*dz + 0.15.*Lmax.*(dz-1).*((dz-1)>=0);
        Ft(dz<0) = -(6+0.4.*Lmax).*abs(dz(dz<0)) - (11+0.4.*Lmax).*(abs(dz(dz<0))-1).*((abs(dz(dz<0))-1)>=0);
        Ft = Lmax + Ft - 6.025 - 0.275*zc;
        SF = SF + 10.^(Ft/10);
    elseif strcmp(model, 'mpeg2')
        % MPEG model 2
        Ft = 15.8111389 + 7.5*(1.05*dz+0.474) - 17.5*(1+(1.05*dz+0.474).^2).^0.5 + 8*min([0 (1.05*dz-0.5).^2-2*(1.05*dz-0.5)]);
        Ft = Lmax + Ft - 6.025 - 0.275*zc;
        SF = SF + 10.^(Ft/10);
    elseif strcmp(model, 'schro')
        % modified Schroeder
        I = min([5*10.^((Lmax-96)./10) 2]);
        Ft = (15.81-I) + 7.5 * (dz+0.474) - (((1+(dz+0.474).^2).^0.5)) .* (17.5-I');
        Ft = Lmax + Ft - 14.5 - zc;
        SF = SF + 10.^(Ft/10);
    else
        %2 slope triangle
        Ft = (-27 .* abs(dz));
        Ft(dz>0) = (-27 + 0.37 .* max([Lmax-40 0])) * abs(dz(dz>0));
        Ft = Lmax + Ft - 16;
        SF = max(SF, 10.^(Ft/10));
    end

    Lm((Lm <= 10*log10(SF))) = -10; 

    z_Lmax = 13*atan(0.76.*f_Lmax_index/1000)+3.5*atan((f_Lmax_index/7500).^2);
    Lmax_index_min = find(z > (z_Lmax - 0.5),1,'first');
    Lmax_index_max = find(z < (z_Lmax + 0.5),1,'last');
    
    Lm(max(Lmax_index-2,1):min(Lmax_index+2,length(Lm))) = -10; 
            
    nbr = nbr+1;
    chk = Lm>10*log10(SF);

end

SF = 10*log10(SF);
