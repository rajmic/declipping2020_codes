clear variables;close all;clc;

ltfatstart;

[x,fs] = audioread('starwars_demo_o.wav');

x = resample(x,16000,fs);
T = length(x);

clip_level = 0.2;
x_clip = x;

Clip_plus = x >= clip_level;
Clip_minus = x <= -clip_level;
x_clip(Clip_plus) = clip_level;
x_clip(Clip_minus) = -clip_level;
reliable = abs(x) < clip_level;
x_r = zeros(T,1);
x_r(reliable) = x_clip(reliable);

x_c = zeros(T,1);
x_c(Clip_plus) = clip_level;
x_c(Clip_minus) = -clip_level;

snr_c = SNR(x(:),x_c(:));



M = 1024;
a = 256;
%g = gabwin({'tight', 'hann'}, a, M);
g = gabwin({'tight', {'hann', 1024}}, a, M);
G_clip = dgtreal(x_clip, g, a, M);
[nf,nt] = size(G_clip);

gab.analysis = @(x) dgtreal(x,g,a,M);
gab.synthesis = @(x) idgtreal(x,g,a,M,T);

shrinkL = struct('type', 'l', 'lambda', 1e-4, 'mu', 1, 'glabel', 'time','neigh', 1, 'center', [1 1], 'expo', 2,'orth',0);
shrinkWGL = struct('type', 'l', 'lambda', 1e-4, 'mu', 1, 'glabel', 'time','neigh', ones(1,7), 'center', [1,4], 'expo', 2,'orth',0);


Gdeclip = zeros(nf,nt);
GdeclipZ = Gdeclip;

shrink = shrinkWGL;

for lambda=logspace(0,-4,10)
    shrink.lambda = lambda;
    for k=1:100
        
        % forward step
        GdeclipOLD = Gdeclip;
        
        xdeclipZ = gab.synthesis(GdeclipZ);
        
        r1 = zeros(T,1);
        r1(reliable) = x_r(reliable) - xdeclipZ(reliable);
        
        grad1 = -gab.analysis(r1);
                
        r2 = zeros(T,1);
        r2(Clip_plus) = x_c(Clip_plus) - xdeclipZ(Clip_plus);
        r2(Clip_minus) = x_c(Clip_minus) - xdeclipZ(Clip_minus);
        
        r2(abs(xdeclipZ)>clip_level) = 0;
        grad2 = -gab.analysis(r2);
        
          
        GdeclipZ = GdeclipZ - grad1 - grad2;
        
        
        % thresholding step
        Gdeclip = gen_thresh(GdeclipZ, shrink);
        
        
        % relaxation step
        GdeclipZ = Gdeclip + (k-1)/(k+5) * (Gdeclip - GdeclipOLD);
        
        
        
        xdeclip = gab.synthesis(Gdeclip);
        
        snr_d = SNR(x,xdeclip);
        
        fprintf('lambda = %f -- it = %d -- SNR = %f -- gSNR = %f \n',lambda,k,snr_d,snr_d-snr_c);
         
    end
end


