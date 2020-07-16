function [A,alg_log] = consIHT(Y,D,alg_param,distortion_param)

%% Initialize parameters

if ~isfield(alg_param, 'target_cost')
    alg_param.target_cost = 0;
end

Nit = alg_param.Nit;

A_init = alg_param.A_init;

mu = 1/norm(D)^2; % gradient descent step

alg_log = [];
if alg_param.save_log
    alg_log.cost = NaN*zeros(Nit,1);
end

%% Declip

% initialize sparse coefficient matrix:
A = A_init;

X_tmp = D*A;
ResidualMat = projection(X_tmp, Y, distortion_param)-X_tmp;

cost = 1/2*sum(sum(ResidualMat.^2));
alg_log.cost(1) = cost;

it = 0;

if alg_param.rate>0
    K = 1;
else
    K = alg_param.K;
end


while it < alg_param.Nit && cost > alg_param.target_cost
    
    it = it+1;
    
    % gradient descent:
    A = A + mu * D'*ResidualMat;

    if alg_param.rate>0 && mod(it, alg_param.rate) == 0
        K = K + 1;
    end
    
    
    % thresholding:
    A = block_hard_thresholding(A, K, 'column_constraint');

    % update residual:
    X_tmp = D*A;
    ResidualMat = projection(X_tmp, Y, distortion_param)-X_tmp;
    % compute cost:
    cost = 1/2*sum(sum(ResidualMat.^2));
    alg_log.cost(it+1) = cost;
   
end

alg_log.final_cost = 1/2*sum(sum(ResidualMat.^2));
 
end
