function [D,A,alg_log] = consDL(Y,paramDL,paramSC,paramDictUpdate, distortion_param)

    %% Global parameters dictionary learning

    A = paramDL.A_init;
    D = paramDL.D_init;

    %% DL parameters

    Nit = paramDL.Nit;

    alg_log.cost = NaN * zeros(1,2*Nit);

    %% DL iterations

    it = 0;

    while it < paramDL.Nit

        it = it+1;

        %% Sparse coding

        if paramDL.warm_start
            paramSC.A_init = A; % warm_start
        else
            paramSC.A_init = zeros(size(A));
        end

        [A, SC_log] = paramSC.alg(Y,D,paramSC,distortion_param);

        cost = SC_log.final_cost;

        alg_log.cost(2*it-1) = cost;

        if paramDL.loud
            fprintf('it %d, SC: cost: %.3f\n',it, cost)
        end

        %% Prune unused atoms:

        unused = (sum(A.^2,2) == 0);

        if sum(unused)>0 && paramDL.loud
            fprintf('  %d atoms pruned!\n', sum(unused))
        end

        D(:,unused) = [];
        A(unused,:) = [];

        %% Dictionary Update

        paramDictUpdate.D_init = D;

        [D, DictUpdate_log] = consDictUpdate(Y,A,paramDictUpdate,distortion_param);

        cost = DictUpdate_log.final_cost;
        alg_log.cost(2*it) = cost;

        if paramDL.loud
            fprintf('it %d, DL: cost: %.3f\n',it, DictUpdate_log.final_cost)
        end


    end % end DL

end

function [D,alg_log] = consDictUpdate(Y,A,paramDictUpdate,distortion_param)


    %% initialization:
    D = paramDictUpdate.D_init;

    % residual:
    X_tmp = D*A;
    ResidualMat = projection(X_tmp, Y, distortion_param)-X_tmp;

    alg_log = [];
    alg_log.cost(1) = 1/2 * sum(sum(ResidualMat.^2)); 

    mu_opt = 0.99*1/norm(A)^2;

    %% Gradient descent:
    
    if paramDictUpdate.accelerate == 0

        for it = 1:paramDictUpdate.Nit

            % gradient descent:
            D = D + mu_opt * ResidualMat*A';

            % normalize atoms:
            D = bsxfun(@times,D,1./max(sqrt(sum(D.^2,1)),1));

            % update residual: 
            X_tmp = D*A;
            ResidualMat = projection(X_tmp, Y, distortion_param)-X_tmp;

            alg_log.cost(it+1) = 1/2 * sum(sum(ResidualMat.^2)); 

        end
        
    else
        
        % Fista parameters:
        U_old = D;
        t_old = 1;
        
        for it = 1:paramDictUpdate.Nit
                  
            U = D + mu_opt * ResidualMat*A';
            U = bsxfun(@times,U,1./max(sqrt(sum(U.^2,1)),1));

            t = (1+sqrt(1+4*t_old^2))/2;

            D = U + (t_old-1)/t*(U-U_old);

            U_old = U;
            t_old = t; 

            % update residual:
            X_tmp = D*A;
            ResidualMat = projection(X_tmp, Y, distortion_param)-X_tmp;
        
            alg_log.cost(it+1) = 1/2 * sum(sum(ResidualMat.^2));
    
        end
        
    end
            
    alg_log.final_cost = 1/2 * sum(sum(ResidualMat.^2));

end




