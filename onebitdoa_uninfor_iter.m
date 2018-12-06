function [ thetadegree_est, NMSE_SBL,X_debiased ] = onebitdoa_uninfor_iter( N, L, M, K, X, Y, wvar, maxit_outer )
% One bit DOA estimator utilizing SBL
% Written by Jiang Zhu and Xiangming Meng
% 2017, Nov. 26

% Input: A: Steering matrix
% X: True Complex coefficients
% N: Grid size
% L: Number of snapshots
% M: Number of attennas
% K: Number of sources. For DOA estimation, it is not needed.    
    ct = 0;
    % ULA-horizontal
    d = 1/2;                % intersensor spacing
    q = 0:1:(M-1);        % sensor numbering
    xq = (q-(M-1)/2)*d;     % sensor locations

    % Bearing grid
    theta = (-90:180/(N-1):90);
    thetadeg = theta;
    theta_r = theta*pi/180;
    u = sin(theta_r);
    % Represenation matrix (steering matrix)
    A = exp(-1i*2*pi*xq'*u)/sqrt(N); % M*N

    NMSE_SBL = zeros(maxit_outer,1);
    % input initialization 
    z_A_ext = zeros(M,L);
    Lar_num = 1e20;
    v_A_ext = Lar_num*ones(M,L);
    X_est = zeros(size(X));
    X_est_old = zeros(size(X));
    var_min = 1e-20;
    var_max = 1e20;
    
    % parameter initialization
    alpha = 1e-1*ones(N,1);
    alpha_all = zeros(N,L);
    a = 1+eps;
    b = 0+eps;

    for i = 1:maxit_outer
%         waitbar(i/maxit_outer)
        v_A_ext = var_max*(v_A_ext<=0)+v_A_ext.*(v_A_ext>0);
        v_A_ext = min(v_A_ext,var_max);
        v_A_ext = max(v_A_ext,var_min);
        % set up the parallel computing toolbox
        for snap = 1:L
            v_A_ext_snap = v_A_ext(:,snap);
            z_A_ext_snap = z_A_ext(:,snap);
            y = Y(:,snap);
            % transforming to real observations
            z_A_ext_real = [real(z_A_ext_snap);imag(z_A_ext_snap)];
            y_real = [real(y);imag(y)];
            v_A_ext_snap_real = [v_A_ext_snap;v_A_ext_snap]/2;
            [z_B_post_real, v_B_post_real] = GaussianMomentsComputation(y_real, 0, z_A_ext_real, v_A_ext_snap_real, wvar/2);
            v_B_post_snap = v_B_post_real(1:M)+v_B_post_real(M+1:end);
            z_B_post_snap = z_B_post_real(1:M)+1j*z_B_post_real(M+1:end);
            v_B_ext_snap = v_B_post_snap.*v_A_ext_snap./(v_A_ext_snap-v_B_post_snap+eps);
            z_B_ext_snap = v_B_ext_snap.*(z_B_post_snap./v_B_post_snap-z_A_ext_snap./v_A_ext_snap);
    %         sum(isnan(z_B_ext_snap))
            v_B_ext_snap = var_max*(v_B_ext_snap<=0)+v_B_ext_snap.*(v_B_ext_snap>0);
            v_B_ext_snap = min(v_B_ext_snap,var_max);
            v_B_ext_snap = max(v_B_ext_snap,var_min);

            beta = 1./v_B_ext_snap;
            y_tilde = z_B_ext_snap;
            Sigma = inv((A'*diag(beta)*A)+diag(alpha));   %  Sigma: posterior variance
            Sigma = (Sigma+Sigma')/2;
            mu = Sigma*A'*diag(beta)*y_tilde;
            alpha_update = a./(mu.*conj(mu)+diag(Sigma)+b);
            X_est(:,snap) = mu;
            alpha_all(:,snap) = alpha_update;
            z_A_post_snap = A*mu;
            v_A_post_snap = diag(A*Sigma*A');
            v_A_post_snap = (v_A_post_snap+conj(v_A_post_snap))/2;

            v_A_ext_snap = v_A_post_snap.*v_B_ext_snap./(v_B_ext_snap-v_A_post_snap+eps);
            z_A_ext_snap = v_A_ext_snap.*(z_A_post_snap./v_A_post_snap-y_tilde./v_B_ext_snap); 
            v_A_ext_snap = var_max*(v_A_ext_snap<=0)+v_A_ext_snap.*(v_A_ext_snap>0);

            v_A_ext_snap = min(v_A_ext_snap,var_max);
            v_A_ext_snap = max(v_A_ext_snap,var_min);

            v_A_ext(:,snap) = v_A_ext_snap;
            z_A_ext(:,snap) = z_A_ext_snap;
        end

        alpha = mean(alpha_all,2);
        c = diag(X_est'*X)./diag(X_est'*X_est);
        NMSE_SBL(i) = 20*log10(norm(X_est*diag(c)-X,'fro')/norm(X,'fro'));
        if(norm(X_est_old-X_est,'fro')/norm(X_est)<1e-3)
            ct = ct+1;
            if(ct==3)
                break;
            end
        end
        X_est_old = X_est;
    end
    X_debiased = X_est*diag(c);
    [~,index_X] = sort(abs(X_est*X_est'),'descend');
    thetadegree_est = thetadeg(index_X(1:K));
end

