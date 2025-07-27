%% hmc
%**************************************************************************
% author: dalton durant
% email: ddurant@utexas.edu, thedaltondurant@gmail.com
%
% USE AT YOU OWN DISCRETION
%
% This code is free to use, all I ask is you cite me when reasonable.
%**************************************************************************
% Markov Chain Monte Carlo 
% Hamiltonian Monte Carlo
% Can handle time history of measurements such that y is a 2D matrix from 
% t0 --> tf. Will return samples at t0 like a batch filter.
%
% SOURCE: 
% 
function [varargout] = hmc(cfig,model,y,m,P)

    %--- parameters
    L       = 30;    % lags/steps in Markov Chain 
    Mburn   = 10;    % number of burned MCMC samples
    maxiter = 1e3;   % maximum number of MCMC iterations 
    LL      = 10;    % pseudotime leapfrog steps
    epsilon = 1e-1;  % HMC step size

    %--- initialize the Markov Chain
    acceptances = 0;
    iter        = 0;
    xbar0 = m;
    Pbar0 = P;
    nx      = model.x_dim;
    X_temp      = repmat(xbar0(2:end), 1, L+1);
    q0 = xbar0(2:end); 

    %--- generate the samples
    Mass = pinv(Pbar0); % constant Mass matrix   
    while acceptances < L && iter < maxiter
        %--- propose a new sample using leapfrog
        q = q0;
        p = sqrtm(Mass)*randn(model.stream, [nx,1]); % independent standard normal variates
        p0 = p;
        % evaluate potential energy at start of trajectory
        [current_U,gradU] = U_func(cfig, model, q0, y, xbar0);
        % evaluate kinetic energy at start of trajectory
        current_K  = K_func(p0, Mass); 
        % make a half step for momentum at the beginning
        p = p - epsilon * gradU / 2;
        % alternate full steps for position and momentum
        for ell = 1:LL
            % make a full step for the position
            [~,gradK] = K_func(p, Mass);
            q = q + epsilon * gradK;
            % make a full step for the momentum, except at end of trajectory
            if (ell < LL)
                [~,gradU] = U_func(cfig, model, q, y, xbar0);
                p = p - epsilon * gradU;
            end
        end
        % evaluate potential energy at end of trajectory
        [proposed_U,gradU] = U_func(cfig, model, q, y, xbar0);
        % make a half step for momentum at the end.
        p = p - epsilon * gradU / 2;  
        % negate momentum at end of trajectory to make the proposal symmetric
        p = -p;
        % evaluate kinetic energy at end of trajectory
        proposed_K = K_func(p, Mass);

        %--- compute the acceptance ratio
        accept_prob = min(1,exp(current_U-proposed_U+current_K-proposed_K));

        %--- accept or reject the proposed sample at end of trajectory
        if rand(model.stream, 1) <= accept_prob
            % accepted 
            q0 = q;  
            acceptances = acceptances + 1;
            epsilon = epsilon * 1.2; % inc by 20%
            X_temp(:,acceptances+1) = q0;
        else
            % rejected
            epsilon = epsilon / 1.2; % dec by 20%
        end
        iter = iter + 1;

        if iter == maxiter
            fprintf(['\nMax iterations reached. Batch MCMC failed. \n'])
            break
        end

    end % end while loop

    %--- outputs
    Msamp   = L - Mburn; % number of requested MCMC samples
    betaS = (4/(Msamp*(nx+2))) ^ (2/(nx+4)); % Silverman's rule of thumb
    Xmcmc   = X_temp(:,Mburn+2:end);
    X_bar = mean(Xmcmc,2); eX = squeeze(Xmcmc - X_bar); P_temp  = (eX * eX') / (Msamp-1);  % sample covariance
    Pmcmc   = repmat(betaS * P_temp, 1, 1, Msamp);
    accept_rate        = acceptances*100/iter; 
    
    varargout{1} = Xmcmc;
    varargout{2} = Pmcmc;
    varargout{3} = accept_rate;
%
end % end mcmc_func


%% Helper functions
% target distribution
function varargout = U_func(cfig, model, q, y, xbar)

    T0 = xbar(1);
    Ty = y(1,:);
    y  = y(2:end,:);
    nx    = size(q,1);
    ny    = size(y,1);
    nmeas = size(y,2);

    % --- precompute inversion
    det_R = prod(eig(model.R)); iR = pinv(model.R);

    qprop = q;
    phi_prop = eye(nx);
    lambda = zeros(nx,nx); % purely diffuse prior
    %lambda = iP;
    grad_pygx = zeros(nx,1);  
    dy_array = zeros(ny,nmeas);
    Tkk = T0;
    for kk = 1:nmeas
        % propagate
        model.T = Ty(kk)-Tkk;
        phi_prop = expm(cfig.F(qprop,model)*model.T) * phi_prop;
        dT = 1; h  = model.T / dT;
        for hh = 1:dT
            k1 = cfig.f(0, qprop, model);
            k2 = cfig.f(0, qprop+1/2.*h.*k1, model);
            k3 = cfig.f(0, qprop+1/2.*h.*k2, model);
            k4 = cfig.f(0, qprop+h.*k3, model);
            deltaq = 1/6*h*(k1+2*k2+2*k3+k4);
            qprop = qprop + deltaq;
        end
        %qprop = cfig.fprop(model, qprop, 'noiseless');

        % update
        yk    = y(:,kk);
        yhat = cfig.h(model, real(qprop), 'noiseless');
        deltay = yk - yhat;  
        dy_array(:,kk)  = deltay; 
        H    = cfig.H(model,qprop);
        HP = H * phi_prop; % H(start) = H(k)*Phi(k,start) 
        lambda = lambda + HP' * (iR * HP); % Lambda(start) 
        grad_pygx = grad_pygx + HP' * (iR * deltay); % N(start)

        Tkk = Ty(kk);
    end

    % --- computing -log of likelihood
    utemp = zeros(ny, nmeas);
    for kk = 1:nmeas
        utemp(:,kk) = iR * dy_array(:,kk);
    end

    nlog_pygx =   dy_array(:)' * utemp(:) / 2 ...
                + nmeas*log(det_R) / 2 ...
                + nmeas*log(2*pi) * ny / 2; 

    % --- output
    varargout{1}     = nlog_pygx;  % -log of likelihood
    varargout{2}     = -grad_pygx; % gradiant   
    varargout{3}     = lambda;     % Fisher Information Matrix (FIM)

end % end function

% proposal distribution
function varargout = K_func(p, Mass)
    nx = size(p,1);
    det_M = prod(eig(Mass)); iM = pinv(Mass);

    % --- computing -log of proposal
    nlog_proposal =   dot(p, iM * p)' / 2 ...
                    + log(det_M) / 2 ...
                    + log(2*pi) * nx / 2;

    % --- output
    varargout{1} = nlog_proposal; % -log of proposal 
    varargout{2} = iM * p;

end % end function