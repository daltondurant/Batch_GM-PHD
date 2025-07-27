%% batchls
%**************************************************************************
% author: dalton durant
% email: ddurant@utexas.edu, thedaltondurant@gmail.com
%
% USE AT YOU OWN DISCRETION
%
% This code is free to use, all I ask is you cite me when reasonable.
%**************************************************************************
% Batch Least Squares Information Filter
%
% SOURCE: 
% [1] Bob Schutz, Byron Tapley, and George H. Born. Statistical Orbit 
% Determination. Elsevier. 2004.
function [Zk, Rk] = batchls(cfig,model,Y,xbar0,varargin)
%
    %--- parameters
    epsilon = 1e-6; % convergence threshold
    maxiter = 1e3;  % maximum number of Batch LS while loop iterations 

    %--- initialize
    T0 = xbar0(1);
    xbar0 = xbar0(2:end);
    TY = Y(1,:);
    Y  = Y(2:end,:);
    nx = size(xbar0,1);
    ny = size(Y,1);
    nmeas = size(Y,2);
    iter = 0;
    xhat = xbar0;   
    deltaxbar = zeros(nx,1);
    RMS_normalized = 1000; 
    RMS_normalized_PREV = 0;
    dRMS = abs(RMS_normalized - RMS_normalized_PREV);
    Vr = chol(model.R); inv_sqrt_R= inv(Vr); iR= inv_sqrt_R*inv_sqrt_R'; % stable inversion
    
    %--- default is diffuse 
    lambda0 = zeros(nx,nx); % default, purely diffuse prior
    %lambda0 = 1e-3*eye(nx); % diffuse prior (relative to the problem, but keeps sim from breaking)
    ilambda0 = 1e6 * eye(nx);

    %--- loading optional covariance argument
    while ~isempty(varargin)
        if isnumeric(varargin{1})
            Pbar0 = varargin{1};
            if sum(Pbar0, 'all') ~= 0 && size(Pbar0,1) == nx && size(Pbar0,2) == nx
                Vp = chol(Pbar0); inv_sqrt_P= inv(Vp); iP= inv_sqrt_P*inv_sqrt_P'; % stable inversion
                lambda0 = iP; % regularized maximum likelihood 
                ilambda0 = Pbar0;
            end
        end
        varargin(1) = [];
    end
    
    while dRMS > epsilon % LS while loop
        % initialize 
        xprop = xhat;
        phi_prop = eye(nx);
        lambda = lambda0; 
        N = lambda*deltaxbar;
        dY_array = zeros(ny,nmeas);
        Tjj = T0;
        for jj = 1:nmeas
            % propagate
            model.T = TY(jj)-Tjj;
            phi_prop = expm(cfig.F(xprop,model)*model.T) * phi_prop;
            dT = 1; h  = model.T / dT; 
            for hh = 1:dT % rk4
                k1 = cfig.f(0, xprop, model);
                k2 = cfig.f(0, xprop+1/2.*h.*k1, model);
                k3 = cfig.f(0, xprop+1/2.*h.*k2, model);
                k4 = cfig.f(0, xprop+h.*k3, model);
                deltax = 1/6*h*(k1+2*k2+2*k3+k4);
                xprop = xprop + deltax;
            end
            %xprop = cfig.fprop(model, xprop, 'noiseless');

            % update
            y    = Y(:,jj);
            yhat = cfig.h(model, real(xprop), 'noiseless');
            deltay = y - yhat;  
            dY_array(:,jj)  = deltay; 
            H    = cfig.H(model,xprop);
            HP = H * phi_prop; % H(start) = H(k)*Phi(k,start) 
            lambda = lambda + HP' * (iR * HP); % Lambda(start) 
            N = N + HP' * (iR * deltay); % N(start)

            Tjj = TY(jj);
        end

        % solve normal equations
        ilambda = pinv(lambda);
        deltaxhat = ilambda * N; 
        xhat = xhat + deltaxhat;      
        deltaxbar = deltaxbar - deltaxhat;   
        iter = iter + 1; % update iteration value
        
        % test convergence
        rms_temp = zeros(ny, nmeas);
        for jj = 1:nmeas
            rms_temp(:,jj) = iR * dY_array(:,jj);
        end
        RMS_normalized = sqrt(dY_array(:)' * rms_temp(:) / (ny*nmeas));
        dRMS = abs(RMS_normalized - RMS_normalized_PREV); % change in RMS
        RMS_normalized_PREV = RMS_normalized; % set the previous RMS as the current RMS

        if iter == maxiter
            fprintf(['\nMax iterations reached. Batch LS failed to converge. \n'])
            break
        end
        if any(isnan(xhat)) || ~all(isreal(xhat))
            ilambda = ilambda0;
            xhat = xbar0 + sqrtm(ilambda)*randn(model.stream, [model.x_dim,1]);
            fprintf(['\nNaNs. Batch LS failed to converge. \n'])
            break
        end

    end % end LS while loop

    Zk   = xhat;
    Rk   = (ilambda+ilambda')./2;
%
end 

