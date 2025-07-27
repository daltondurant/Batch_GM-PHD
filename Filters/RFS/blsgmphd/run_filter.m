%% blsgmphd
%**************************************************************************
% author: dalton durant
% email: ddurant@utexas.edu, thedaltondurant@gmail.com
%
% USE AT YOU OWN DISCRETION
%
% This code is free to use, all I ask is you cite me when reasonable.
%**************************************************************************
% A Batch Least Squares Gaussian Mixture Probability Hypothesis Density 
% filter. It uses Least Squares to batch process measurements into 
% Gaussians in the state space.
%
% SOURCE: 
% 
function [model,meas,est] = run_filter(stream,cfig,model,meas,est)
%
    % --- filter parameters
    run_flag= 'disp';             % 'disp' or 'silence' for on the fly output
    J_birth = 10;                 % number of births per time step 
    elim_threshold  = 1e-5;       % pruning threshold
    merge_threshold = 4;          % merging threshold
    cap_limit = 250;              % capping limit
    beta_c   = 1e8; % clutter intensity penalty parameter (penalizes small tracklets)
    vareps_c = 10;  % clutter decay parameter, =1 standard scaling, <1 slower decay, >1 faster decay
        
    % --- output variables
    est.X = cell(meas.K,1);  % states
    est.N = zeros(meas.K,1); % cardinality 
    est.G = zeros(meas.K,1); % number of Gaussians
    est.Z = cell(meas.K,1);  % Batched states
    
    model.name  = 'LS GM-PHD';
    model.stream = stream;

    % --- initial prior
    w_update = 1e-16;
    m_update = model.m0 + sqrtm(model.P0)*randn(stream,[model.x_dim,1]);
    P_update = model.P0;

    % --- known list of target identification numbers
    unique_ids = unique(vertcat(meas.ID{:}))'; 

    % --- recursive filtering
    for k = 1:meas.K % time loop
        % 1. Time handling
        model = cfig.time(model, k);
        % --- test to see if inside exposure time
        k_exposure = meas.K; 

        % 2. Predict
        % --- survivors
        [m_predict,P_predict] = predict(cfig,model,m_update,P_update); % surviving components
        w_predict = model.P_S*w_update; % surviving weights

        % 3. Births
        [m_birth, P_birth, w_birth] = gen_gms(model,model.w_birth,model.m_birth,model.P_birth,J_birth); 
        % --- append
        m_predict = cat(2,m_predict,m_birth); 
        P_predict = cat(3,P_predict,P_birth);
        w_predict = cat(1,w_predict,w_birth); 

        % 4. Batch measurement processing
        % --- generate tracklets if there are any
        Zk_cell = cell(length(unique_ids),1);
        k_step = k; 
        while k_step <= k_exposure
            [~,idxs,tidxs] = intersect(meas.ID{k_step},unique_ids);
            if ~isempty(idxs)
                % store measurements with same ID's into a batch
                for ii = 1:length(idxs)
                    Zk_cell{tidxs(ii)} = cat(2, Zk_cell{tidxs(ii)}, [meas.T{k_step}(idxs(ii)); meas.Z{k_step}(:,idxs(ii))]);
                end
                % ensures we only batch these measurements once
                meas.ID{k_step}(idxs)   = []; 
            end
            % iterate to next time step
            if k_step < k_exposure
                k_step = k_step + 1;
            else
                break
            end
        end
        % --- generate batch processed measurements
        if any(~cellfun(@isempty, Zk_cell))
            Zk_ss  = []; Rk_ss  = []; Tk_ss  = [];
            for ii = 1:length(Zk_cell)
                % if a tracklet exists and meets observability requirements
                if ~isempty(Zk_cell{ii}) 
                    % proposal distribution 
                    m_proposal = [model.Time(k); model.m_birth]; 
                    % batch process
                    [Zii_ss, Rii_ss] = batchls(cfig, model, Zk_cell{ii}, m_proposal);
                    % store
                    Zk_ss = cat(3, Zk_ss, Zii_ss);
                    Rk_ss = cat(4, Rk_ss, Rii_ss);
                    Tk_ss = cat(2, Tk_ss, [Zk_cell{ii}(1,1); Zk_cell{ii}(1,end); size(Zk_cell{ii},2)]);
                else
                    continue
                end
            end
            est.Z{k} = Zk_ss; 
        end

        % 5. Update
        % --- missed detections
        w_update = model.Q_D*w_predict;
        m_update = m_predict;
        P_update = P_predict;    
        % --- KF update and compute hypotheses
        if any(~cellfun(@isempty, Zk_cell))
            for ii=1:size(Zk_ss,3) % loop through each measurement
                m_temp_array = []; P_temp_array = []; w_temp_array = [];
                for jj = 1:size(m_predict,2) % loop through each hypothesis
                    m_temp = zeros(model.x_dim,size(Zk_ss,2)); P_temp = zeros(model.x_dim,model.x_dim,size(Zk_ss,2)); U_temp = zeros(size(Zk_ss,2),1);
                    for ll = 1:size(Zk_ss,2) % loop through each measurement Gaussian mixture component
                        [m_temp(:,ll),P_temp(:,:,ll),U_temp(ll)] = update(cfig,model,Zk_ss(:,ll,ii),Rk_ss(:,:,ll,ii),m_predict(:,jj),P_predict(:,:,jj));
                    end
                    % Gaussian sum update
                    U_temp = exp(U_temp);
                    m_temp_bar = m_temp * U_temp ./ sum(U_temp);
                    P_temp_bar = zeros(model.x_dim);
                    for ll = 1:size(Zk_ss,2)
                        P_temp_bar = P_temp_bar + U_temp(ll) * (P_temp(:,:,ll) + (m_temp(:,ll)- m_temp_bar)*(m_temp(:,ll)- m_temp_bar)');
                    end
                    P_temp_bar = P_temp_bar ./ sum(U_temp);
                    w_temp_bar = mean(U_temp);
                    %--- store
                    m_temp_array = cat(2,m_temp_array,m_temp_bar); 
                    P_temp_array = cat(3,P_temp_array,P_temp_bar); 
                    w_temp_array = cat(1,w_temp_array,w_temp_bar); 
                end
                % clutter intensity switching function
                if Tk_ss(3,ii) < ceil(model.x_dim / model.z_dim) 
                    kappa_c = beta_c;
                else
                    kappa_c = model.lambda_c*model.pdf_c/(Tk_ss(3,ii)^(vareps_c*model.z_dim/2));
                end
                % this is a numerically stable way to do this using the log-sum-exp trick
                log_w_temp = log(model.P_D) + log(w_predict(:)) + log(w_temp_array);
                log_sum_w_temp = max(log_w_temp) + log(sum(exp(log_w_temp - max(log_w_temp))));
                log_denom_w_temp = log(kappa_c + exp(log_sum_w_temp));
                log_w_temp = log_w_temp - log_denom_w_temp;
                w_temp = exp(log_w_temp); 

                % append
                w_update = cat(1, w_update, w_temp);
                m_update = cat(2, m_update, m_temp_array);
                P_update = cat(3, P_update, P_temp_array); 
            end
        end

        % 6. Pruning, merging, and capping
        [w_update,m_update,P_update]= gaus_prune(w_update,m_update,P_update,...
            'ElimThreshold',elim_threshold);
        [w_update,m_update,P_update]= gaus_merge(w_update,m_update,P_update,...
            'MergeThreshold',merge_threshold);   
        [w_update,m_update,P_update]= gaus_cap(w_update,m_update,P_update,...
            'Lmax',cap_limit);

        % 7. State extraction [1] (not used by filter)
        idx= find(w_update > 0.5 );
        for ii = 1:length(idx)
            repeat_num_targets= round(w_update(idx(ii)));
            est.X{k} = [ est.X{k} repmat(m_update(:,idx(ii)),[1,repeat_num_targets]) ];
            est.N(k) = est.N(k) + repeat_num_targets;
        end
        est.G(k) = length(w_update);
        
        % 8. Diagnostics
        if ~strcmp(run_flag,'silence')
            disp([' time= ',num2str(k),...
             ' #est mean=' num2str(sum(w_update),4),... % estimated mean number of targets
             ' #est card=' num2str(est.N(k),4)]);       % estimated cardinality 
        end
    end % end time loop
%
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Helper functions
% Predict States and Covariances for the Extended Kalman Filter (EKF)
function [Xkp1,Pkp1] = predict(cfig,model,Xk,Pk)   
% 
    Xkp1 = zeros(size(Xk));
    Pkp1 = zeros(size(Pk));
    for idxp=1:size(Xk,2)
        Xkp1(:,idxp) = cfig.fprop(model, Xk(:,idxp), 'noiseless');
        Phi = expm(cfig.F(Xk(:,idxp), model)*model.T); % approximation for STM for small time steps
        %-- add scaled noise 
        Pkp1(:,:,idxp) = Phi*Pk(:,:,idxp)*Phi' + model.Q;
        Pkp1(:,:,idxp) = (Pkp1(:,:,idxp) + Pkp1(:,:,idxp)')/2;
    end
%
end

% Update States and Covariances for the Extended Kalman Filter (EKF)
function [m_update,P_update,U_update] = update(cfig,model,y,R,m,P)
%      
    U_update = zeros(size(m,2),1);
    m_update = zeros(size(m,1),size(m,2));
    P_update = zeros(size(m,1),size(m,1),size(m,2));
    for jj=1:size(m,2)
        % individual EKF update
        ybar = m(:,jj);
        Hj   = eye(model.x_dim);
        Pxxj = P(:,:,jj);
        Pxyj = Pxxj * Hj';
        Pyyj= Hj*Pxxj*Hj' + R; 
        Pyyj= (Pyyj + Pyyj')/2;   % additional step to avoid numerical problems
        det_Pyyj = prod(eig(Pyyj)); iPyyj = pinv(Pyyj);
        Kj = Pxyj * iPyyj;
        m_update(:,jj)  = m(:,jj) + Kj*(y-ybar);
        Ij = eye(size(m,1));
        P_update(:,:,jj) = (Ij - Kj*Hj) * Pxxj * (Ij - Kj*Hj)' + Kj * R * Kj'; % Joseph form
        
        % weight update
        U_update(jj) = -(y-ybar)' * (iPyyj * (y-ybar)) / 2 ...
                       - log(det_Pyyj) / 2 ...
                       - log(2*pi) * size(y,1) / 2;
    end
%
end