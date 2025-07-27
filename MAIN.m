%% A simple 3D example.
% 08-01-2025 %
%**************************************************************************
% author: dalton durant
% email: ddurant@utexas.edu, thedaltondurant@gmail.com
%
% USE AT YOU OWN DISCRETION
%
% This code is free to use, all I ask is you cite me when reasonable.
%**************************************************************************
% Batch Processing Tracklets in the GM-PHD Filter using Least Squares and 
% Markov Chain Monte Carlo.
%
% SOURCE: 
%
clc; clear all; close all; 
w = warning ('off','all');
%parfevalOnAll(@warning,0,'off','all');
addpath('Filters\');
addpath('Support\');
seed = 42; % random seed value

set(0,'DefaultTextFontSize',8);
set(0,'DefaultAxesFontSize',8);
set(0,'DefaultLegendFontSize',8);
set(0,'DefaultLineLineWidth',2);
set(0,'DefaultAxesLineWidth',1);
set(0,'DefaultLegendInterpreter','Latex');
set(0,'DefaultAxesTickLabelInterpreter','Latex');
set(0,'DefaultTextInterpreter','Latex');

fprintf('Simple multi-target tracking... \n')

options = odeset('reltol',1e-12,'abstol',1e-12); % tolerances for ode solver
cfig = Config();
cfig.options = options;

run('parameters.m')

cfig.nMonte = 1; % number of Monte Carlo simulations
toggle_loadmat.truths      = false; % load truth.mat from /.data
toggle_loadmat.results     = false; % load filter .mat's from /.data
toggle_cleardata           = true;  % overwrite the .mat's in /.data

filters = {
           % GM-PHD
           'RFS\gmphd'; ...                                                % 1.
           % (Batch) LS GM-PHD 
           'RFS\blsgmphd'; ...                                             % 2.
           % (Batch) MCMC GM-PHD 
           'RFS\bmcmcgmphd'; ...                                           % 3.
            };

colors = {% colorblind colors
          '#D81B60'; ... % red
          '#009E73'; ... % green
          '#0072B2'; ... % blue
          };

run_this = [1 2 3]; % run these filters


% . Loading MATs or Generating Truths
truths_fname    = append('.data/','truth','.mat');
if toggle_loadmat.truths
    fprintf('\n Loading Truths... \n')
    %--- loading truth mats
    if isfile(truths_fname)
        load(truths_fname, "truth_STORE", "model_STORE", "meas_STORE")
    end

else
    % --- generating noisy measurements
    fprintf('\n Generating True States and Noisy Measurements... \n')
    truth_STORE    = truth;
    model_STORE    = cell(cfig.nMonte,1);
    meas_STORE     = cell(cfig.nMonte,1);
    %parfor iMonte = 1:cfig.nMonte % run in parallel
    for iMonte = 1:cfig.nMonte
        stream = RandStream('mrg32k3a','Seed', seed + iMonte-1);
        model_STORE{iMonte} = model;
        model_STORE{iMonte}.stream = stream;
        fprintf('\n Monte Carlo Iteration: %d \n', iMonte)
        starttime = tic;
        % --- run sim
        meas_STORE{iMonte}.K  = model.len_time; 
        meas_STORE{iMonte}.Z  = cell(model.len_time,1); % measurements
        meas_STORE{iMonte}.C  = cell(model.len_time,1); % state space representation of clutter measurements (helpful for plotting)
        meas_STORE{iMonte}.T  = cell(model.len_time,1); % times at which measurements occur (for batch processing filters only)
        meas_STORE{iMonte}.ID = cell(model.len_time,1); % target identification numbers (for batch processing filters only)
        for tk = 1:model_STORE{iMonte}.len_time
            Zk = []; Ck = []; Tk = []; IDk = [];
            % generate measurements
            if truth.total_tracks > 0
                idx = find( rand(stream, [truth.total_tracks,1]) <= model.P_D ); % detected target indices
                if ~isempty(idx)
                    z = cfig.h(model_STORE{iMonte}, truth.X{tk}(:,idx), 'noise');
                    rho = z(1,:); az = z(2,:); el = z(3,:);
                    Zk   = cat(2, Zk, z);
                    Ck   = cat(2, Ck, [rho.*cos(az).*cos(el); rho.*sin(az).*cos(el); rho.*sin(el)]); % Cartesian detection points (for plotting)
                    Tk   = cat(1, Tk, repmat(time_array(tk),length(idx),1));
                    IDk  = cat(1, IDk, idx);
                end
            end
            % generate clutter 
            if model.lambda_c > 0
                N_c = poissrnd(model.lambda_c); % number of clutter points
                c = repmat(model.range_c(:,1),[1 N_c]) + diag(model.range_c*[ -1; 1])*rand(3,N_c); % Cartesian clutter points (for plotting)   
                z_c = cfig.h(model,c,'noiseless'); 
                Zk = cat(2, Zk, z_c); 
                Ck   = cat(2, Ck, c); % Cartesian clutter points (for plotting)                                                     
                Tk   = cat(1, Tk, repmat(time_array(tk),1,N_c)');
                IDk  = cat(1, IDk, truth.total_tracks + randperm(stream, 1e5, N_c)');
            end
            % store
            meas_STORE{iMonte}.Z{tk}  = Zk;
            meas_STORE{iMonte}.C{tk}  = Ck;
            meas_STORE{iMonte}.T{tk}  = Tk;
            meas_STORE{iMonte}.ID{tk} = IDk;
        end
    end
    %--- if mat doesn't exist, then make one
    if ~isfile(truths_fname) || toggle_cleardata 
        save(truths_fname, "truth_STORE", "model_STORE", "meas_STORE")
    end
end



% . Loading MATs or Running Filters
Time_STORE      = zeros(length(run_this),1);
est_STORE_all   = cell(cfig.nMonte,length(run_this));
model_STORE_all = cell(cfig.nMonte,length(run_this));
meas_STORE_all  = cell(cfig.nMonte,length(run_this));
if toggle_loadmat.results 
    fprintf('\n Loading Results... \n')
    %--- loading results mats
    for ii = 1:length(run_this)
        [~,D] = fileparts(filters{run_this(ii)});
        results_fname  = append('.data/',D,'.mat');
        if isfile(results_fname)
            load(results_fname, "Time_i", "est_i", "meas_i", "model_i")
            %--- store data
            Time_STORE(ii)              = Time_i;
            est_STORE_all(:,ii)         = est_i(1:cfig.nMonte,:);
            meas_STORE_all(:,ii)        = meas_i(1:cfig.nMonte,:);
            model_STORE_all(:,ii)       = model_i(1:cfig.nMonte,:);
            for iMonte = 1:cfig.nMonte
                if ii <= length(colors)
                    model_STORE_all{iMonte,ii}.color = colors{ii};
                end
            end
        end
    end % end filter loop
else
    fprintf('\n Running Filters... \n')
    %--- running filters
    for ii = 1:length(run_this)
        addpath(append('Filters\',filters{run_this(ii)}))
        [~,D] = fileparts(filters{run_this(ii)});
        results_fname  = append('.data/',D,'.mat');
    
        fprintf('\n Filter: '); disp([filters{run_this(ii)}]); fprintf(' \n');
        stoptime = zeros(cfig.nMonte,1);
        %parfor iMonte = 1:cfig.nMonte % run in parallel
        for iMonte = 1:cfig.nMonte
            stream = RandStream('mrg32k3a','Seed', seed + iMonte-1);
            fprintf(' . \n . \n . \n Monte Carlo Iteration: %d \n', iMonte)
            starttime = tic;  
            [model_STORE_all{iMonte,ii},meas_STORE_all{iMonte,ii},est_STORE_all{iMonte,ii}] = run_filter(stream,cfig,model_STORE{iMonte},meas_STORE{iMonte});
            stoptime(iMonte) = toc(starttime); % wall-clock-time [s]
            if ii <= length(colors)
                model_STORE_all{iMonte,ii}.color = colors{ii};
            end
        end
        Time_i = mean(stoptime);
        %--- store data
        Time_STORE(ii) = Time_i;
        est_i          = est_STORE_all(:,ii);
        meas_i         = meas_STORE_all(:,ii);
        model_i        = model_STORE_all(:,ii);

        %--- if mat doesn't exist, then make one
        if ~isfile(results_fname) || toggle_cleardata       
            save(results_fname, "Time_i", "est_i", "meas_i", "model_i")
        end
    
        %--- report results
        fprintf('\n Report of: '); disp([filters{run_this(ii)}]); fprintf(' \n');
        fprintf('    Wall-Clock-Time: %d [s] \n', mean(stoptime)) 
    
        %--- reset things for next filter
        rmpath(append('Filters\',filters{run_this(ii)}))
    end % end filter loop
end




% . Plot Results
fprintf('\n Plotting... \n')
plot_results_MC(cfig,model_STORE_all,truth_STORE,meas_STORE_all,est_STORE_all,Time_STORE);