clc; clear;

load('nature.mat');
load('BC_2000_2017_sa_dt.mat');
load('BCC_2000_2017_sa_dt.mat');
load('FO_2000_2017_sa_dt.mat');
load('S_2000_2017_sa_dt.mat');
load('FO_no_W_2000_2017_sa_dt.mat');
load('FO_no_tau_eta_2000_2017_sa_dt.mat');

% initialize variables
num_iterations = 50000;
states = BC_2000_2017_sa_dt.Properties.VariableNames(3:end-1);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Compute mutual information %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% define symbolized time series
X1 = double(diff(nature.Background_checks(ismember(nature.Year,[2000:2017])))>0);
X2 = double(diff(BCC_2000_2017_sa_dt.USA)>0);
X3 = double(diff(S_2000_2017_sa_dt.USA)>0);
X4 = double(diff(FO_2000_2017_sa_dt.USA)>0);
Y = nature.Mass_shooting(ismember(nature.Year,[2000:2017]));
Y = double(Y(1:end-1)>0);
Z = double(diff(nature.Firearm_laws_and_regulations(ismember(nature.Year,[2000:2017])))>0);


% test time series' dependence on their previous time step
[I_X1,~,~,p_value_X1] = test_past_independence(X1,50000); % 
[I_X2,~,~,p_value_X2] = test_past_independence(X2,50000); % 
[I_X3,~,~,p_value_X3] = test_past_independence(X3,50000); % 
[I_X4,~,~,p_value_X4] = test_past_independence(X4,50000); % 
[I_Y,~,~,p_value_Y] = test_past_independence(Y,50000); % 
[I_Z,~,~,p_value_Z] = test_past_independence(Z,50000); % 

% % test for contemporanoues effects in BC/MS/MO triad
% [I_X1_Y_Z,~,~,p_value_X1_Y_Z] = test_contemporaneous(X1,Y,Z,50000); % BC and MS, conditioned on MO
% [I_X1_Z_Y,~,~,p_value_X1_Z_Y] = test_contemporaneous(X1,Z,Y,50000); % BC and MO, conditioned on MS
% [I_Y_Z_X1,~,~,p_value_Y_Z_X1] = test_contemporaneous(Y,Z,X1,50000); % MS and MO, conditioned on BC
% 
% % test for contemporanoues effects in BCC/MS/MO triad
% [I_X2_Y_Z,~,~,p_value_X2_Y_Z] = test_contemporaneous(X2,Y,Z,50000); % BCC and MS, conditioned on MO
% [I_X2_Z_Y,~,~,p_value_X2_Z_Y] = test_contemporaneous(X2,Z,Y,50000); % BCC and MO, conditioned on MS
% [I_Y_Z_X2,~,~,p_value_Y_Z_X2] = test_contemporaneous(Y,Z,X2,50000); % MS and MO, conditioned on BCC
% 
% % test for contemporanoues effects in SF/MS/MO triad
% [I_X3_Y_Z,~,~,p_value_X3_Y_Z] = test_contemporaneous(X3,Y,Z,50000); % SF and MS, conditioned on MO
% [I_X3_Z_Y,~,~,p_value_X3_Z_Y] = test_contemporaneous(X3,Z,Y,50000); % SF and MO, conditioned on MS
% [I_Y_Z_X3,~,~,p_value_Y_Z_X3] = test_contemporaneous(Y,Z,X3,50000); % MS and MO, conditioned on SF
% 
% % test for contemporanoues effects in FO/MS/MO triad
% [I_X4_Y_Z,~,~,p_value_X4_Y_Z] = test_contemporaneous(X4,Y,Z,50000); % FO and MS, conditioned on MO
% [I_X4_Z_Y,~,~,p_value_X4_Z_Y] = test_contemporaneous(X4,Z,Y,50000); % FO and MO, conditioned on MS
% [I_Y_Z_X4,~,~,p_value_Y_Z_X4] = test_contemporaneous(Y,Z,X4,50000); % MS and MO, conditioned on FO



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Background checks as a proxy %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% national level %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% define time series
X = nature.Background_checks(ismember(nature.Year,[2000:2017]));
Y = nature.Mass_shooting(ismember(nature.Year,[2000:2017]));
Z = nature.Firearm_laws_and_regulations(ismember(nature.Year,[2000:2017]));

% symbolize data
XYZ = [diff(X)>0 Y(1:end-1)>0 diff(Z)>0];
XYZX1Y1Z1 = array2table([XYZ(1:end-1,:) XYZ(2:end,:)],'VariableNames',{'X','Y','Z','X1','Y1','Z1'}); % time series for t and t+1

% compute te, surrogate distributions,and 95th percentile
[TE,distributions,percentiles_95,p_values] = compute_TE(XYZX1Y1Z1,num_iterations);
TE_results = array2table([TE; p_values; percentiles_95],'VariableNames',{'TE_X_YZ','TE_X_ZY','TE_Y_XZ','TE_Y_ZX','TE_Z_XY','TE_Z_YX'},'RowNames', {'TE','p_values','95th_percentile'});



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% state level %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

state_TE = array2table(nan(length(states),7),'VariableNames',{'state','TE_X_YZ','TE_X_ZY','TE_Y_XZ','TE_Y_ZX','TE_Z_XY','TE_Z_YX',});
state_percentiles_95 = array2table(nan(length(states),7),'VariableNames',{'state','TE_X_YZ','TE_X_ZY','TE_Y_XZ','TE_Y_ZX','TE_Z_XY','TE_Z_YX',});
state_p_values = array2table(nan(length(states),7),'VariableNames',{'state','TE_X_YZ','TE_X_ZY','TE_Y_XZ','TE_Y_ZX','TE_Z_XY','TE_Z_YX',});
state_significance = array2table(nan(length(states),7),'VariableNames',{'state','TE_X_YZ','TE_X_ZY','TE_Y_XZ','TE_Y_ZX','TE_Z_XY','TE_Z_YX',});
for s = states
    
    % define time series
    X = BC_2000_2017_sa_dt{:,strcmp(BC_2000_2017_sa_dt.Properties.VariableNames,s)};
    Y = nature.Mass_shooting(ismember(nature.Year,[2000:2017]));
    Z = nature.Firearm_laws_and_regulations(ismember(nature.Year,[2000:2017]));

    % symbolize data
    XYZ = [diff(X)>0 Y(1:end-1)>0 diff(Z)>0];
    XYZX1Y1Z1 = array2table([XYZ(1:end-1,:) XYZ(2:end,:)],'VariableNames',{'X','Y','Z','X1','Y1','Z1'}); % time series for t and t+1

    % compute te, surrogate distributions,and 95th percentile
    [TE,distributions,percentiles_95,p_values] = compute_TE(XYZX1Y1Z1,num_iterations);
 
    % summarize results
    state_TE(strcmp(states,s),2:end) = array2table(TE);
    state_percentiles_95(strcmp(states,s),2:end) = array2table(percentiles_95);
    state_p_values(strcmp(states,s),2:end) = array2table(p_values);
    
end
state_TE.state = transpose(states);
state_percentiles_95.state = transpose(states);
state_p_values.state = transpose(states);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Background checks per capita as a proxy %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% national level %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% define time series for 2000-2017
X = BCC_2000_2017_sa_dt.USA;
Y = nature.Mass_shooting(ismember(nature.Year,[2000:2017]));
Z = nature.Firearm_laws_and_regulations(ismember(nature.Year,[2000:2017]));

% symbolize data
XYZ = [diff(X)>0 Y(1:end-1)>0 diff(Z)>0];
XYZX1Y1Z1 = array2table([XYZ(1:end-1,:) XYZ(2:end,:)],'VariableNames',{'X','Y','Z','X1','Y1','Z1'}); % time series for t and t+1

% compute te, surrogate distributions,and 95th percentile
[TE,distributions,percentiles_95,p_values] = compute_TE(XYZX1Y1Z1,num_iterations);
TE_results = array2table([TE; p_values; percentiles_95],'VariableNames',{'TE_X_YZ','TE_X_ZY','TE_Y_XZ','TE_Y_ZX','TE_Z_XY','TE_Z_YX'},'RowNames', {'TE','p_values','95th_percentile'});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% state level %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% 2000-2017
state_TE = array2table(nan(length(states),7),'VariableNames',{'state','TE_X_YZ','TE_X_ZY','TE_Y_XZ','TE_Y_ZX','TE_Z_XY','TE_Z_YX',});
state_percentiles_95 = array2table(nan(length(states),7),'VariableNames',{'state','TE_X_YZ','TE_X_ZY','TE_Y_XZ','TE_Y_ZX','TE_Z_XY','TE_Z_YX',});
state_p_values = array2table(nan(length(states),7),'VariableNames',{'state','TE_X_YZ','TE_X_ZY','TE_Y_XZ','TE_Y_ZX','TE_Z_XY','TE_Z_YX',});
state_significance = array2table(nan(length(states),7),'VariableNames',{'state','TE_X_YZ','TE_X_ZY','TE_Y_XZ','TE_Y_ZX','TE_Z_XY','TE_Z_YX',});
for s = states
    
    % define time series
    X = BCC_2000_2017_sa_dt{:,strcmp(BCC_2000_2017_sa_dt.Properties.VariableNames,s)};
    Y = nature.Mass_shooting(ismember(nature.Year,[2000:2017]));
    Z = nature.Firearm_laws_and_regulations(ismember(nature.Year,[2000:2017]));

    % symbolize data
    XYZ = [diff(X)>0 Y(1:end-1)>0 diff(Z)>0];
    XYZX1Y1Z1 = array2table([XYZ(1:end-1,:) XYZ(2:end,:)],'VariableNames',{'X','Y','Z','X1','Y1','Z1'}); % time series for t and t+1

    % compute te, surrogate distributions,and 95th percentile
    [TE,distributions,percentiles_95,p_values] = compute_TE(XYZX1Y1Z1,num_iterations);
    
    % summarize results
    state_TE(strcmp(states,s),2:end) = array2table(TE);
    state_percentiles_95(strcmp(states,s),2:end) = array2table(percentiles_95);
    state_p_values(strcmp(states,s),2:end) = array2table(p_values);
 
end
state_TE.state = transpose(states);
state_percentiles_95.state = transpose(states);
state_p_values.state = transpose(states);



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Suicides as a proxy %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% national level %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% define time series
X = S_2000_2017_sa_dt.USA;
Y = nature.Mass_shooting(ismember(nature.Year,[2000:2017]));
Z = nature.Firearm_laws_and_regulations(ismember(nature.Year,[2000:2017]));

% symbolize data
XYZ = [diff(X)>0 Y(1:end-1)>0 diff(Z)>0];
XYZX1Y1Z1 = array2table([XYZ(1:end-1,:) XYZ(2:end,:)],'VariableNames',{'X','Y','Z','X1','Y1','Z1'}); % time series for t and t+1

% compute te, surrogate distributions,and 95th percentile
[TE,distributions,percentiles_95,p_values] = compute_TE(XYZX1Y1Z1,num_iterations);
TE_results = array2table([TE; p_values; percentiles_95],'VariableNames',{'TE_X_YZ','TE_X_ZY','TE_Y_XZ','TE_Y_ZX','TE_Z_XY','TE_Z_YX'},'RowNames', {'TE','p_values','95th_percentile'});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% state level %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

state_TE = array2table(nan(length(states),7),'VariableNames',{'state','TE_X_YZ','TE_X_ZY','TE_Y_XZ','TE_Y_ZX','TE_Z_XY','TE_Z_YX',});
state_percentiles_95 = array2table(nan(length(states),7),'VariableNames',{'state','TE_X_YZ','TE_X_ZY','TE_Y_XZ','TE_Y_ZX','TE_Z_XY','TE_Z_YX',});
state_p_values = array2table(nan(length(states),7),'VariableNames',{'state','TE_X_YZ','TE_X_ZY','TE_Y_XZ','TE_Y_ZX','TE_Z_XY','TE_Z_YX',});
state_significance = array2table(nan(length(states),7),'VariableNames',{'state','TE_X_YZ','TE_X_ZY','TE_Y_XZ','TE_Y_ZX','TE_Z_XY','TE_Z_YX',});
for s = states
    
    % define time series
    X = S_2000_2017_sa_dt{:,strcmp(S_2000_2017_sa_dt.Properties.VariableNames,s)};
    Y = nature.Mass_shooting(ismember(nature.Year,[2000:2017]));
    Z = nature.Firearm_laws_and_regulations(ismember(nature.Year,[2000:2017]));

    % symbolize data
    XYZ = [diff(X)>0 Y(1:end-1)>0 diff(Z)>0];
    XYZX1Y1Z1 = array2table([XYZ(1:end-1,:) XYZ(2:end,:)],'VariableNames',{'X','Y','Z','X1','Y1','Z1'}); % time series for t and t+1

    % compute te, surrogate distributions,and 95th percentile
    [TE,distributions,percentiles_95,p_values] = compute_TE(XYZX1Y1Z1,num_iterations);
    
    % organize and save results
    state_TE(strcmp(states,s),2:end) = array2table(TE);
    state_percentiles_95(strcmp(states,s),2:end) = array2table(percentiles_95);
    state_p_values(strcmp(states,s),2:end) = array2table(p_values);
   
end
state_TE.state = transpose(states);
state_percentiles_95.state = transpose(states);
state_p_values.state = transpose(states);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Gun ownership as a proxy %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% national level %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% define time series for 2000-2017
X = FO_2000_2017_sa_dt.USA;
Y = nature.Mass_shooting(ismember(nature.Year,[2000:2017]));
Z = nature.Firearm_laws_and_regulations(ismember(nature.Year,[2000:2017]));

% symbolize data
XYZ = [diff(X)>0 Y(1:end-1)>0 diff(Z)>0];
XYZX1Y1Z1 = array2table([XYZ(1:end-1,:) XYZ(2:end,:)],'VariableNames',{'X','Y','Z','X1','Y1','Z1'}); % time series for t and t+1
XYZX1Y1Z1 = array2table([XYZ(1:end-3,:) XYZ.X(4:end,:) XYZ.Y(2:end-2,:) XYZ.Z(4:end,:)],'VariableNames',{'X','Y','Z','X1','Y1','Z1'}); % time series for t and t+1

% compute te, surrogate distributions,and 95th percentile
[TE,distributions,percentiles_95,p_values] = compute_TE(XYZX1Y1Z1,num_iterations);
TE_results = array2table([TE; p_values; percentiles_95],'VariableNames',{'TE_X_YZ','TE_X_ZY','TE_Y_XZ','TE_Y_ZX','TE_Z_XY','TE_Z_YX'},'RowNames', {'TE','p_values','95th_percentile'});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% state level %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

state_TE = array2table(nan(length(states),7),'VariableNames',{'state','TE_X_YZ','TE_X_ZY','TE_Y_XZ','TE_Y_ZX','TE_Z_XY','TE_Z_YX',});
state_percentiles_95 = array2table(nan(length(states),7),'VariableNames',{'state','TE_X_YZ','TE_X_ZY','TE_Y_XZ','TE_Y_ZX','TE_Z_XY','TE_Z_YX',});
state_p_values = array2table(nan(length(states),7),'VariableNames',{'state','TE_X_YZ','TE_X_ZY','TE_Y_XZ','TE_Y_ZX','TE_Z_XY','TE_Z_YX',});
state_significance = array2table(nan(length(states),7),'VariableNames',{'state','TE_X_YZ','TE_X_ZY','TE_Y_XZ','TE_Y_ZX','TE_Z_XY','TE_Z_YX',});
for s = states
    
    % define time series
    X = FO_2000_2017_sa_dt{:,strcmp(FO_2000_2017_sa_dt.Properties.VariableNames,s)};
    Y = nature.Mass_shooting(ismember(nature.Year,[2000:2017]));
    Z = nature.Firearm_laws_and_regulations(ismember(nature.Year,[2000:2017]));

    % symbolize data
    XYZ = [diff(X)>0 Y(1:end-1)>0 diff(Z)>0];
    XYZX1Y1Z1 = array2table([XYZ(1:end-1,:) XYZ(2:end,:)],'VariableNames',{'X','Y','Z','X1','Y1','Z1'}); % time series for t and t+1

    % compute te, surrogate distributions,and 95th percentile
    [TE,distributions,percentiles_95,p_values] = compute_TE(XYZX1Y1Z1,num_iterations);
    
    % summarize and save results
    state_TE(strcmp(states,s),2:end) = array2table(TE);
    state_percentiles_95(strcmp(states,s),2:end) = array2table(percentiles_95);
    state_p_values(strcmp(states,s),2:end) = array2table(p_values);

end
state_TE.state = transpose(states);
state_percentiles_95.state = transpose(states);
state_p_values.state = transpose(states);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Firearm ownership without spatial interactions (W=0) as a proxy %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% national level %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% define time series for 2000-2017
X = FO_no_W_2000_2017_sa_dt.USA;
Y = nature.Mass_shooting(ismember(nature.Year,[2000:2017]));
Z = nature.Firearm_laws_and_regulations(ismember(nature.Year,[2000:2017]));

% symbolize data
XYZ = [diff(X)>0 Y(1:end-1)>0 diff(Z)>0];
XYZX1Y1Z1 = array2table([XYZ(1:end-1,:) XYZ(2:end,:)],'VariableNames',{'X','Y','Z','X1','Y1','Z1'}); % time series for t and t+1

% compute te, surrogate distributions,and 95th percentile
[TE,distributions,percentiles_95,p_values] = compute_TE(XYZX1Y1Z1,num_iterations);
TE_results = array2table([TE; p_values; percentiles_95],'VariableNames',{'TE_X_YZ','TE_X_ZY','TE_Y_XZ','TE_Y_ZX','TE_Z_XY','TE_Z_YX'},'RowNames', {'TE','p_values','95th_percentile'});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% state level %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

state_TE = array2table(nan(length(states),7),'VariableNames',{'state','TE_X_YZ','TE_X_ZY','TE_Y_XZ','TE_Y_ZX','TE_Z_XY','TE_Z_YX',});
state_percentiles_95 = array2table(nan(length(states),7),'VariableNames',{'state','TE_X_YZ','TE_X_ZY','TE_Y_XZ','TE_Y_ZX','TE_Z_XY','TE_Z_YX',});
state_p_values = array2table(nan(length(states),7),'VariableNames',{'state','TE_X_YZ','TE_X_ZY','TE_Y_XZ','TE_Y_ZX','TE_Z_XY','TE_Z_YX',});
state_significance = array2table(nan(length(states),7),'VariableNames',{'state','TE_X_YZ','TE_X_ZY','TE_Y_XZ','TE_Y_ZX','TE_Z_XY','TE_Z_YX',});
for s = states
    
    % define time series
    X = FO_2000_2017_sa_dt{:,strcmp(FO_2000_2017_sa_dt.Properties.VariableNames,s)};
    Y = nature.Mass_shooting(ismember(nature.Year,[2000:2017]));
    Z = nature.Firearm_laws_and_regulations(ismember(nature.Year,[2000:2017]));

    % symbolize data
    XYZ = [diff(X)>0 Y(1:end-1)>0 diff(Z)>0];
    XYZX1Y1Z1 = array2table([XYZ(1:end-1,:) XYZ(2:end,:)],'VariableNames',{'X','Y','Z','X1','Y1','Z1'}); % time series for t and t+1

    % compute te, surrogate distributions,and 95th percentile
    [TE,distributions,percentiles_95,p_values] = compute_TE(XYZX1Y1Z1,num_iterations);
    
    % summarize and save results
    state_TE(strcmp(states,s),2:end) = array2table(TE);
    state_percentiles_95(strcmp(states,s),2:end) = array2table(percentiles_95);
    state_p_values(strcmp(states,s),2:end) = array2table(p_values);

end
state_TE.state = transpose(states);
state_percentiles_95.state = transpose(states);
state_p_values.state = transpose(states);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Firearm ownership estimated without tau and eta as a proxy %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% national level %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% define time series for 2000-2017
X = FO_no_tau_eta_2000_2017_sa_dt.USA;
Y = nature.Mass_shooting(ismember(nature.Year,[2000:2017]));
Z = nature.Firearm_laws_and_regulations(ismember(nature.Year,[2000:2017]));

% symbolize data
XYZ = [diff(X)>0 Y(1:end-1)>0 diff(Z)>0];
XYZX1Y1Z1 = array2table([XYZ(1:end-1,:) XYZ(2:end,:)],'VariableNames',{'X','Y','Z','X1','Y1','Z1'}); % time series for t and t+1

% compute te, surrogate distributions,and 95th percentile
[TE,distributions,percentiles_95,p_values] = compute_TE(XYZX1Y1Z1,num_iterations);
TE_results = array2table([TE; p_values; percentiles_95],'VariableNames',{'TE_X_YZ','TE_X_ZY','TE_Y_XZ','TE_Y_ZX','TE_Z_XY','TE_Z_YX'},'RowNames', {'TE','p_values','95th_percentile'});


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% state level %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

state_TE = array2table(nan(length(states),7),'VariableNames',{'state','TE_X_YZ','TE_X_ZY','TE_Y_XZ','TE_Y_ZX','TE_Z_XY','TE_Z_YX',});
state_percentiles_95 = array2table(nan(length(states),7),'VariableNames',{'state','TE_X_YZ','TE_X_ZY','TE_Y_XZ','TE_Y_ZX','TE_Z_XY','TE_Z_YX',});
state_p_values = array2table(nan(length(states),7),'VariableNames',{'state','TE_X_YZ','TE_X_ZY','TE_Y_XZ','TE_Y_ZX','TE_Z_XY','TE_Z_YX',});
state_significance = array2table(nan(length(states),7),'VariableNames',{'state','TE_X_YZ','TE_X_ZY','TE_Y_XZ','TE_Y_ZX','TE_Z_XY','TE_Z_YX',});
for s = states
    
    % define time series
    X = FO_2000_2017_sa_dt{:,strcmp(FO_2000_2017_sa_dt.Properties.VariableNames,s)};
    Y = nature.Mass_shooting(ismember(nature.Year,[2000:2017]));
    Z = nature.Firearm_laws_and_regulations(ismember(nature.Year,[2000:2017]));

    % symbolize data
    XYZ = [diff(X)>0 Y(1:end-1)>0 diff(Z)>0];
    XYZX1Y1Z1 = array2table([XYZ(1:end-1,:) XYZ(2:end,:)],'VariableNames',{'X','Y','Z','X1','Y1','Z1'}); % time series for t and t+1

    % compute te, surrogate distributions,and 95th percentile
    [TE,distributions,percentiles_95,p_values] = compute_TE(XYZX1Y1Z1,num_iterations);
    
    % summarize and save results
    state_TE(strcmp(states,s),2:end) = array2table(TE);
    state_percentiles_95(strcmp(states,s),2:end) = array2table(percentiles_95);
    state_p_values(strcmp(states,s),2:end) = array2table(p_values);

end
state_TE.state = transpose(states);
state_percentiles_95.state = transpose(states);
state_p_values.state = transpose(states);
