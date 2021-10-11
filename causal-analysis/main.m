%% Transfer entropy analysis on a national level

clc; clear;

load('nature.mat');
load('BC_2000_2017_sa_dt.mat');
load('BCC_2000_2017_sa_dt.mat');
load('FO_2000_2017_sa_dt.mat');
load('S_2000_2017_sa_dt.mat');

% initialize variables
num_iterations = 50000;
states = BC_2000_2017_sa_dt.Properties.VariableNames(3:end-1);


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
[TE,distributions,percentiles_95,p_values,~] = compute_TE(XYZX1Y1Z1,num_iterations);
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
    [TE,distributions,percentiles_95,p_values,~] = compute_TE(XYZX1Y1Z1,num_iterations);
 
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
[TE,distributions,percentiles_95,p_values,~] = compute_TE(XYZX1Y1Z1,num_iterations);
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
    [TE,distributions,percentiles_95,p_values,~] = compute_TE(XYZX1Y1Z1,num_iterations);
    
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
[TE,distributions,percentiles_95,p_values,significance] = compute_TE(XYZX1Y1Z1,num_iterations);
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
    [TE,distributions,percentiles_95,p_values,~] = compute_TE(XYZX1Y1Z1,num_iterations);
    
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

% compute te, surrogate distributions,and 95th percentile
[TE,distributions,percentiles_95,p_values,significance] = compute_TE(XYZX1Y1Z1,num_iterations);
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
    [TE,distributions,percentiles_95,p_values,~] = compute_TE(XYZX1Y1Z1,num_iterations);
    
    % summarize and save results
    state_TE(strcmp(states,s),2:end) = array2table(TE);
    state_percentiles_95(strcmp(states,s),2:end) = array2table(percentiles_95);
    state_p_values(strcmp(states,s),2:end) = array2table(p_values);

end
state_TE.state = transpose(states);
state_percentiles_95.state = transpose(states);
state_p_values.state = transpose(states);

