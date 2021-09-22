%% initialize variables

num_iterations = 50000;
results_directory = ''; % input where data should be saved


load('nature.mat');
load('restrictiveness.mat');
load('BC_2000_2017_sa_dt.mat');
load('BCC_2000_2017_sa_dt.mat');
load('GO_2000_2017_sa_dt.mat');
load('S_2000_2017_sa_dt.mat');

%% Test for stationarity with raw data

display('before processing')
for s = transpose(unique(states))
    
    [h,p0] = adftest(gun_ownership_monthly{1:216,strcmp(gun_ownership_monthly.Properties.VariableNames,s)});
    [h,p1] = adftest(BC_state_monthly{3:218,strcmp(BC_state_monthly.Properties.VariableNames,s)});
    [h,p2] = adftest(BCC_state_monthly{3:218,strcmp(BCC_state_monthly.Properties.VariableNames,s)});
    [h,p3] = adftest(suicides_monthly{3:218,strcmp(suicides_monthly.Properties.VariableNames,s)});
    display(strcat(string(s),{' '},string(p0)))%,{' '},string(p1),{' '},string(p2),{' '},string(p3)))

end

display('after processing')
for s = transpose(unique(states))

    [h,p1] = adftest(GO_2000_2017_sa_dt{:,strcmp(GO_2000_2017_sa_dt.Properties.VariableNames,s)});
    [h,p2] = adftest(BCC_2000_2017_sa_dt{:,strcmp(BCC_2000_2017_sa_dt.Properties.VariableNames,s)});
    [h,p3] = adftest(S_2000_2017_sa_dt{:,strcmp(S_2000_2017_sa_dt.Properties.VariableNames,s)});

    display(strcat(string(s),{' '},string(p1),{' '},string(p2),{' '},string(p3)))

end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Background checks as a proxy %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
display('Background checks as proxy');
results_directory_BC = strcat(results_directory,'/background_checks');

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

% summarize and save results
TE_results = array2table([TE; p_values; percentiles_95],'VariableNames',{'TE_X_YZ','TE_X_ZY','TE_Y_XZ','TE_Y_ZX','TE_Z_XY','TE_Z_YX'},'RowNames', {'TE','p_values','95th_percentile'});
writetable(TE_results,strcat(results_directory_BC,'/national_te.csv'),'WriteRowNames',true);
writetable(distributions,strcat(results_directory_BC,'/national_distributions.csv'));

display('analysis on national level completed');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% state level %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

state_TE = array2table(nan(length(states),7),'VariableNames',{'state','TE_X_YZ','TE_X_ZY','TE_Y_XZ','TE_Y_ZX','TE_Z_XY','TE_Z_YX',});
state_percentiles_95 = array2table(nan(length(states),7),'VariableNames',{'state','TE_X_YZ','TE_X_ZY','TE_Y_XZ','TE_Y_ZX','TE_Z_XY','TE_Z_YX',});
state_p_values = array2table(nan(length(states),7),'VariableNames',{'state','TE_X_YZ','TE_X_ZY','TE_Y_XZ','TE_Y_ZX','TE_Z_XY','TE_Z_YX',});
state_significance = array2table(nan(length(states),7),'VariableNames',{'state','TE_X_YZ','TE_X_ZY','TE_Y_XZ','TE_Y_ZX','TE_Z_XY','TE_Z_YX',});
for s = transpose(states)
    
    % define time series
    X = BC_1999_2017_sa_dt{ismember(BC_1999_2017_sa_dt.Year,[2000:2017]),strcmp(BC_1999_2017_sa_dt.Properties.VariableNames,s)};
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
 
    writetable(distributions,strcat(results_directory_BC,'/distributions_',string(s),'.csv'));
    writetable(state_TE,strcat(results_directory_BC,'/state_te.csv'));
    writetable(state_percentiles_95,strcat(results_directory_BC,'/state_percentiles_95.csv'));
    writetable(state_p_values,strcat(results_directory_BC,'/state_p_values.csv'));
    
    display(string(s));
end
state_TE.state = states;
state_percentiles_95.state = states;
state_p_values.state = states;

writetable(state_TE,strcat(results_directory_BC,'/state_te.csv'));
writetable(state_percentiles_95,strcat(results_directory_BC,'/state_percentiles_95.csv'));
writetable(state_p_values,strcat(results_directory_BC,'/state_p_values.csv'));
    
% % plot results
% results_directory = strcat(results_directory,'/');
% plot_state_TE(state_data,states,restrictiveness,results_directory)


display('analysis on state level completed');


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Background checks per capita as a proxy %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
display('Background checks per capita as proxy');
results_directory_BCC = strcat(results_directory,'/background_checks_per_capita');

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

% summarize and save results
TE_results = array2table([TE; p_values; percentiles_95],'VariableNames',{'TE_X_YZ','TE_X_ZY','TE_Y_XZ','TE_Y_ZX','TE_Z_XY','TE_Z_YX'},'RowNames', {'TE','p_values','95th_percentile'});
writetable(TE_results,strcat(results_directory_BCC,'/national_te.csv'),'WriteRowNames',true);
writetable(distributions,strcat(results_directory_BCC,'/national_distributions.csv'));

display('analysis on national level completed');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% state level %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% 2000-2017
state_TE = array2table(nan(length(states),7),'VariableNames',{'state','TE_X_YZ','TE_X_ZY','TE_Y_XZ','TE_Y_ZX','TE_Z_XY','TE_Z_YX',});
state_percentiles_95 = array2table(nan(length(states),7),'VariableNames',{'state','TE_X_YZ','TE_X_ZY','TE_Y_XZ','TE_Y_ZX','TE_Z_XY','TE_Z_YX',});
state_p_values = array2table(nan(length(states),7),'VariableNames',{'state','TE_X_YZ','TE_X_ZY','TE_Y_XZ','TE_Y_ZX','TE_Z_XY','TE_Z_YX',});
state_significance = array2table(nan(length(states),7),'VariableNames',{'state','TE_X_YZ','TE_X_ZY','TE_Y_XZ','TE_Y_ZX','TE_Z_XY','TE_Z_YX',});
for s = transpose(states)
    
    % define time series
    X = BCC_2000_2017_sa_dt{:,strcmp(BCC_2000_2017_sa_dt.Properties.VariableNames,s)};
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
 
    writetable(distributions,strcat(results_directory_BCC,'/distributions_',string(s),'.csv'));
    writetable(state_TE,strcat(results_directory_BCC,'/state_te.csv'));
    writetable(state_percentiles_95,strcat(results_directory_BCC,'/state_percentiles_95.csv'));
    writetable(state_p_values,strcat(results_directory_BCC,'/state_p_values.csv'));
    
    display(string(s));
 end
state_TE.state = states;
state_percentiles_95.state = states;
state_p_values.state = states;

writetable(state_TE,strcat(results_directory_BCC,'/state_te.csv'));
writetable(state_percentiles_95,strcat(results_directory_BCC,'/state_percentiles_95.csv'));
writetable(state_p_values,strcat(results_directory_BCC,'/state_p_values.csv'));

% % plot results
% results_directory = strcat(results_directory,'/');
% plot_state_TE(state_data,states,restrictiveness,results_directory)


display('analysis on state level completed');


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Suicides as a proxy %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
display('suicides as proxy');
results_directory_S = strcat(results_directory,'/suicides');

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

% summarize and save results
TE_results = array2table([TE; p_values; percentiles_95],'VariableNames',{'TE_X_YZ','TE_X_ZY','TE_Y_XZ','TE_Y_ZX','TE_Z_XY','TE_Z_YX'},'RowNames', {'TE','p_values','95th_percentile'});
writetable(TE_results,strcat(results_directory_S,'/national_te.csv'),'WriteRowNames',true);
writetable(distributions,strcat(results_directory_S,'/national_distributions.csv'));


display('analysis on national level completed');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% state level %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

state_TE = array2table(nan(length(states),7),'VariableNames',{'state','TE_X_YZ','TE_X_ZY','TE_Y_XZ','TE_Y_ZX','TE_Z_XY','TE_Z_YX',});
state_percentiles_95 = array2table(nan(length(states),7),'VariableNames',{'state','TE_X_YZ','TE_X_ZY','TE_Y_XZ','TE_Y_ZX','TE_Z_XY','TE_Z_YX',});
state_p_values = array2table(nan(length(states),7),'VariableNames',{'state','TE_X_YZ','TE_X_ZY','TE_Y_XZ','TE_Y_ZX','TE_Z_XY','TE_Z_YX',});
state_significance = array2table(nan(length(states),7),'VariableNames',{'state','TE_X_YZ','TE_X_ZY','TE_Y_XZ','TE_Y_ZX','TE_Z_XY','TE_Z_YX',});
for s = transpose(states)
    
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
   
    writetable(distributions,strcat(results_directory_S,'/distributions_',string(s),'.csv'));
    writetable(state_TE,strcat(results_directory_S,'/state_te.csv'));
    writetable(state_percentiles_95,strcat(results_directory_S,'/state_percentiles_95.csv'));
    writetable(state_p_values,strcat(results_directory_S,'/state_p_values.csv'));
    
    display(string(s));
end
state_TE.state = states;
state_percentiles_95.state = states;
state_p_values.state = states;

% save results
writetable(state_TE,strcat(results_directory_S,'/state_te.csv'));
writetable(state_percentiles_95,strcat(results_directory_S,'/state_percentiles_95.csv'));
writetable(state_p_values,strcat(results_directory_S,'/state_p_values.csv'));

% % plot results
% results_directory = strcat(results_directory,'suicides_2000-2017/');
% plot_state_TE(state_data,states,restrictiveness,results_directory)


display('analysis on state level completed');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Gun ownership as a proxy %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
display('gun ownership as proxy');
results_directory_GO = strcat(results_directory,'/gun_ownership');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% national level %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% define time series for 2000-2017
X = GO_2000_2017_sa_dt.USA;
Y = nature.Mass_shooting(ismember(nature.Year,[2000:2017]));
Z = nature.Firearm_laws_and_regulations(ismember(nature.Year,[2000:2017]));

% symbolize data
XYZ = [diff(X)>0 Y(1:end-1)>0 diff(Z)>0];
XYZX1Y1Z1 = array2table([XYZ(1:end-1,:) XYZ(2:end,:)],'VariableNames',{'X','Y','Z','X1','Y1','Z1'}); % time series for t and t+1

% compute te, surrogate distributions,and 95th percentile
[TE,distributions,percentiles_95,p_values,significance] = compute_TE(XYZX1Y1Z1,num_iterations);

% summarize and save results
TE_results = array2table([TE; p_values; percentiles_95],'VariableNames',{'TE_X_YZ','TE_X_ZY','TE_Y_XZ','TE_Y_ZX','TE_Z_XY','TE_Z_YX'},'RowNames', {'TE','p_values','95th_percentile'});
writetable(TE_results,strcat(results_directory_GO,'/national_te.csv'),'WriteRowNames',true);
writetable(distributions,strcat(results_directory_GO,'/national_distributions.csv'));


display('analysis on national level completed');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% state level %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

state_TE = array2table(nan(length(states),7),'VariableNames',{'state','TE_X_YZ','TE_X_ZY','TE_Y_XZ','TE_Y_ZX','TE_Z_XY','TE_Z_YX',});
state_percentiles_95 = array2table(nan(length(states),7),'VariableNames',{'state','TE_X_YZ','TE_X_ZY','TE_Y_XZ','TE_Y_ZX','TE_Z_XY','TE_Z_YX',});
state_p_values = array2table(nan(length(states),7),'VariableNames',{'state','TE_X_YZ','TE_X_ZY','TE_Y_XZ','TE_Y_ZX','TE_Z_XY','TE_Z_YX',});
state_significance = array2table(nan(length(states),7),'VariableNames',{'state','TE_X_YZ','TE_X_ZY','TE_Y_XZ','TE_Y_ZX','TE_Z_XY','TE_Z_YX',});
for s = transpose(states)
    
    % define time series
    X = GO_2000_2017_sa_dt{:,strcmp(GO_2000_2017_sa_dt.Properties.VariableNames,s)};
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
 
    writetable(distributions,strcat(results_directory_GO,'/distributions_',string(s),'.csv'));
    writetable(state_TE,strcat(results_directory_GO,'/state_te.csv'));
    writetable(state_percentiles_95,strcat(results_directory_GO,'/state_percentiles_95.csv'));
    writetable(state_p_values,strcat(results_directory_GO,'/state_p_values.csv'));
    
    display(string(s));

end
state_TE.state = states;
state_percentiles_95.state = states;
state_p_values.state = states;

writetable(state_TE,strcat(results_directory_GO,'/state_te.csv'));
writetable(state_percentiles_95,strcat(results_directory_GO,'/state_percentiles_95.csv'));
writetable(state_p_values,strcat(results_directory_GO,'/state_p_values.csv'));
 

display('analysis on state level completed');



%% plots

plot_state_TE(state_data,states,restrictiveness,results_directory_BC)
plot_state_TE(state_data,states,restrictiveness,results_directory_BCC)
plot_state_TE(state_data,states,restrictiveness,results_directory_S)
plot_state_TE(state_data,states,restrictiveness,results_directory_GO)
