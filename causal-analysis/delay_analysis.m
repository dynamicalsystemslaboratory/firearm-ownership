%% delay analysis

clc; clear;

load('nature.mat');
load('BC_2000_2017_sa_dt.mat');
load('BCC_2000_2017_sa_dt.mat');
load('FO_2000_2017_sa_dt.mat');
load('S_2000_2017_sa_dt.mat');

% initialize variables
states = BC_2000_2017_sa_dt.Properties.VariableNames(3:end-1);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Background checks as a proxy %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% define time series
X = nature.Background_checks(ismember(nature.Year,[2000:2017]));
Y = nature.Mass_shooting(ismember(nature.Year,[2000:2017]));
Z = nature.Firearm_laws_and_regulations(ismember(nature.Year,[2000:2017]));

% symbolize data
XYZ = array2table(double([diff(X)>0 Y(1:end-1)>0 diff(Z)>0]),'VariableNames',{'X','Y','Z'}); % time series for t and t+1

% compute transfer entropy for delays ranging from 0 to 11 months
TE_delays = [];
for d1 = 0:11 % delay on source variable
    for d2 = 0:11 % delay on conditional variable
        TE_delays = [TE_delays; d1 d2 compute_delay_TE(XYZ,d1,d2)];
    end
end
TE_delays = array2table(TE_delays,'VariableNames',{'delay_source','delay_condition','TE_X_YZ','TE_X_ZY','TE_Y_XZ','TE_Y_ZX','TE_Z_XY','TE_Z_YX'});


%%%%%%%%%%%%%%%%%%%% plot the transfer entropy values %%%%%%%%%%%%%%%%%%%%

figure
t=tiledlayout(6,2);

%%% BC to MS, conditioned on MO %%%
nexttile 

title('BC to MS, conditioned on MO');
hold on
for d = 0:11
    m = mean(TE_delays.TE_X_YZ(TE_delays.delay_source==d));
    s = std(TE_delays.TE_X_YZ(TE_delays.delay_source==d));
    errorbar(d,m,s,'Color','black','LineWidth',0.5,'CapSize',0)
    scatter(d,m,'filled','blue');
end
hold off
ylim([0 0.05])
yticks([0:0.01:0.05])
xlim([-0.5 11.5])
xticks([0:11])
xticklabels({'0','1','2','3','4','5','6','7','8','9','10','11'})
xlabel('Delay on source variable','FontSize',10,'Color','black')
ylabel({'Conditional'; 'transfer entropy'},'FontSize',10,'Color','black')


nexttile 

title('BC to MS, conditioned on MO');
hold on
for d = 0:11
    m = mean(TE_delays.TE_X_YZ(TE_delays.delay_condition==d));
    s = std(TE_delays.TE_X_YZ(TE_delays.delay_condition==d));
    errorbar(d,m,s,'Color','black','LineWidth',0.5,'CapSize',0)
    scatter(d,m,'filled','red');
end
hold off
ylim([0 0.05])
yticks([0:0.01:0.05])
xlim([-0.5 11.5])
xticks([0:11])
xticklabels({'0','1','2','3','4','5','6','7','8','9','10','11'})
xlabel('Delay on conditional variable','FontSize',10,'Color','black')
ylabel({'Conditional'; 'transfer entropy'},'FontSize',10,'Color','black')


%%% BC to MO, conditioned on MS %%%
nexttile 

title('BC to MO, conditioned on MS');
hold on
for d = 0:11
    m = mean(TE_delays.TE_X_ZY(TE_delays.delay_source==d));
    s = std(TE_delays.TE_X_ZY(TE_delays.delay_source==d));
    errorbar(d,m,s,'Color','black','LineWidth',0.5,'CapSize',0)
    scatter(d,m,'filled','blue');
end
hold off
ylim([0 0.05])
yticks([0:0.01:0.05])
xlim([-0.5 11.5])
xticks([0:11])
xticklabels({'0','1','2','3','4','5','6','7','8','9','10','11'})
xlabel('Delay on source variable','FontSize',10,'Color','black')
ylabel({'Conditional'; 'transfer entropy'},'FontSize',10,'Color','black')


nexttile 

title('BC to MO, conditioned on MS');
hold on
for d = 0:11
    m = mean(TE_delays.TE_X_ZY(TE_delays.delay_condition==d));
    s = std(TE_delays.TE_X_ZY(TE_delays.delay_condition==d));
    errorbar(d,m,s,'Color','black','LineWidth',0.5,'CapSize',0)
    scatter(d,m,'filled','red');
end
hold off
ylim([0 0.05])
yticks([0:0.01:0.05])
xlim([-0.5 11.5])
xticks([0:11])
xticklabels({'0','1','2','3','4','5','6','7','8','9','10','11'})
xlabel('Delay on conditional variable','FontSize',10,'Color','black')
ylabel({'Conditional'; 'transfer entropy'},'FontSize',10,'Color','black')

%%% MS to BC, conditioned on MO %%%
nexttile 

title('MS to BC, conditioned on MO');
hold on
for d = 0:11
    m = mean(TE_delays.TE_Y_XZ(TE_delays.delay_source==d));
    s = std(TE_delays.TE_Y_XZ(TE_delays.delay_source==d));
    errorbar(d,m,s,'Color','black','LineWidth',0.5,'CapSize',0)
    scatter(d,m,'filled','blue');
end
hold off
ylim([0 0.05])
yticks([0:0.01:0.05])
xlim([-0.5 11.5])
xticks([0:11])
xticklabels({'0','1','2','3','4','5','6','7','8','9','10','11'})
xlabel('Delay on source variable','FontSize',10,'Color','black')
ylabel({'Conditional'; 'transfer entropy'},'FontSize',10,'Color','black')


nexttile 

title('MS to BC, conditioned on MO');
hold on
for d = 0:11
    m = mean(TE_delays.TE_Y_XZ(TE_delays.delay_condition==d));
    s = std(TE_delays.TE_Y_XZ(TE_delays.delay_condition==d));
    errorbar(d,m,s,'Color','black','LineWidth',0.5,'CapSize',0)
    scatter(d,m,'filled','red');
end
hold off
ylim([0 0.05])
yticks([0:0.01:0.05])
xlim([-0.5 11.5])
xticks([0:11])
xticklabels({'0','1','2','3','4','5','6','7','8','9','10','11'})
xlabel('Delay on conditional variable','FontSize',10,'Color','black')
ylabel({'Conditional'; 'transfer entropy'},'FontSize',10,'Color','black')


%%% MS to MO, conditioned on BC %%%
nexttile 

title('MS to MO, conditioned on BC');
hold on
for d = 0:11
    m = mean(TE_delays.TE_Y_ZX(TE_delays.delay_source==d));
    s = std(TE_delays.TE_Y_ZX(TE_delays.delay_source==d));
    errorbar(d,m,s,'Color','black','LineWidth',0.5,'CapSize',0)
    scatter(d,m,'filled','blue');
end
hold off
ylim([0 0.05])
yticks([0:0.01:0.05])
xlim([-0.5 11.5])
xticks([0:11])
xticklabels({'0','1','2','3','4','5','6','7','8','9','10','11'})
xlabel('Delay on source variable','FontSize',10,'Color','black')
ylabel({'Conditional'; 'transfer entropy'},'FontSize',10,'Color','black')


nexttile 

title('MS to MO, conditioned on BC');
hold on
for d = 0:11
    m = mean(TE_delays.TE_Y_ZX(TE_delays.delay_condition==d));
    s = std(TE_delays.TE_Y_ZX(TE_delays.delay_condition==d));
    errorbar(d,m,s,'Color','black','LineWidth',0.5,'CapSize',0)
    scatter(d,m,'filled','red');
end
hold off
ylim([0 0.05])
yticks([0:0.01:0.05])
xlim([-0.5 11.5])
xticks([0:11])
xticklabels({'0','1','2','3','4','5','6','7','8','9','10','11'})
xlabel('Delay on conditional variable','FontSize',10,'Color','black')
ylabel({'Conditional'; 'transfer entropy'},'FontSize',10,'Color','black')



%%% MO to BC, conditioned on MS %%%
nexttile 

title('MO to BC, conditioned on MS');
hold on
for d = 0:11
    m = mean(TE_delays.TE_Z_XY(TE_delays.delay_source==d));
    s = std(TE_delays.TE_Z_XY(TE_delays.delay_source==d));
    errorbar(d,m,s,'Color','black','LineWidth',0.5,'CapSize',0)
    scatter(d,m,'filled','blue');
end
hold off
ylim([0 0.05])
yticks([0:0.01:0.05])
xlim([-0.5 11.5])
xticks([0:11])
xticklabels({'0','1','2','3','4','5','6','7','8','9','10','11'})
xlabel('Delay on source variable','FontSize',10,'Color','black')
ylabel({'Conditional'; 'transfer entropy'},'FontSize',10,'Color','black')


nexttile 

title('MO to BC, conditioned on MS');
hold on
for d = 0:11
    m = mean(TE_delays.TE_Z_XY(TE_delays.delay_condition==d));
    s = std(TE_delays.TE_Z_XY(TE_delays.delay_condition==d));
    errorbar(d,m,s,'Color','black','LineWidth',0.5,'CapSize',0)
    scatter(d,m,'filled','red');
end
hold off
ylim([0 0.05])
yticks([0:0.01:0.05])
xlim([-0.5 11.5])
xticks([0:11])
xticklabels({'0','1','2','3','4','5','6','7','8','9','10','11'})
xlabel('Delay on conditional variable','FontSize',10,'Color','black')
ylabel({'Conditional'; 'transfer entropy'},'FontSize',10,'Color','black')


%%% MO to MS, conditioned on BC %%%
nexttile 

title('MO to MS, conditioned on BC');
hold on
for d = 0:11
    m = mean(TE_delays.TE_Z_YX(TE_delays.delay_source==d));
    s = std(TE_delays.TE_Z_YX(TE_delays.delay_source==d));
    errorbar(d,m,s,'Color','black','LineWidth',0.5,'CapSize',0)
    scatter(d,m,'filled','blue');
end
hold off
ylim([0 0.05])
yticks([0:0.01:0.05])
xlim([-0.5 11.5])
xticks([0:11])
xticklabels({'0','1','2','3','4','5','6','7','8','9','10','11'})
xlabel('Delay on source variable','FontSize',10,'Color','black')
ylabel({'Conditional'; 'transfer entropy'},'FontSize',10,'Color','black')


nexttile 

title('MO to MS, conditioned on BC');
hold on
for d = 0:11
    m = mean(TE_delays.TE_Z_YX(TE_delays.delay_condition==d));
    s = std(TE_delays.TE_Z_YX(TE_delays.delay_condition==d));
    errorbar(d,m,s,'Color','black','LineWidth',0.5,'CapSize',0)
    scatter(d,m,'filled','red');
end
hold off
ylim([0 0.05])
yticks([0:0.01:0.05])
xlim([-0.5 11.5])
xticks([0:11])
xticklabels({'0','1','2','3','4','5','6','7','8','9','10','11'})
xlabel('Delay on conditional variable','FontSize',10,'Color','black')
ylabel({'Conditional'; 'transfer entropy'},'FontSize',10,'Color','black')


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Background checks per capita as a proxy %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% define time series
X = BCC_2000_2017_sa_dt.USA;
Y = nature.Mass_shooting(ismember(nature.Year,[2000:2017]));
Z = nature.Firearm_laws_and_regulations(ismember(nature.Year,[2000:2017]));

% symbolize data
XYZ = array2table(double([diff(X)>0 Y(1:end-1)>0 diff(Z)>0]),'VariableNames',{'X','Y','Z'}); % time series for t and t+1

% compute transfer entropy for delays ranging from 0 to 11 months
TE_delays = [];
for d1 = 0:11 % delay on source variable
    for d2 = 0:11 % delay on conditional variable
        TE_delays = [TE_delays; d1 d2 compute_delay_TE(XYZ,d1,d2)];
    end
end
TE_delays = array2table(TE_delays,'VariableNames',{'delay_source','delay_condition','TE_X_YZ','TE_X_ZY','TE_Y_XZ','TE_Y_ZX','TE_Z_XY','TE_Z_YX'});


%%%%%%%%%%%%%%%%%%%% plot the transfer entropy values %%%%%%%%%%%%%%%%%%%%

figure
t=tiledlayout(6,2);

%%% BCC to MS, conditioned on MO %%%
nexttile 

title('BCC to MS, conditioned on MO');
hold on
for d = 0:11
    m = mean(TE_delays.TE_X_YZ(TE_delays.delay_source==d));
    s = std(TE_delays.TE_X_YZ(TE_delays.delay_source==d));
    errorbar(d,m,s,'Color','black','LineWidth',0.5,'CapSize',0)
    scatter(d,m,'filled','blue');
end
hold off
ylim([0 0.05])
yticks([0:0.01:0.05])
xlim([-0.5 11.5])
xticks([0:11])
xticklabels({'0','1','2','3','4','5','6','7','8','9','10','11'})
xlabel('Delay on source variable','FontSize',10,'Color','black')
ylabel({'Conditional'; 'transfer entropy'},'FontSize',10,'Color','black')


nexttile 

title('BCC to MS, conditioned on MO');
hold on
for d = 0:11
    m = mean(TE_delays.TE_X_YZ(TE_delays.delay_condition==d));
    s = std(TE_delays.TE_X_YZ(TE_delays.delay_condition==d));
    errorbar(d,m,s,'Color','black','LineWidth',0.5,'CapSize',0)
    scatter(d,m,'filled','red');
end
hold off
ylim([0 0.05])
yticks([0:0.01:0.05])
xlim([-0.5 11.5])
xticks([0:11])
xticklabels({'0','1','2','3','4','5','6','7','8','9','10','11'})
xlabel('Delay on conditional variable','FontSize',10,'Color','black')
ylabel({'Conditional'; 'transfer entropy'},'FontSize',10,'Color','black')


%%% BCC to MO, conditioned on MS %%%
nexttile 

title('BCC to MO, conditioned on MS');
hold on
for d = 0:11
    m = mean(TE_delays.TE_X_ZY(TE_delays.delay_source==d));
    s = std(TE_delays.TE_X_ZY(TE_delays.delay_source==d));
    errorbar(d,m,s,'Color','black','LineWidth',0.5,'CapSize',0)
    scatter(d,m,'filled','blue');
end
hold off
ylim([0 0.05])
yticks([0:0.01:0.05])
xlim([-0.5 11.5])
xticks([0:11])
xticklabels({'0','1','2','3','4','5','6','7','8','9','10','11'})
xlabel('Delay on source variable','FontSize',10,'Color','black')
ylabel({'Conditional'; 'transfer entropy'},'FontSize',10,'Color','black')


nexttile 

title('BCC to MO, conditioned on MS');
hold on
for d = 0:11
    m = mean(TE_delays.TE_X_ZY(TE_delays.delay_condition==d));
    s = std(TE_delays.TE_X_ZY(TE_delays.delay_condition==d));
    errorbar(d,m,s,'Color','black','LineWidth',0.5,'CapSize',0)
    scatter(d,m,'filled','red');
end
hold off
ylim([0 0.05])
yticks([0:0.01:0.05])
xlim([-0.5 11.5])
xticks([0:11])
xticklabels({'0','1','2','3','4','5','6','7','8','9','10','11'})
xlabel('Delay on conditional variable','FontSize',10,'Color','black')
ylabel({'Conditional'; 'transfer entropy'},'FontSize',10,'Color','black')

%%% MS to BCC, conditioned on MO %%%
nexttile 

title('MS to BCC, conditioned on MO');
hold on
for d = 0:11
    m = mean(TE_delays.TE_Y_XZ(TE_delays.delay_source==d));
    s = std(TE_delays.TE_Y_XZ(TE_delays.delay_source==d));
    errorbar(d,m,s,'Color','black','LineWidth',0.5,'CapSize',0)
    scatter(d,m,'filled','blue');
end
hold off
ylim([0 0.05])
yticks([0:0.01:0.05])
xlim([-0.5 11.5])
xticks([0:11])
xticklabels({'0','1','2','3','4','5','6','7','8','9','10','11'})
xlabel('Delay on source variable','FontSize',10,'Color','black')
ylabel({'Conditional'; 'transfer entropy'},'FontSize',10,'Color','black')


nexttile 

title('MS to BCC, conditioned on MO');
hold on
for d = 0:11
    m = mean(TE_delays.TE_Y_XZ(TE_delays.delay_condition==d));
    s = std(TE_delays.TE_Y_XZ(TE_delays.delay_condition==d));
    errorbar(d,m,s,'Color','black','LineWidth',0.5,'CapSize',0)
    scatter(d,m,'filled','red');
end
hold off
ylim([0 0.05])
yticks([0:0.01:0.05])
xlim([-0.5 11.5])
xticks([0:11])
xticklabels({'0','1','2','3','4','5','6','7','8','9','10','11'})
xlabel('Delay on conditional variable','FontSize',10,'Color','black')
ylabel({'Conditional'; 'transfer entropy'},'FontSize',10,'Color','black')


%%% MS to MO, conditioned on BCC %%%
nexttile 

title('MS to MO, conditioned on BCC');
hold on
for d = 0:11
    m = mean(TE_delays.TE_Y_ZX(TE_delays.delay_source==d));
    s = std(TE_delays.TE_Y_ZX(TE_delays.delay_source==d));
    errorbar(d,m,s,'Color','black','LineWidth',0.5,'CapSize',0)
    scatter(d,m,'filled','blue');
end
hold off
ylim([0 0.05])
yticks([0:0.01:0.05])
xlim([-0.5 11.5])
xticks([0:11])
xticklabels({'0','1','2','3','4','5','6','7','8','9','10','11'})
xlabel('Delay on source variable','FontSize',10,'Color','black')
ylabel({'Conditional'; 'transfer entropy'},'FontSize',10,'Color','black')


nexttile 

title('MS to MO, conditioned on BCC');
hold on
for d = 0:11
    m = mean(TE_delays.TE_Y_ZX(TE_delays.delay_condition==d));
    s = std(TE_delays.TE_Y_ZX(TE_delays.delay_condition==d));
    errorbar(d,m,s,'Color','black','LineWidth',0.5,'CapSize',0)
    scatter(d,m,'filled','red');
end
hold off
ylim([0 0.05])
yticks([0:0.01:0.05])
xlim([-0.5 11.5])
xticks([0:11])
xticklabels({'0','1','2','3','4','5','6','7','8','9','10','11'})
xlabel('Delay on conditional variable','FontSize',10,'Color','black')
ylabel({'Conditional'; 'transfer entropy'},'FontSize',10,'Color','black')



%%% MO to BCC, conditioned on MS %%%
nexttile 

title('MO to BCC, conditioned on MS');
hold on
for d = 0:11
    m = mean(TE_delays.TE_Z_XY(TE_delays.delay_source==d));
    s = std(TE_delays.TE_Z_XY(TE_delays.delay_source==d));
    errorbar(d,m,s,'Color','black','LineWidth',0.5,'CapSize',0)
    scatter(d,m,'filled','blue');
end
hold off
ylim([0 0.05])
yticks([0:0.01:0.05])
xlim([-0.5 11.5])
xticks([0:11])
xticklabels({'0','1','2','3','4','5','6','7','8','9','10','11'})
xlabel('Delay on source variable','FontSize',10,'Color','black')
ylabel({'Conditional'; 'transfer entropy'},'FontSize',10,'Color','black')


nexttile 

title('MO to BCC, conditioned on MS');
hold on
for d = 0:11
    m = mean(TE_delays.TE_Z_XY(TE_delays.delay_condition==d));
    s = std(TE_delays.TE_Z_XY(TE_delays.delay_condition==d));
    errorbar(d,m,s,'Color','black','LineWidth',0.5,'CapSize',0)
    scatter(d,m,'filled','red');
end
hold off
ylim([0 0.05])
yticks([0:0.01:0.05])
xlim([-0.5 11.5])
xticks([0:11])
xticklabels({'0','1','2','3','4','5','6','7','8','9','10','11'})
xlabel('Delay on conditional variable','FontSize',10,'Color','black')
ylabel({'Conditional'; 'transfer entropy'},'FontSize',10,'Color','black')


%%% MO to MS, conditioned on BCC %%%
nexttile 

title('MO to MS, conditioned on BCC');
hold on
for d = 0:11
    m = mean(TE_delays.TE_Z_YX(TE_delays.delay_source==d));
    s = std(TE_delays.TE_Z_YX(TE_delays.delay_source==d));
    errorbar(d,m,s,'Color','black','LineWidth',0.5,'CapSize',0)
    scatter(d,m,'filled','blue');
end
hold off
ylim([0 0.05])
yticks([0:0.01:0.05])
xlim([-0.5 11.5])
xticks([0:11])
xticklabels({'0','1','2','3','4','5','6','7','8','9','10','11'})
xlabel('Delay on source variable','FontSize',10,'Color','black')
ylabel({'Conditional'; 'transfer entropy'},'FontSize',10,'Color','black')


nexttile 

title('MO to MS, conditioned on BCC');
hold on
for d = 0:11
    m = mean(TE_delays.TE_Z_YX(TE_delays.delay_condition==d));
    s = std(TE_delays.TE_Z_YX(TE_delays.delay_condition==d));
    errorbar(d,m,s,'Color','black','LineWidth',0.5,'CapSize',0)
    scatter(d,m,'filled','red');
end
hold off
ylim([0 0.05])
yticks([0:0.01:0.05])
xlim([-0.5 11.5])
xticks([0:11])
xticklabels({'0','1','2','3','4','5','6','7','8','9','10','11'})
xlabel('Delay on conditional variable','FontSize',10,'Color','black')
ylabel({'Conditional'; 'transfer entropy'},'FontSize',10,'Color','black')


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Suicides with firearms as a proxy %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% define time series
X = S_2000_2017_sa_dt.USA;
Y = nature.Mass_shooting(ismember(nature.Year,[2000:2017]));
Z = nature.Firearm_laws_and_regulations(ismember(nature.Year,[2000:2017]));

% symbolize data
XYZ = array2table(double([diff(X)>0 Y(1:end-1)>0 diff(Z)>0]),'VariableNames',{'X','Y','Z'}); % time series for t and t+1

% compute transfer entropy for delays ranging from 0 to 11 months
TE_delays = [];
for d1 = 0:11 % delay on source variable
    for d2 = 0:11 % delay on conditional variable
        TE_delays = [TE_delays; d1 d2 compute_delay_TE(XYZ,d1,d2)];
    end
end
TE_delays = array2table(TE_delays,'VariableNames',{'delay_source','delay_condition','TE_X_YZ','TE_X_ZY','TE_Y_XZ','TE_Y_ZX','TE_Z_XY','TE_Z_YX'});


%%%%%%%%%%%%%%%%%%%% plot the transfer entropy values %%%%%%%%%%%%%%%%%%%%

figure
t=tiledlayout(6,2);

%%% S to MS, conditioned on MO %%%
nexttile 

title('S to MS, conditioned on MO');
hold on
for d = 0:11
    m = mean(TE_delays.TE_X_YZ(TE_delays.delay_source==d));
    s = std(TE_delays.TE_X_YZ(TE_delays.delay_source==d));
    errorbar(d,m,s,'Color','black','LineWidth',0.5,'CapSize',0)
    scatter(d,m,'filled','blue');
end
hold off
ylim([0 0.05])
yticks([0:0.01:0.05])
xlim([-0.5 11.5])
xticks([0:11])
xticklabels({'0','1','2','3','4','5','6','7','8','9','10','11'})
xlabel('Delay on source variable','FontSize',10,'Color','black')
ylabel({'Conditional'; 'transfer entropy'},'FontSize',10,'Color','black')


nexttile 

title('S to MS, conditioned on MO');
hold on
for d = 0:11
    m = mean(TE_delays.TE_X_YZ(TE_delays.delay_condition==d));
    s = std(TE_delays.TE_X_YZ(TE_delays.delay_condition==d));
    errorbar(d,m,s,'Color','black','LineWidth',0.5,'CapSize',0)
    scatter(d,m,'filled','red');
end
hold off
ylim([0 0.05])
yticks([0:0.01:0.05])
xlim([-0.5 11.5])
xticks([0:11])
xticklabels({'0','1','2','3','4','5','6','7','8','9','10','11'})
xlabel('Delay on conditional variable','FontSize',10,'Color','black')
ylabel({'Conditional'; 'transfer entropy'},'FontSize',10,'Color','black')


%%% S to MO, conditioned on MS %%%
nexttile 

title('S to MO, conditioned on MS');
hold on
for d = 0:11
    m = mean(TE_delays.TE_X_ZY(TE_delays.delay_source==d));
    s = std(TE_delays.TE_X_ZY(TE_delays.delay_source==d));
    errorbar(d,m,s,'Color','black','LineWidth',0.5,'CapSize',0)
    scatter(d,m,'filled','blue');
end
hold off
ylim([0 0.05])
yticks([0:0.01:0.05])
xlim([-0.5 11.5])
xticks([0:11])
xticklabels({'0','1','2','3','4','5','6','7','8','9','10','11'})
xlabel('Delay on source variable','FontSize',10,'Color','black')
ylabel({'Conditional'; 'transfer entropy'},'FontSize',10,'Color','black')


nexttile 

title('S to MO, conditioned on MS');
hold on
for d = 0:11
    m = mean(TE_delays.TE_X_ZY(TE_delays.delay_condition==d));
    s = std(TE_delays.TE_X_ZY(TE_delays.delay_condition==d));
    errorbar(d,m,s,'Color','black','LineWidth',0.5,'CapSize',0)
    scatter(d,m,'filled','red');
end
hold off
ylim([0 0.05])
yticks([0:0.01:0.05])
xlim([-0.5 11.5])
xticks([0:11])
xticklabels({'0','1','2','3','4','5','6','7','8','9','10','11'})
xlabel('Delay on conditional variable','FontSize',10,'Color','black')
ylabel({'Conditional'; 'transfer entropy'},'FontSize',10,'Color','black')

%%% MS to S, conditioned on MO %%%
nexttile 

title('MS to S, conditioned on MO');
hold on
for d = 0:11
    m = mean(TE_delays.TE_Y_XZ(TE_delays.delay_source==d));
    s = std(TE_delays.TE_Y_XZ(TE_delays.delay_source==d));
    errorbar(d,m,s,'Color','black','LineWidth',0.5,'CapSize',0)
    scatter(d,m,'filled','blue');
end
hold off
ylim([0 0.05])
yticks([0:0.01:0.05])
xlim([-0.5 11.5])
xticks([0:11])
xticklabels({'0','1','2','3','4','5','6','7','8','9','10','11'})
xlabel('Delay on source variable','FontSize',10,'Color','black')
ylabel({'Conditional'; 'transfer entropy'},'FontSize',10,'Color','black')


nexttile 

title('MS to S, conditioned on MO');
hold on
for d = 0:11
    m = mean(TE_delays.TE_Y_XZ(TE_delays.delay_condition==d));
    s = std(TE_delays.TE_Y_XZ(TE_delays.delay_condition==d));
    errorbar(d,m,s,'Color','black','LineWidth',0.5,'CapSize',0)
    scatter(d,m,'filled','red');
end
hold off
ylim([0 0.05])
yticks([0:0.01:0.05])
xlim([-0.5 11.5])
xticks([0:11])
xticklabels({'0','1','2','3','4','5','6','7','8','9','10','11'})
xlabel('Delay on conditional variable','FontSize',10,'Color','black')
ylabel({'Conditional'; 'transfer entropy'},'FontSize',10,'Color','black')


%%% MS to MO, conditioned on S %%%
nexttile 

title('MS to MO, conditioned on S');
hold on
for d = 0:11
    m = mean(TE_delays.TE_Y_ZX(TE_delays.delay_source==d));
    s = std(TE_delays.TE_Y_ZX(TE_delays.delay_source==d));
    errorbar(d,m,s,'Color','black','LineWidth',0.5,'CapSize',0)
    scatter(d,m,'filled','blue');
end
hold off
ylim([0 0.05])
yticks([0:0.01:0.05])
xlim([-0.5 11.5])
xticks([0:11])
xticklabels({'0','1','2','3','4','5','6','7','8','9','10','11'})
xlabel('Delay on source variable','FontSize',10,'Color','black')
ylabel({'Conditional'; 'transfer entropy'},'FontSize',10,'Color','black')


nexttile 

title('MS to MO, conditioned on S');
hold on
for d = 0:11
    m = mean(TE_delays.TE_Y_ZX(TE_delays.delay_condition==d));
    s = std(TE_delays.TE_Y_ZX(TE_delays.delay_condition==d));
    errorbar(d,m,s,'Color','black','LineWidth',0.5,'CapSize',0)
    scatter(d,m,'filled','red');
end
hold off
ylim([0 0.05])
yticks([0:0.01:0.05])
xlim([-0.5 11.5])
xticks([0:11])
xticklabels({'0','1','2','3','4','5','6','7','8','9','10','11'})
xlabel('Delay on conditional variable','FontSize',10,'Color','black')
ylabel({'Conditional'; 'transfer entropy'},'FontSize',10,'Color','black')



%%% MO to S, conditioned on MS %%%
nexttile 

title('MO to S, conditioned on MS');
hold on
for d = 0:11
    m = mean(TE_delays.TE_Z_XY(TE_delays.delay_source==d));
    s = std(TE_delays.TE_Z_XY(TE_delays.delay_source==d));
    errorbar(d,m,s,'Color','black','LineWidth',0.5,'CapSize',0)
    scatter(d,m,'filled','blue');
end
hold off
ylim([0 0.05])
yticks([0:0.01:0.05])
xlim([-0.5 11.5])
xticks([0:11])
xticklabels({'0','1','2','3','4','5','6','7','8','9','10','11'})
xlabel('Delay on source variable','FontSize',10,'Color','black')
ylabel({'Conditional'; 'transfer entropy'},'FontSize',10,'Color','black')


nexttile 

title('MO to S, conditioned on MS');
hold on
for d = 0:11
    m = mean(TE_delays.TE_Z_XY(TE_delays.delay_condition==d));
    s = std(TE_delays.TE_Z_XY(TE_delays.delay_condition==d));
    errorbar(d,m,s,'Color','black','LineWidth',0.5,'CapSize',0)
    scatter(d,m,'filled','red');
end
hold off
ylim([0 0.05])
yticks([0:0.01:0.05])
xlim([-0.5 11.5])
xticks([0:11])
xticklabels({'0','1','2','3','4','5','6','7','8','9','10','11'})
xlabel('Delay on conditional variable','FontSize',10,'Color','black')
ylabel({'Conditional'; 'transfer entropy'},'FontSize',10,'Color','black')


%%% MO to MS, conditioned on S %%%
nexttile 

title('MO to MS, conditioned on S');
hold on
for d = 0:11
    m = mean(TE_delays.TE_Z_YX(TE_delays.delay_source==d));
    s = std(TE_delays.TE_Z_YX(TE_delays.delay_source==d));
    errorbar(d,m,s,'Color','black','LineWidth',0.5,'CapSize',0)
    scatter(d,m,'filled','blue');
end
hold off
ylim([0 0.05])
yticks([0:0.01:0.05])
xlim([-0.5 11.5])
xticks([0:11])
xticklabels({'0','1','2','3','4','5','6','7','8','9','10','11'})
xlabel('Delay on source variable','FontSize',10,'Color','black')
ylabel({'Conditional'; 'transfer entropy'},'FontSize',10,'Color','black')


nexttile 

title('MO to MS, conditioned on S');
hold on
for d = 0:11
    m = mean(TE_delays.TE_Z_YX(TE_delays.delay_condition==d));
    s = std(TE_delays.TE_Z_YX(TE_delays.delay_condition==d));
    errorbar(d,m,s,'Color','black','LineWidth',0.5,'CapSize',0)
    scatter(d,m,'filled','red');
end
hold off
ylim([0 0.05])
yticks([0:0.01:0.05])
xlim([-0.5 11.5])
xticks([0:11])
xticklabels({'0','1','2','3','4','5','6','7','8','9','10','11'})
xlabel('Delay on conditional variable','FontSize',10,'Color','black')
ylabel({'Conditional'; 'transfer entropy'},'FontSize',10,'Color','black')


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% Firearm ownership as a proxy %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% define time series
X = FO_2000_2017_sa_dt.USA;
Y = nature.Mass_shooting(ismember(nature.Year,[2000:2017]));
Z = nature.Firearm_laws_and_regulations(ismember(nature.Year,[2000:2017]));

% symbolize data
XYZ = array2table(double([diff(X)>0 Y(1:end-1)>0 diff(Z)>0]),'VariableNames',{'X','Y','Z'}); % time series for t and t+1

% compute transfer entropy for delays ranging from 0 to 11 months
TE_delays = [];
for d1 = 0:11 % delay on source variable
    for d2 = 0:11 % delay on conditional variable
        TE_delays = [TE_delays; d1 d2 compute_delay_TE(XYZ,d1,d2)];
    end
end
TE_delays = array2table(TE_delays,'VariableNames',{'delay_source','delay_condition','TE_X_YZ','TE_X_ZY','TE_Y_XZ','TE_Y_ZX','TE_Z_XY','TE_Z_YX'});


%%%%%%%%%%%%%%%%%%%% plot the transfer entropy values %%%%%%%%%%%%%%%%%%%%

figure
t=tiledlayout(6,2);

%%% FO to MS, conditioned on MO %%%
nexttile 

title('FO to MS, conditioned on MO');
hold on
for d = 0:11
    m = mean(TE_delays.TE_X_YZ(TE_delays.delay_source==d));
    s = std(TE_delays.TE_X_YZ(TE_delays.delay_source==d));
    errorbar(d,m,s,'Color','black','LineWidth',0.5,'CapSize',0)
    scatter(d,m,'filled','blue');
end
hold off
ylim([0 0.05])
yticks([0:0.01:0.05])
xlim([-0.5 11.5])
xticks([0:11])
xticklabels({'0','1','2','3','4','5','6','7','8','9','10','11'})
xlabel('Delay on source variable','FontSize',10,'Color','black')
ylabel({'Conditional'; 'transfer entropy'},'FontSize',10,'Color','black')


nexttile 

title('FO to MS, conditioned on MO');
hold on
for d = 0:11
    m = mean(TE_delays.TE_X_YZ(TE_delays.delay_condition==d));
    s = std(TE_delays.TE_X_YZ(TE_delays.delay_condition==d));
    errorbar(d,m,s,'Color','black','LineWidth',0.5,'CapSize',0)
    scatter(d,m,'filled','red');
end
hold off
ylim([0 0.05])
yticks([0:0.01:0.05])
xlim([-0.5 11.5])
xticks([0:11])
xticklabels({'0','1','2','3','4','5','6','7','8','9','10','11'})
xlabel('Delay on conditional variable','FontSize',10,'Color','black')
ylabel({'Conditional'; 'transfer entropy'},'FontSize',10,'Color','black')


%%% FO to MO, conditioned on MS %%%
nexttile 

title('FO to MO, conditioned on MS');
hold on
for d = 0:11
    m = mean(TE_delays.TE_X_ZY(TE_delays.delay_source==d));
    s = std(TE_delays.TE_X_ZY(TE_delays.delay_source==d));
    errorbar(d,m,s,'Color','black','LineWidth',0.5,'CapSize',0)
    scatter(d,m,'filled','blue');
end
hold off
ylim([0 0.1])
yticks([0:0.02:0.1])
xlim([-0.5 11.5])
xticks([0:11])
xticklabels({'0','1','2','3','4','5','6','7','8','9','10','11'})
xlabel('Delay on source variable','FontSize',10,'Color','black')
ylabel({'Conditional'; 'transfer entropy'},'FontSize',10,'Color','black')


nexttile 

title('FO to MO, conditioned on MS');
hold on
for d = 0:11
    m = mean(TE_delays.TE_X_ZY(TE_delays.delay_condition==d));
    s = std(TE_delays.TE_X_ZY(TE_delays.delay_condition==d));
    errorbar(d,m,s,'Color','black','LineWidth',0.5,'CapSize',0)
    scatter(d,m,'filled','red');
end
hold off
ylim([0 0.05])
yticks([0:0.01:0.05])
xlim([-0.5 11.5])
xticks([0:11])
xticklabels({'0','1','2','3','4','5','6','7','8','9','10','11'})
xlabel('Delay on conditional variable','FontSize',10,'Color','black')
ylabel({'Conditional'; 'transfer entropy'},'FontSize',10,'Color','black')

%%% MS to FO, conditioned on MO %%%
nexttile 

title('MS to FO, conditioned on MO');
hold on
for d = 0:11
    m = mean(TE_delays.TE_Y_XZ(TE_delays.delay_source==d));
    s = std(TE_delays.TE_Y_XZ(TE_delays.delay_source==d));
    errorbar(d,m,s,'Color','black','LineWidth',0.5,'CapSize',0)
    scatter(d,m,'filled','blue');
end
hold off
ylim([0 0.05])
yticks([0:0.01:0.05])
xlim([-0.5 11.5])
xticks([0:11])
xticklabels({'0','1','2','3','4','5','6','7','8','9','10','11'})
xlabel('Delay on source variable','FontSize',10,'Color','black')
ylabel({'Conditional'; 'transfer entropy'},'FontSize',10,'Color','black')


nexttile 

title('MS to FO, conditioned on MO');
hold on
for d = 0:11
    m = mean(TE_delays.TE_Y_XZ(TE_delays.delay_condition==d));
    s = std(TE_delays.TE_Y_XZ(TE_delays.delay_condition==d));
    errorbar(d,m,s,'Color','black','LineWidth',0.5,'CapSize',0)
    scatter(d,m,'filled','red');
end
hold off
ylim([0 0.05])
yticks([0:0.01:0.05])
xlim([-0.5 11.5])
xticks([0:11])
xticklabels({'0','1','2','3','4','5','6','7','8','9','10','11'})
xlabel('Delay on conditional variable','FontSize',10,'Color','black')
ylabel({'Conditional'; 'transfer entropy'},'FontSize',10,'Color','black')


%%% MS to MO, conditioned on FO %%%
nexttile 

title('MS to MO, conditioned on FO');
hold on
for d = 0:11
    m = mean(TE_delays.TE_Y_ZX(TE_delays.delay_source==d));
    s = std(TE_delays.TE_Y_ZX(TE_delays.delay_source==d));
    errorbar(d,m,s,'Color','black','LineWidth',0.5,'CapSize',0)
    scatter(d,m,'filled','blue');
end
hold off
ylim([0 0.05])
yticks([0:0.01:0.05])
xlim([-0.5 11.5])
xticks([0:11])
xticklabels({'0','1','2','3','4','5','6','7','8','9','10','11'})
xlabel('Delay on source variable','FontSize',10,'Color','black')
ylabel({'Conditional'; 'transfer entropy'},'FontSize',10,'Color','black')


nexttile 

title('MS to MO, conditioned on FO');
hold on
for d = 0:11
    m = mean(TE_delays.TE_Y_ZX(TE_delays.delay_condition==d));
    s = std(TE_delays.TE_Y_ZX(TE_delays.delay_condition==d));
    errorbar(d,m,s,'Color','black','LineWidth',0.5,'CapSize',0)
    scatter(d,m,'filled','red');
end
hold off
ylim([0 0.05])
yticks([0:0.01:0.05])
xlim([-0.5 11.5])
xticks([0:11])
xticklabels({'0','1','2','3','4','5','6','7','8','9','10','11'})
xlabel('Delay on conditional variable','FontSize',10,'Color','black')
ylabel({'Conditional'; 'transfer entropy'},'FontSize',10,'Color','black')



%%% MO to FO, conditioned on MS %%%
nexttile 

title('MO to FO, conditioned on MS');
hold on
for d = 0:11
    m = mean(TE_delays.TE_Z_XY(TE_delays.delay_source==d));
    s = std(TE_delays.TE_Z_XY(TE_delays.delay_source==d));
    errorbar(d,m,s,'Color','black','LineWidth',0.5,'CapSize',0)
    scatter(d,m,'filled','blue');
end
hold off
ylim([0 0.05])
yticks([0:0.01:0.05])
xlim([-0.5 11.5])
xticks([0:11])
xticklabels({'0','1','2','3','4','5','6','7','8','9','10','11'})
xlabel('Delay on source variable','FontSize',10,'Color','black')
ylabel({'Conditional'; 'transfer entropy'},'FontSize',10,'Color','black')


nexttile 

title('MO to FO, conditioned on MS');
hold on
for d = 0:11
    m = mean(TE_delays.TE_Z_XY(TE_delays.delay_condition==d));
    s = std(TE_delays.TE_Z_XY(TE_delays.delay_condition==d));
    errorbar(d,m,s,'Color','black','LineWidth',0.5,'CapSize',0)
    scatter(d,m,'filled','red');
end
hold off
ylim([0 0.05])
yticks([0:0.01:0.05])
xlim([-0.5 11.5])
xticks([0:11])
xticklabels({'0','1','2','3','4','5','6','7','8','9','10','11'})
xlabel('Delay on conditional variable','FontSize',10,'Color','black')
ylabel({'Conditional'; 'transfer entropy'},'FontSize',10,'Color','black')


%%% MO to MS, conditioned on FO %%%
nexttile 

title('MO to MS, conditioned on FO');
hold on
for d = 0:11
    m = mean(TE_delays.TE_Z_YX(TE_delays.delay_source==d));
    s = std(TE_delays.TE_Z_YX(TE_delays.delay_source==d));
    errorbar(d,m,s,'Color','black','LineWidth',0.5,'CapSize',0)
    scatter(d,m,'filled','blue');
end
hold off
ylim([0 0.05])
yticks([0:0.01:0.05])
xlim([-0.5 11.5])
xticks([0:11])
xticklabels({'0','1','2','3','4','5','6','7','8','9','10','11'})
xlabel('Delay on source variable','FontSize',10,'Color','black')
ylabel({'Conditional'; 'transfer entropy'},'FontSize',10,'Color','black')


nexttile 

title('MO to MS, conditioned on FO');
hold on
for d = 0:11
    m = mean(TE_delays.TE_Z_YX(TE_delays.delay_condition==d));
    s = std(TE_delays.TE_Z_YX(TE_delays.delay_condition==d));
    errorbar(d,m,s,'Color','black','LineWidth',0.5,'CapSize',0)
    scatter(d,m,'filled','red');
end
hold off
ylim([0 0.05])
yticks([0:0.01:0.05])
xlim([-0.5 11.5])
xticks([0:11])
xticklabels({'0','1','2','3','4','5','6','7','8','9','10','11'})
xlabel('Delay on conditional variable','FontSize',10,'Color','black')
ylabel({'Conditional'; 'transfer entropy'},'FontSize',10,'Color','black')


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Firearm ownership as a proxy, state level %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('data.mat');

% compute delays for each state
TE_delays = array2table(nan(144*48,10),'VariableNames',{'state','population','delay_source','delay_condition','TE_X_YZ','TE_X_ZY','TE_Y_XZ','TE_Y_ZX','TE_Z_XY','TE_Z_YX'});
TE_delays.state = transpose(repelem(states,144));
for s = states
    
    % define time series
    X = FO_2000_2017_sa_dt{:,strcmp(FO_2000_2017_sa_dt.Properties.VariableNames,s)};
    Y = nature.Mass_shooting(ismember(nature.Year,[2000:2017]));
    Z = nature.Firearm_laws_and_regulations(ismember(nature.Year,[2000:2017]));

    % symbolize data
    XYZ = array2table(double([diff(X)>0 Y(1:end-1)>0 diff(Z)>0]),'VariableNames',{'X','Y','Z'}); % time series for t and t+1

    % compute transfer entropy for delays ranging from 0 to 11 months
    for d1 = 0:11 % delay on source variable
        for d2 = 0:11 % delay on conditional variable
            TE_delays{(find(strcmp(states,s))-1)*144+d1*12+d2+1,2:end} = [mean(data.population(strcmp(data.state,s))) d1 d2 compute_delay_TE(XYZ,d1,d2)];
        end
    end    
end


% plot delays for permissive and restrictive states
permissive_states = states([1:3 5 7:10 12:17 20:27 29 31:48]);
restrictive_states = states([4 6 11 18 19 28 30]);

permissive_TE_delays = TE_delays(ismember(TE_delays.state,permissive_states),:);
restrictive_TE_delays = TE_delays(ismember(TE_delays.state,restrictive_states),:);

figure
hold on
for d = 0:11
    
    m = (mean(reshape(permissive_TE_delays.TE_X_ZY(logical((ismember(permissive_TE_delays.state,permissive_states)).*(permissive_TE_delays.delay_source==d))),12,size(permissive_states,2)),1))*(permissive_TE_delays.population(1:144:end)/sum(permissive_TE_delays.population(1:144:end)));
    sd = mean((std(reshape(permissive_TE_delays.TE_X_ZY(logical((ismember(permissive_TE_delays.state,permissive_states)).*(permissive_TE_delays.delay_source==d))),12,size(permissive_states,2)),1)));
    errorbar(d,m,sd,'Color','black','LineWidth',0.5,'CapSize',0)
    scatter(d,m,'filled','blue');

end
hold off
ylim([0 0.05])
yticks([0:0.01:0.05])
xlim([-0.5 11.5])
xticks([0:11])
xticklabels({'0','1','2','3','4','5','6','7','8','9','10','11'})
xlabel('Delay on source variable','FontSize',10,'Color','black')
ylabel({'Weighted conditional'; 'transfer entropy'},'FontSize',10,'Color','black')
% title('Delay analysis for permissive states')

figure
hold on
for d = 0:11
    
    m = (mean(reshape(restrictive_TE_delays.TE_X_ZY(logical((ismember(restrictive_TE_delays.state,restrictive_states)).*(restrictive_TE_delays.delay_source==d))),12,size(restrictive_states,2)),1))*(restrictive_TE_delays.population(1:144:end)/sum(restrictive_TE_delays.population(1:144:end)));
    sd = mean((std(reshape(restrictive_TE_delays.TE_X_ZY(logical((ismember(restrictive_TE_delays.state,restrictive_states)).*(restrictive_TE_delays.delay_source==d))),12,size(restrictive_states,2)),1)));
    errorbar(d,m,sd,'Color','black','LineWidth',0.5,'CapSize',0)
    scatter(d,m,'filled','blue');

end
hold off
ylim([0 0.05])
yticks([0:0.01:0.05])
xlim([-0.5 11.5])
xticks([0:11])
xticklabels({'0','1','2','3','4','5','6','7','8','9','10','11'})
xlabel('Delay on source variable','FontSize',10,'Color','black')
ylabel({'Weighted conditional'; 'transfer entropy'},'FontSize',10,'Color','black')
% title('Delay analysis for restrictive states')