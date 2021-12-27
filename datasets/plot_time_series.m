%% Plotting the raw time series used in transfer entropy analysis

load('nature.mat');
load('BC.mat');
load('BCC.mat');
load('S.mat');
load('FO.mat');
load('BC_2000_2017_sa_dt.mat');
load('BCC_2000_2017_sa_dt.mat');
load('S_2000_2017_sa_dt.mat');
load('FO_2000_2017_sa_dt.mat');

%%%%%%%%%%%%%%%% plot raw time series on a national level %%%%%%%%%%%%%%%%%

% plot background checks
figure
plot(2000:(1/12):2017+(11/12),BC.USA(ismember(BC.Year,[2000:2017])),'Color','black','LineWidth',1)
ylabel({'Background checks'})
ylim([0,3000000])
set(gca,'box','off')
set(gca, 'TickDir', 'out')

% plot background checks per capita
figure
plot(2000:(1/12):2017+(11/12),BCC.USA(ismember(BCC.Year,[2000:2017])),'Color','black','LineWidth',1)
ylabel({'Background checks'; 'per capita'})
ylim([0,0.008])
set(gca,'box','off')
set(gca, 'TickDir', 'out')

% plot suicides
figure
plot(2000:(1/12):2017+(11/12),S.USA(ismember(S.Year,[2000:2017])),'Color','black','LineWidth',1)
ylabel({'Fraction of suicides'; 'committed with firearms'})
ylim([0,1])
set(gca,'box','off')
set(gca, 'TickDir', 'out')

% plot firearm ownership
figure
plot(2000:(1/12):2017+(11/12),FO.USA(ismember(FO.Year,[2000:2017])),'Color','black','LineWidth',1)
ylabel({'Fraction of'; 'firearm owners'})
ylim([0,1])
set(gca,'box','off')
set(gca, 'TickDir', 'out')

% plot mass shootings
figure
bar(2000:(1/12):2017+(11/12),nature.Mass_shooting(ismember(nature.Year,[2000:2017])),'FaceColor','black','LineWidth',0.1)
ylabel('Mass shootings')
ylim([0,5])
set(gca,'box','off')
set(gca, 'TickDir', 'out')

% plot media output
figure
mo_ts = log10(nature.Firearm_laws_and_regulations(ismember(nature.Year,[2000:2017])));
mo_ts(69) = 0;
plot(2000:(1/12):2017+(11/12),mo_ts,'Color','black','LineWidth',1)
yticks([0:4])
yticklabels({'0','10^1','10^2','10^3','10^4'})
ylabel({'Media output on firearm'; 'laws and regulations'})
set(gca,'box','off')
set(gca, 'TickDir', 'out')


%%%%%%%%%%%%% plot processed time series on a national level %%%%%%%%%%%%%%

% plot sa/dt background checks
figure
plot(2000:(1/12):2017+(11/12),BC_2000_2017_sa_dt.USA(ismember(BC_2000_2017_sa_dt.Year,[2000:2017])),'Color','black','LineWidth',1)
ylabel({'Background checks'})
ylim([-200000,1000000])
set(gca,'box','off')
set(gca, 'TickDir', 'out')

% plot sa/dt background checks per capita
figure
plot(2000:(1/12):2017+(11/12),BCC_2000_2017_sa_dt.USA(ismember(BCC_2000_2017_sa_dt.Year,[2000:2017])),'Color','black','LineWidth',1)
ylabel({'Background checks'; 'per capita'})
ylim([-0.001,0.003])
set(gca,'box','off')
set(gca, 'TickDir', 'out')

% plot sa/dt suicides
figure
plot(2000:(1/12):2017+(11/12),S_2000_2017_sa_dt.USA(ismember(S_2000_2017_sa_dt.Year,[2000:2017])),'Color','black','LineWidth',1)
ylabel({'Fraction of suicides'; 'committed with firearms'})
ylim([-0.1,0.1])
set(gca,'box','off')
set(gca, 'TickDir', 'out')

% plot firearm ownership
figure
plot(2000:(1/12):2017+(11/12),FO_2000_2017_sa_dt.USA(ismember(FO_2000_2017_sa_dt.Year,[2000:2017])),'Color','black','LineWidth',1)
ylabel({'Fraction of'; 'firearm owners'})
ylim([-0.3,0.3])
yticks([-0.3:0.15:0.3])
set(gca,'box','off')
set(gca, 'TickDir', 'out')

