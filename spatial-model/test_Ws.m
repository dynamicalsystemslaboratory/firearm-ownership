% This script explores models that use different spatial weight schemes

clc; clear;

load('data.mat')
load('A.mat')
load('D.mat')
load('S.mat')
load('BC.mat')
load('GDP.mat')
load('B.mat')

%% define data sets

national_data=data((strcmp(data.state, 'United States')),:); % get national-level data

state_data=data(~(strcmp(data.state, 'United States')),:); % remove USA
state_data=state_data(~(strcmp(state_data.state, 'Alaska')),:); % remove AK (no gun ownership data)
state_data=state_data(~(strcmp(state_data.state, 'Hawaii')),:); % remove HI (no BC data)

BC.Alaska = [];
BC.Hawaii = [];

S(strcmp(S.State, 'Alaska'),:) =[];
S(strcmp(S.State, 'Hawaii'),:) =[];

states = unique(state_data.state);
years = unique(state_data.year);
areas = A.total_area_km2(ismember(A.state,states));  


%% arrange monthly data from November until October

% background checks
BC_state_monthly = BC(intersect(find(BC.Year==min(years)-1),find(BC.Month==11)):intersect(find(BC.Year==max(years)),find(BC.Month==10)),:); % background checks in each state, every month
BC_nation_monthly = sum(BC_state_monthly{:,3:end},2); % background checks in the entire country, every month

BCC_nation_monthly = BC_nation_monthly./repelem(national_data.population,12); % convert to background checks per capita (BC -> BCC)
BCC_state_monthly = BC_state_monthly; % will be converted to background checks per capita in the loop below
for s = 1:length(states)
    BCC_state_monthly(:,s+2) = array2table(BCC_state_monthly{:,s+2}./repelem(state_data.population(strcmp(state_data.state,states(s))),12)); % convert to background checks per capita (BC -> BCC)
end 

% suicides
suicides_monthly = reshape(S.FractionOfSuicidesWithFirearms,(length(years)+1)*12,length(states));
suicides_monthly = suicides_monthly(intersect(find(BC.Year==min(years)-1),find(BC.Month==11)):intersect(find(BC.Year==max(years)),find(BC.Month==10)),:);
suicides_monthly = array2table([BC_state_monthly.Year BC_state_monthly.Month suicides_monthly]);
suicides_monthly.Properties.VariableNames = BC_state_monthly.Properties.VariableNames;

suicides_national = [];
for y = transpose(unique(S.Year))
    for m = transpose(unique(S.Month))
        suicides_national = [suicides_national; y m sum(S.SuicidesWithFirearms(intersect(find(S.Year==y),find(S.Month==m))))/sum(S.Suicides(intersect(find(S.Year==y),find(S.Month==m))))];
    end
end

%% define time series for all years

% time series t
fraction_owners_t = state_data.fraction_firearm_owners;
fraction_owners_t(1:(max(years)-min(years)+1):end) = []; 

fraction_suicides_t = reshape(suicides_monthly{:,3:end},length(years)*12*length(states),1);
fraction_suicides_t = fraction_suicides_t(12:12:end);
fraction_suicides_t(1:(max(years)-min(years)+1):end) = [];

BCC_state_t = reshape(BCC_state_monthly{:,3:end},length(years)*12*length(states),1);
BCC_state_t = BCC_state_t(12:12:end);
BCC_state_t(1:(max(years)-min(years)+1):end) = [];

% time series t-12
fraction_owners_t_12 = state_data.fraction_firearm_owners;
fraction_owners_t_12((max(years)-min(years)+1):(max(years)-min(years)+1):end) = [];


% time series t-1
fraction_suicides_t_1 = reshape(suicides_monthly{:,3:end},length(years)*12*length(states),1);
fraction_suicides_t_1 = fraction_suicides_t_1(11:12:end);
fraction_suicides_t_1(1:(max(years)-min(years)+1):end) = [];

BCC_state_t_1 = reshape(BCC_state_monthly{:,3:end},length(years)*12*length(states),1);
BCC_state_t_1 = BCC_state_t_1(11:12:end);
BCC_state_t_1(1:(max(years)-min(years)+1):end) = [];


% create a vector of booleans based on response rate (high/low)
RR = 10; % threshold response rate to surveys
LRR = ismember(state_data.state,unique(state_data.state(state_data.num_respondents<RR))); 
LRR = LRR(ismember(state_data.year,[2000:2018]));
HRR = 1-LRR;

%% prepare time series and matrices for calibration of models

y_t = fraction_owners_t;
y_t_12 = fraction_owners_t_12;

x1_t = BCC_state_t;
x1_t_1 = BCC_state_t_1;    

x2_t = fraction_suicides_t;
x2_t_1 = fraction_suicides_t_1;

% define w and W matrices
w = 1./D{:,:};
w(isinf(w)) = 0; % replace diagonal with 0
w(:,strcmp(D.Properties.RowNames, 'HI')) = []; % remove HI (no BC data)
w(strcmp(D.Properties.RowNames, 'HI'),:) = []; 
w(:,strcmp(D.Properties.RowNames, 'AK')) = []; % remove AK (no BC data)
w(strcmp(D.Properties.RowNames, 'AK'),:) = []; 
w = normw(w); % row-normalize the matrix

% define borders matrix
b = B;
b(:,strcmp(B.Properties.RowNames, 'HI')) = []; % remove HI (no BC data)
b(strcmp(B.Properties.RowNames, 'HI'),:) = []; 
b(:,strcmp(B.Properties.RowNames, 'AK')) = []; % remove AK (no BC data)
b(strcmp(B.Properties.RowNames, 'AK'),:) = []; 

% create a vector of dummy variables to capture a time trend
d = transpose(repmat([1:(length(years)-1)],1,length(states))); 


%% %%%%%%%%%%%%%%%%%%%%%%%%%%% W = 1/distance %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
W = zeros(length(states)*(length(years)-1));
for s = 1:length(states)
    for y = 1:(length(years)-1)
        W((s-1)*(length(years)-1)+y,y:(size(y_t,1)/length(states)):end) = w(s,:);  
    end
end
W = normw(W); % row-normalize the matrix

% run model
results = calibrate_model(y_t,[y_t_12...
    W*y_t_12...
    x1_t_1.*HRR...
    x1_t_1.*LRR...
    x2_t_1.*HRR...
    x2_t_1.*LRR...
    W*(x1_t_1)...
    W*(x2_t_1)...
    HRR...
    LRR...
    d],W);

rho = results.rho; % coefficient for W*y_t
tau = results.beta(1); % coefficient for y_t_12
eta = results.beta(2); % coefficient for W*y_t_12
phi1_hrr = results.beta(3); % coefficient for x1
phi1_lrr = results.beta(4); % coefficient for x1
phi2_hrr = results.beta(5); % coefficient for x2
phi2_lrr = results.beta(6); % coefficient for x2
psi1 = results.beta(7); % coefficient for W*x2
psi2 = results.beta(8); % coefficient for W*x1
alpha_hrr = results.beta(9);
alpha_lrr = results.beta(10);
gamma = results.beta(11);
sig = results.sige;

t_rho = results.tstat(12);
t_tau = results.tstat(1);
t_eta = results.tstat(2);
t_phi1_hrr = results.tstat(3);
t_phi1_lrr = results.tstat(4);
t_phi2_hrr = results.tstat(5);
t_phi2_lrr = results.tstat(6);
t_psi1 = results.tstat(7);
t_psi2 = results.tstat(8);
t_alpha_hrr = results.tstat(9);
t_alpha_lrr = results.tstat(10);
t_gamma = results.tstat(11);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% W = area %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
W = zeros(length(states)*(length(years)-1));
for s = 1:length(states)
    for y = 1:(length(years)-1)
        state_areas = areas;
        state_areas(s)=0;
        W((s-1)*(length(years)-1)+y,y:(size(y_t,1)/length(states)):end) = state_areas;  
    end
end
W = normw(W); % row-normalize the matrix

% run model
results = calibrate_model(y_t,[y_t_12...
    W*y_t_12...
    x1_t_1.*HRR...
    x1_t_1.*LRR...
    x2_t_1.*HRR...
    x2_t_1.*LRR...
    W*(x1_t_1)...
    W*(x2_t_1)...
    HRR...
    LRR...
    d],W);

rho = results.rho; % coefficient for W*y_t
tau = results.beta(1); % coefficient for y_t_12
eta = results.beta(2); % coefficient for W*y_t_12
phi1_hrr = results.beta(3); % coefficient for x1
phi1_lrr = results.beta(4); % coefficient for x1
phi2_hrr = results.beta(5); % coefficient for x2
phi2_lrr = results.beta(6); % coefficient for x2
psi1 = results.beta(7); % coefficient for W*x2
psi2 = results.beta(8); % coefficient for W*x1
alpha_hrr = results.beta(9);
alpha_lrr = results.beta(10);
gamma = results.beta(11);
sig = results.sige;

t_rho = results.tstat(12);
t_tau = results.tstat(1);
t_eta = results.tstat(2);
t_phi1_hrr = results.tstat(3);
t_phi1_lrr = results.tstat(4);
t_phi2_hrr = results.tstat(5);
t_phi2_lrr = results.tstat(6);
t_psi1 = results.tstat(7);
t_psi2 = results.tstat(8);
t_alpha_hrr = results.tstat(9);
t_alpha_lrr = results.tstat(10);
t_gamma = results.tstat(11);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%% W = population %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
W = zeros(length(states)*(length(years)-1));
for s = 1:length(states)
    for y = 1:(length(years)-1)
        population_size = state_data.population(state_data.year==years(y));
        population_size(s) = 0;
        W((s-1)*(length(years)-1)+y,y:(size(y_t,1)/length(states)):end) = (population_size);
    end
end
W = normw(W); % row-normalize the matrix

% run model
results = calibrate_model(y_t,[y_t_12...
    W*y_t_12...
    x1_t_1.*HRR...
    x1_t_1.*LRR...
    x2_t_1.*HRR...
    x2_t_1.*LRR...
    W*(x1_t_1)...
    W*(x2_t_1)...
    HRR...
    LRR...
    d],W);

rho = results.rho; % coefficient for W*y_t
tau = results.beta(1); % coefficient for y_t_12
eta = results.beta(2); % coefficient for W*y_t_12
phi1_hrr = results.beta(3); % coefficient for x1
phi1_lrr = results.beta(4); % coefficient for x1
phi2_hrr = results.beta(5); % coefficient for x2
phi2_lrr = results.beta(6); % coefficient for x2
psi1 = results.beta(7); % coefficient for W*x2
psi2 = results.beta(8); % coefficient for W*x1
alpha_hrr = results.beta(9);
alpha_lrr = results.beta(10);
gamma = results.beta(11);
sig = results.sige;

t_rho = results.tstat(12);
t_tau = results.tstat(1);
t_eta = results.tstat(2);
t_phi1_hrr = results.tstat(3);
t_phi1_lrr = results.tstat(4);
t_phi2_hrr = results.tstat(5);
t_phi2_lrr = results.tstat(6);
t_psi1 = results.tstat(7);
t_psi2 = results.tstat(8);
t_alpha_hrr = results.tstat(9);
t_alpha_lrr = results.tstat(10);
t_gamma = results.tstat(11);



%% %%%%%%%%%%%%%%%%%%%%%%%% W = population/area %%%%%%%%%%%%%%%%%%%%%%%%
W = zeros(length(states)*(length(years)-1));
for s = 1:length(states)
    for y = 1:(length(years)-1)
        population_density = state_data.population(state_data.year==years(y))./areas;
        population_density(s) = 0;
        W((s-1)*(length(years)-1)+y,y:(size(y_t,1)/length(states)):end) = transpose(population_density);  
    end
end
W = normw(W); % row-normalize the matrix

% run model
results = calibrate_model(y_t,[y_t_12...
    W*y_t_12...
    x1_t_1.*HRR...
    x1_t_1.*LRR...
    x2_t_1.*HRR...
    x2_t_1.*LRR...
    W*(x1_t_1)...
    W*(x2_t_1)...
    HRR...
    LRR...
    d],W);

rho = results.rho; % coefficient for W*y_t
tau = results.beta(1); % coefficient for y_t_12
eta = results.beta(2); % coefficient for W*y_t_12
phi1_hrr = results.beta(3); % coefficient for x1
phi1_lrr = results.beta(4); % coefficient for x1
phi2_hrr = results.beta(5); % coefficient for x2
phi2_lrr = results.beta(6); % coefficient for x2
psi1 = results.beta(7); % coefficient for W*x2
psi2 = results.beta(8); % coefficient for W*x1
alpha_hrr = results.beta(9);
alpha_lrr = results.beta(10);
gamma = results.beta(11);
sig = results.sige;

t_rho = results.tstat(12);
t_tau = results.tstat(1);
t_eta = results.tstat(2);
t_phi1_hrr = results.tstat(3);
t_phi1_lrr = results.tstat(4);
t_phi2_hrr = results.tstat(5);
t_phi2_lrr = results.tstat(6);
t_psi1 = results.tstat(7);
t_psi2 = results.tstat(8);
t_alpha_hrr = results.tstat(9);
t_alpha_lrr = results.tstat(10);
t_gamma = results.tstat(11);


%% %%%%%%%%%%%%%%%%%%%% W = population/distance/area %%%%%%%%%%%%%%%%%%%
W = zeros(length(states)*(length(years)-1));
for s = 1:length(states)
    for y = 1:(length(years)-1)
        population_density = state_data.population(state_data.year==years(y))./areas;
        population_density(s) = 0;
        W((s-1)*(length(years)-1)+y,y:(size(y_t,1)/length(states)):end) = w(s,:).*transpose(population_density);  
    end
end
W = normw(W); % row-normalize the matrix

% run model
results = calibrate_model(y_t,[y_t_12...
    W*y_t_12...
    x1_t_1.*HRR...
    x1_t_1.*LRR...
    x2_t_1.*HRR...
    x2_t_1.*LRR...
    W*(x1_t_1)...
    W*(x2_t_1)...
    HRR...
    LRR...
    d],W);

rho = results.rho; % coefficient for W*y_t
tau = results.beta(1); % coefficient for y_t_12
eta = results.beta(2); % coefficient for W*y_t_12
phi1_hrr = results.beta(3); % coefficient for x1
phi1_lrr = results.beta(4); % coefficient for x1
phi2_hrr = results.beta(5); % coefficient for x2
phi2_lrr = results.beta(6); % coefficient for x2
psi1 = results.beta(7); % coefficient for W*x2
psi2 = results.beta(8); % coefficient for W*x1
alpha_hrr = results.beta(9);
alpha_lrr = results.beta(10);
gamma = results.beta(11);
sig = results.sige;

t_rho = results.tstat(12);
t_tau = results.tstat(1);
t_eta = results.tstat(2);
t_phi1_hrr = results.tstat(3);
t_phi1_lrr = results.tstat(4);
t_phi2_hrr = results.tstat(5);
t_phi2_lrr = results.tstat(6);
t_psi1 = results.tstat(7);
t_psi2 = results.tstat(8);
t_alpha_hrr = results.tstat(9);
t_alpha_lrr = results.tstat(10);
t_gamma = results.tstat(11);



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% W = GDP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
W = zeros(length(states)*(length(years)-1));
for s = 1:length(states)
    for y = 1:(length(years)-1)
        gdps = GDP{GDP.Year==years(y),ismember(GDP.Properties.VariableNames,states)};
        gdps(s) = 0;
        W((s-1)*(length(years)-1)+y,y:(size(y_t,1)/length(states)):end) = gdps;  
    end
end
W = normw(W); % row-normalize the matrix

% run model
results = calibrate_model(y_t,[y_t_12...
    W*y_t_12...
    x1_t_1.*HRR...
    x1_t_1.*LRR...
    x2_t_1.*HRR...
    x2_t_1.*LRR...
    W*(x1_t_1)...
    W*(x2_t_1)...
    HRR...
    LRR...
    d],W);

rho = results.rho; % coefficient for W*y_t
tau = results.beta(1); % coefficient for y_t_12
eta = results.beta(2); % coefficient for W*y_t_12
phi1_hrr = results.beta(3); % coefficient for x1
phi1_lrr = results.beta(4); % coefficient for x1
phi2_hrr = results.beta(5); % coefficient for x2
phi2_lrr = results.beta(6); % coefficient for x2
psi1 = results.beta(7); % coefficient for W*x2
psi2 = results.beta(8); % coefficient for W*x1
alpha_hrr = results.beta(9);
alpha_lrr = results.beta(10);
gamma = results.beta(11);
sig = results.sige;

t_rho = results.tstat(12);
t_tau = results.tstat(1);
t_eta = results.tstat(2);
t_phi1_hrr = results.tstat(3);
t_phi1_lrr = results.tstat(4);
t_phi2_hrr = results.tstat(5);
t_phi2_lrr = results.tstat(6);
t_psi1 = results.tstat(7);
t_psi2 = results.tstat(8);
t_alpha_hrr = results.tstat(9);
t_alpha_lrr = results.tstat(10);
t_gamma = results.tstat(11);



%% %%%%%%%%%%%%%%%%%%%%%%%%% W = GDP/distance %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
W = zeros(length(states)*(length(years)-1));
for s = 1:length(states)
    for y = 1:(length(years)-1)
        gdps = GDP{GDP.Year==years(y),ismember(GDP.Properties.VariableNames,states)};
        gdps(s) = 0;
        W((s-1)*(length(years)-1)+y,y:(size(y_t,1)/length(states)):end) = w(s,:).*gdps;  
    end
end
W = normw(W); % row-normalize the matrix

% run model
results = calibrate_model(y_t,[y_t_12...
    W*y_t_12...
    x1_t_1.*HRR...
    x1_t_1.*LRR...
    x2_t_1.*HRR...
    x2_t_1.*LRR...
    W*(x1_t_1)...
    W*(x2_t_1)...
    HRR...
    LRR...
    d],W);

rho = results.rho; % coefficient for W*y_t
tau = results.beta(1); % coefficient for y_t_12
eta = results.beta(2); % coefficient for W*y_t_12
phi1_hrr = results.beta(3); % coefficient for x1
phi1_lrr = results.beta(4); % coefficient for x1
phi2_hrr = results.beta(5); % coefficient for x2
phi2_lrr = results.beta(6); % coefficient for x2
psi1 = results.beta(7); % coefficient for W*x2
psi2 = results.beta(8); % coefficient for W*x1
alpha_hrr = results.beta(9);
alpha_lrr = results.beta(10);
gamma = results.beta(11);
sig = results.sige;

t_rho = results.tstat(12);
t_tau = results.tstat(1);
t_eta = results.tstat(2);
t_phi1_hrr = results.tstat(3);
t_phi1_lrr = results.tstat(4);
t_phi2_hrr = results.tstat(5);
t_phi2_lrr = results.tstat(6);
t_psi1 = results.tstat(7);
t_psi2 = results.tstat(8);
t_alpha_hrr = results.tstat(9);
t_alpha_lrr = results.tstat(10);
t_gamma = results.tstat(11);



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% W = GDP/area %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
W = zeros(length(states)*(length(years)-1));
for s = 1:length(states)
    for y = 1:(length(years)-1)
        gdps = GDP{GDP.Year==years(y),ismember(GDP.Properties.VariableNames,states)};
        gdps(s) = 0;
        W((s-1)*(length(years)-1)+y,y:(size(y_t,1)/length(states)):end) = transpose(gdps)./areas;  
    end
end
W = normw(W); % row-normalize the matrix

% run model
results = calibrate_model(y_t,[y_t_12...
    W*y_t_12...
    x1_t_1.*HRR...
    x1_t_1.*LRR...
    x2_t_1.*HRR...
    x2_t_1.*LRR...
    W*(x1_t_1)...
    W*(x2_t_1)...
    HRR...
    LRR...
    d],W);

rho = results.rho; % coefficient for W*y_t
tau = results.beta(1); % coefficient for y_t_12
eta = results.beta(2); % coefficient for W*y_t_12
phi1_hrr = results.beta(3); % coefficient for x1
phi1_lrr = results.beta(4); % coefficient for x1
phi2_hrr = results.beta(5); % coefficient for x2
phi2_lrr = results.beta(6); % coefficient for x2
psi1 = results.beta(7); % coefficient for W*x2
psi2 = results.beta(8); % coefficient for W*x1
alpha_hrr = results.beta(9);
alpha_lrr = results.beta(10);
gamma = results.beta(11);
sig = results.sige;

t_rho = results.tstat(12);
t_tau = results.tstat(1);
t_eta = results.tstat(2);
t_phi1_hrr = results.tstat(3);
t_phi1_lrr = results.tstat(4);
t_phi2_hrr = results.tstat(5);
t_phi2_lrr = results.tstat(6);
t_psi1 = results.tstat(7);
t_psi2 = results.tstat(8);
t_alpha_hrr = results.tstat(9);
t_alpha_lrr = results.tstat(10);
t_gamma = results.tstat(11);


%% %%%%%%%%%%%%%%%%%%%%%%% W = GDP/area/distance %%%%%%%%%%%%%%%%%%%%%%%%%%
W = zeros(length(states)*(length(years)-1));
for s = 1:length(states)
    for y = 1:(length(years)-1)
        gdps = GDP{GDP.Year==years(y),ismember(GDP.Properties.VariableNames,states)};
        gdps(s) = 0;
        W((s-1)*(length(years)-1)+y,y:(size(y_t,1)/length(states)):end) = transpose(w(s,:)).*(transpose(gdps)./areas);  
    end
end
W = normw(W); % row-normalize the matrix

% run model
results = calibrate_model(y_t,[y_t_12...
    W*y_t_12...
    x1_t_1.*HRR...
    x1_t_1.*LRR...
    x2_t_1.*HRR...
    x2_t_1.*LRR...
    W*(x1_t_1)...
    W*(x2_t_1)...
    HRR...
    LRR...
    d],W);

rho = results.rho; % coefficient for W*y_t
tau = results.beta(1); % coefficient for y_t_12
eta = results.beta(2); % coefficient for W*y_t_12
phi1_hrr = results.beta(3); % coefficient for x1
phi1_lrr = results.beta(4); % coefficient for x1
phi2_hrr = results.beta(5); % coefficient for x2
phi2_lrr = results.beta(6); % coefficient for x2
psi1 = results.beta(7); % coefficient for W*x2
psi2 = results.beta(8); % coefficient for W*x1
alpha_hrr = results.beta(9);
alpha_lrr = results.beta(10);
gamma = results.beta(11);
sig = results.sige;

t_rho = results.tstat(12);
t_tau = results.tstat(1);
t_eta = results.tstat(2);
t_phi1_hrr = results.tstat(3);
t_phi1_lrr = results.tstat(4);
t_phi2_hrr = results.tstat(5);
t_phi2_lrr = results.tstat(6);
t_psi1 = results.tstat(7);
t_psi2 = results.tstat(8);
t_alpha_hrr = results.tstat(9);
t_alpha_lrr = results.tstat(10);
t_gamma = results.tstat(11);


%% %%%%%%%%%%%%%%%%%%%%%%%%% W = border present %%%%%%%%%%%%%%%%%%%%%%%%%%%
W = zeros(length(states)*(length(years)-1));
for s = 1:length(states)
    for y = 1:(length(years)-1)
        W((s-1)*(length(years)-1)+y,y:(size(y_t,1)/length(states)):end) = b{s,:}~=0;  
    end
end
W = normw(W); % row-normalize the matrix


% run model
results = calibrate_model(y_t,[y_t_12...
    W*y_t_12...
    x1_t_1.*HRR...
    x1_t_1.*LRR...
    x2_t_1.*HRR...
    x2_t_1.*LRR...
    W*(x1_t_1)...
    W*(x2_t_1)...
    HRR...
    LRR...
    d],W);

rho = results.rho; % coefficient for W*y_t
tau = results.beta(1); % coefficient for y_t_12
eta = results.beta(2); % coefficient for W*y_t_12
phi1_hrr = results.beta(3); % coefficient for x1
phi1_lrr = results.beta(4); % coefficient for x1
phi2_hrr = results.beta(5); % coefficient for x2
phi2_lrr = results.beta(6); % coefficient for x2
psi1 = results.beta(7); % coefficient for W*x2
psi2 = results.beta(8); % coefficient for W*x1
alpha_hrr = results.beta(9);
alpha_lrr = results.beta(10);
gamma = results.beta(11);
sig = results.sige;

t_rho = results.tstat(12);
t_tau = results.tstat(1);
t_eta = results.tstat(2);
t_phi1_hrr = results.tstat(3);
t_phi1_lrr = results.tstat(4);
t_phi2_hrr = results.tstat(5);
t_phi2_lrr = results.tstat(6);
t_psi1 = results.tstat(7);
t_psi2 = results.tstat(8);
t_alpha_hrr = results.tstat(9);
t_alpha_lrr = results.tstat(10);
t_gamma = results.tstat(11);


%% %%%%%%%%%%%%%%% W = population in states with borders %%%%%%%%%%%%%%%%%%
W = zeros(length(states)*(length(years)-1));
for s = 1:length(states)
    for y = 1:(length(years)-1)
        population_size = state_data.population(state_data.year==years(y));
        W((s-1)*(length(years)-1)+y,y:(size(y_t,1)/length(states)):end) = (b{s,:}~=0).*transpose(population_size);  
    end
end
W = normw(W); % row-normalize the matrix

% run model
results = calibrate_model(y_t,[y_t_12...
    W*y_t_12...
    x1_t_1.*HRR...
    x1_t_1.*LRR...
    x2_t_1.*HRR...
    x2_t_1.*LRR...
    W*(x1_t_1)...
    W*(x2_t_1)...
    HRR...
    LRR...
    d],W);

rho = results.rho; % coefficient for W*y_t
tau = results.beta(1); % coefficient for y_t_12
eta = results.beta(2); % coefficient for W*y_t_12
phi1_hrr = results.beta(3); % coefficient for x1
phi1_lrr = results.beta(4); % coefficient for x1
phi2_hrr = results.beta(5); % coefficient for x2
phi2_lrr = results.beta(6); % coefficient for x2
psi1 = results.beta(7); % coefficient for W*x2
psi2 = results.beta(8); % coefficient for W*x1
alpha_hrr = results.beta(9);
alpha_lrr = results.beta(10);
gamma = results.beta(11);
sig = results.sige;

t_rho = results.tstat(12);
t_tau = results.tstat(1);
t_eta = results.tstat(2);
t_phi1_hrr = results.tstat(3);
t_phi1_lrr = results.tstat(4);
t_phi2_hrr = results.tstat(5);
t_phi2_lrr = results.tstat(6);
t_psi1 = results.tstat(7);
t_psi2 = results.tstat(8);
t_alpha_hrr = results.tstat(9);
t_alpha_lrr = results.tstat(10);
t_gamma = results.tstat(11);


%% %%%%%%%%%%%%%%%%%%% W = GDP in states with borders %%%%%%%%%%%%%%%%%%%%%
W = zeros(length(states)*(length(years)-1));
for s = 1:length(states)
    for y = 1:(length(years)-1)
        gdps = GDP{GDP.Year==years(y),ismember(GDP.Properties.VariableNames,states)};
        gdps(s) = 0;
        W((s-1)*(length(years)-1)+y,y:(size(y_t,1)/length(states)):end) = (b{s,:}~=0).*transpose(population_size);  
    end
end
W = normw(W); % row-normalize the matrix

% run model
results = calibrate_model(y_t,[y_t_12...
    W*y_t_12...
    x1_t_1.*HRR...
    x1_t_1.*LRR...
    x2_t_1.*HRR...
    x2_t_1.*LRR...
    W*(x1_t_1)...
    W*(x2_t_1)...
    HRR...
    LRR...
    d],W);

rho = results.rho; % coefficient for W*y_t
tau = results.beta(1); % coefficient for y_t_12
eta = results.beta(2); % coefficient for W*y_t_12
phi1_hrr = results.beta(3); % coefficient for x1
phi1_lrr = results.beta(4); % coefficient for x1
phi2_hrr = results.beta(5); % coefficient for x2
phi2_lrr = results.beta(6); % coefficient for x2
psi1 = results.beta(7); % coefficient for W*x2
psi2 = results.beta(8); % coefficient for W*x1
alpha_hrr = results.beta(9);
alpha_lrr = results.beta(10);
gamma = results.beta(11);
sig = results.sige;

t_rho = results.tstat(12);
t_tau = results.tstat(1);
t_eta = results.tstat(2);
t_phi1_hrr = results.tstat(3);
t_phi1_lrr = results.tstat(4);
t_phi2_hrr = results.tstat(5);
t_phi2_lrr = results.tstat(6);
t_psi1 = results.tstat(7);
t_psi2 = results.tstat(8);
t_alpha_hrr = results.tstat(9);
t_alpha_lrr = results.tstat(10);
t_gamma = results.tstat(11);

