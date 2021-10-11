clc; clear;

load('data.mat')
load('D.mat')
load('S.mat')
load('BC.mat')

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

%% estimate model parameters

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

W = zeros(length(states)*(length(years)-1));
for s = 1:length(states)
    for y = 1:(length(years)-1)
        population_fraction = state_data.population(state_data.year==years(y))/sum(state_data.population(state_data.year==years(y)));
        W((s-1)*(length(years)-1)+y,y:(size(y_t,1)/length(states)):end) = w(s,:).*transpose(population_fraction);  
    end
end
W = normw(W); % row-normalize the matrix

% create a vector of dummy variables to capture a time trend
d = transpose(repmat([1:(length(years)-1)],1,length(states))); 

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

%% infer gun ownership on a monthly resolution

firearm_ownership_monthly = array2table([[1998; 1998; transpose(repelem(1999,12))] [11; 12; transpose(1:12)] repmat(transpose(state_data.fraction_firearm_owners(state_data.year==2000)),14,1); ]);
firearm_ownership_monthly.Properties.VariableNames = [{'Year' 'Month'} transpose(states)];

% infer ownership
for obs = 3:size(BCC_state_monthly,1)
    y = BCC_state_monthly.Year(obs);
    m = BCC_state_monthly.Month(obs);

    % define the variables for the year/month
    if m == 10 % we have a measurement of the ownership from surveys every October
        if y == 2000 % if we do not have a measurement of ownership because this is the first year of survey data
             firearm_ownership_monthly(end+1,:) = array2table([y m transpose(state_data.fraction_firearm_owners(state_data.year==y))]);
        else
            y_t_12 = state_data.fraction_firearm_owners(state_data.year==y-1); % get ownership from previous year          
            x1_t_1 = transpose(BCC_state_monthly{intersect(find(BCC_state_monthly.Year==y-1),find(BCC_state_monthly.Month==12)),3:end}); % get background checks data from previous month
            x2_t_1 = transpose(suicides_monthly{intersect(find(suicides_monthly.Year==y-1),find(suicides_monthly.Month==12)),3:end}); % get suicides data from previous month
            firearm_ownership_monthly(end+1,:) = array2table([y m transpose(inv(eye(size(w))-rho*w)*(...
            tau*y_t_12 + ...
            eta*w*y_t_12 + ...
            phi1_hrr*x1_t_1.*HRR(1:19:end) + ...
            phi1_lrr*x1_t_1.*LRR(1:19:end) + ...
            phi2_hrr*x2_t_1.*HRR(1:19:end) + ...
            phi2_lrr*x2_t_1.*LRR(1:19:end) + ...
            psi1*w*(x1_t_1) + ...
            psi2*w*(x2_t_1) + ...
            alpha_hrr*HRR(1:19:end) + ...
            alpha_lrr*LRR(1:19:end) + ...
            repelem(gamma*(y-min(years)+1),length(states),1)))]); % infer ownershop using the model
        end

    elseif m == 1 % data for the month before is from the previous year
        y_t_12 = transpose(firearm_ownership_monthly{intersect(find(firearm_ownership_monthly.Year==y-1),find(firearm_ownership_monthly.Month==m)),3:end}); % get ownership from previous year's estimate
        x1_t_1 = transpose(BCC_state_monthly{intersect(find(BCC_state_monthly.Year==y-1),find(BCC_state_monthly.Month==12)),3:end}); % get background checks data from previous month
        x2_t_1 = transpose(suicides_monthly{intersect(find(suicides_monthly.Year==y-1),find(suicides_monthly.Month==12)),3:end}); % get suicides data from previous month
        firearm_ownership_monthly(end+1,:) = array2table([y m transpose(inv(eye(size(w))-rho*w)*(...
            tau*y_t_12 + ...
            eta*w*y_t_12 + ...
            phi1_hrr*x1_t_1.*HRR(1:19:end) + ...
            phi1_lrr*x1_t_1.*LRR(1:19:end) + ...
            phi2_hrr*x2_t_1.*HRR(1:19:end) + ...
            phi2_lrr*x2_t_1.*LRR(1:19:end) + ...
            psi1*w*(x1_t_1) + ...
            psi2*w*(x2_t_1) + ...
            alpha_hrr*HRR(1:19:end) + ...
            alpha_lrr*LRR(1:19:end) + ...
            repelem(gamma*(y-min(years)+1),length(states),1)))]); % infer ownershop using the model

    else % data for the month before is from the same year
        y_t_12 = transpose(firearm_ownership_monthly{intersect(find(firearm_ownership_monthly.Year==y-1),find(firearm_ownership_monthly.Month==m)),3:end}); % get ownership from previous year's estimate
        x1_t_1 = transpose(BCC_state_monthly{intersect(find(BCC_state_monthly.Year==y),find(BCC_state_monthly.Month==m-1)),3:end}); % get background checks data from previous month
        x2_t_1 = transpose(suicides_monthly{intersect(find(suicides_monthly.Year==y),find(suicides_monthly.Month==m-1)),3:end}); % get suicides data from previous month
        firearm_ownership_monthly(end+1,:) = array2table([y m transpose(inv(eye(size(w))-rho*w)*(...
            tau*y_t_12 + ...
            eta*w*y_t_12 + ...
            phi1_hrr*x1_t_1.*HRR(1:19:end) + ...
            phi1_lrr*x1_t_1.*LRR(1:19:end) + ...
            phi2_hrr*x2_t_1.*HRR(1:19:end) + ...
            phi2_lrr*x2_t_1.*LRR(1:19:end) + ...
            psi1*w*(x1_t_1) + ...
            psi2*w*(x2_t_1) + ...
            alpha_hrr*HRR(1:19:end) + ...
            alpha_lrr*LRR(1:19:end) + ...
            repelem(gamma*(y-min(years)+1),length(states),1)))]); % infer ownershop using the model
    end
end


firearm_ownership_monthly = firearm_ownership_monthly(ismember(firearm_ownership_monthly.Year,years),:); % remove primer points for initial estimations

% compute national-level firearm ownership
for p = 1:size(firearm_ownership_monthly,1)
    firearm_ownership_monthly.USA(p) = sum(transpose(firearm_ownership_monthly{p,3:(length(states)+2)}).*state_data.population(state_data.year==firearm_ownership_monthly.Year(p)))./sum(state_data.population(state_data.year==firearm_ownership_monthly.Year(p)));
end

figure % reproduce figure S4
t=tiledlayout(7,7);
for s = 1:length(states)   
    nexttile
    plot(2000:(1/12):(max(years)+9/12),firearm_ownership_monthly{:,s+2},'Color','black')
    ylim([0 1])
    hold on 
    sp = scatter((2000+10/12):1:2020,state_data.fraction_firearm_owners((strcmp(state_data.state,states(s)))),20,'filled')
    sp.MarkerEdgeColor = 'black';
    sp.MarkerFaceColor = 'red';
    title(string(states(s)));
    hold off    
end
xlabel(t,'Time','FontSize',16,'Color','black')
ylabel(t,'Fraction of gun owners','FontSize',16,'Color','black')

%% compute sum of square errors (SSE) and mean of square errors (MSE)

errors = array2table(nan(length(states)+1,3),'VariableNames',{'state' 'sse' 'mse'});
errors.state = ['USA'; states];
errors.sse(1) = sum((firearm_ownership_monthly.USA(firearm_ownership_monthly.Month==10)-national_data.fraction_firearm_owners).^2);
errors.mse(1) = sum((firearm_ownership_monthly.USA(firearm_ownership_monthly.Month==10)-national_data.fraction_firearm_owners).^2)/length(national_data.fraction_firearm_owners);
for s = transpose(states)
    
    x = firearm_ownership_monthly{firearm_ownership_monthly.Month==10,strcmp(firearm_ownership_monthly.Properties.VariableNames,s)};
    y = state_data.fraction_firearm_owners(strcmp(state_data.state,s));
    e = (x-y).^2;
    
    errors.sse(strcmp(errors.state,s)) = sum(e);
    errors.mse(strcmp(errors.state,s)) = sum(e)/length(e);
    
end

sse = sum((state_data.fraction_firearm_owners-reshape(firearm_ownership_monthly{firearm_ownership_monthly.Month==10,3:50},960,1).^2));
mse = sum((state_data.fraction_firearm_owners-reshape(firearm_ownership_monthly{firearm_ownership_monthly.Month==10,3:50},960,1).^2))/960;
