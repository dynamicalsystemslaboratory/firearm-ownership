function [TE,perm_tests] = compute_conditional_TE(ts,n_iter)

n_vars = size(ts,2); % number of variables
var_pairs = [nchoosek([1:n_vars],2); nchoosek([n_vars:-1:1],2)]; % possible pairs of variables

% set up tables for results
TE = array2table(nan(size(var_pairs,1),6),'VariableNames',[{'source'} {'sink'} {'sourcesink'} {'TE'} {'percentile_95'} {'p_value'}]);
TE.source = transpose(string(ts.Properties.VariableNames(var_pairs(:,1))));
TE.sink = transpose(string(ts.Properties.VariableNames(var_pairs(:,2))));
TE.sourcesink = strcat(TE.source,"-",TE.sink);
perm_tests = array2table(nan(n_iter,size(var_pairs,1)),'VariableNames',TE.sourcesink); % vector to be filled with n_iter values of surrogate TE
for vp = 1:size(var_pairs,1) % iterate through variable pairs
    
    source_ts = ts{:,strcmp(ts.Properties.VariableNames,TE.source(vp))}; % time series of source
    sink_ts = ts{:,strcmp(ts.Properties.VariableNames,TE.sink(vp))}; % time series of sink
    conditioned_ts = ts{:,~(strcmp(ts.Properties.VariableNames,TE.sink(vp))+strcmp(ts.Properties.VariableNames,TE.source(vp)))}; % time series of variables to be conditioned upon
    
    TE.TE(vp) = compute_entropy([sink_ts(2:end) sink_ts(1:end-1) conditioned_ts(1:end-1,:)])-...
    compute_entropy([sink_ts(2:end) sink_ts(1:end-1) conditioned_ts(1:end-1,:) source_ts(1:end-1)])-...
    compute_entropy([sink_ts(1:end-1) conditioned_ts(1:end-1,:)])+...
    compute_entropy([sink_ts(1:end-1) conditioned_ts(1:end-1,:) source_ts(1:end-1)]); % compute observed TE

    for it = 1:n_iter
        
        cond_comb = unique([sink_ts conditioned_ts],'rows'); % possible combinations that need to be preserved upon shuffling
        shuffled_source_ts = source_ts; % source time series to be shuffled
        for c_v1v2 = 1:size(cond_comb,1)
            shuffle_subset = shuffled_source_ts(ismember([sink_ts conditioned_ts], cond_comb(c_v1v2,:), 'rows'));
            shuffled_source_ts(ismember([sink_ts conditioned_ts], cond_comb(c_v1v2,:), 'rows')) = shuffle_subset(randperm(length(shuffle_subset)));
        end
        
        perm_tests{it,vp} = compute_entropy([sink_ts(2:end) sink_ts(1:end-1) conditioned_ts(1:end-1,:)])-...
        compute_entropy([sink_ts(2:end) sink_ts(1:end-1) conditioned_ts(1:end-1,:) shuffled_source_ts(1:end-1)])-...
        compute_entropy([sink_ts(1:end-1) conditioned_ts(1:end-1,:)])+...
        compute_entropy([sink_ts(1:end-1) conditioned_ts(1:end-1,:) shuffled_source_ts(1:end-1)]); % compute surrogate TE
    
    end 

    TE.percentile_95(vp) = prctile(perm_tests{:,vp},95); % compute the 95th percentile
    prcnt = [1-transpose([1:0.0001:100])/100 prctile(perm_tests{:,vp},[1:0.0001:100],'all')];
    TE.p_value(vp) = prcnt(find(abs(prcnt(:,2)-TE.TE(vp))==min(abs(prcnt(:,2)-TE.TE(vp))),1,'first'),1);
end
    
    
    
end