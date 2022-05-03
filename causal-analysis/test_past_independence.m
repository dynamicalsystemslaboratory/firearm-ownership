function [I,perms,percentile,p_value] = test_past_independence(ts,n_iter)

T = ts(1:end-3); % ts at time step t
T1 = ts(2:end-2); % ts at time step t+1 
T2 = ts(3:end-1); % ts at time step t+2 
T3 = ts(4:end); % ts at time step t+3

I = compute_entropy([T3 T2 T])-compute_entropy([T2 T])-compute_entropy([T3 T2 T1 T])+compute_entropy([T2 T1 T]);

perms = [];
for it = 1:n_iter

    % shuffle time series 
    sT = T(randperm(length(T))); 
    sT1 = T1(randperm(length(T1)));   
    
    % compute I with shuffled time series
    sI = compute_entropy([T3 T2 T])-compute_entropy([T2 T])-compute_entropy([T3 T2 sT1 T])+compute_entropy([T2 sT1 T]);
    perms = [perms; sI];
    
end
percentile = prctile(perms,95); % compute the 95th percentile
prcnt = [1-transpose([1:0.0001:100])/100 prctile(perms,[1:0.0001:100],'all')];
p_value = prcnt(find(abs(prcnt(:,2)-I)==min(abs(prcnt(:,2)-I)),1,'first'),1); % compute the p-value



end