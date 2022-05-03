function [I,perms,percentile,p_value] = test_contemporaneous(X,Y,Z,n_iter)

Xt = X(1:end-1); % X at time step t
Xt1 = X(2:end); % X at time step t+1
Yt = Y(1:end-1); % Y at time step t
Yt1 = Y(2:end); % Y at time step t+1
Zt = Z(1:end-1); % Z at time step t

I = compute_entropy([Xt1 Xt Yt Zt])-compute_entropy([Xt Yt Zt])+compute_entropy([Yt1 Xt Yt Zt])-compute_entropy([Xt1 Yt1 Xt Yt Zt]);

perms = [];
for it = 1:n_iter

    sX = X(randperm(length(X))); % shuffle time series X
    sY = Y(randperm(length(Y))); % shuffle time series Y
    sXt = sX(1:end-1); % shuffled time series of X at time step t
    sXt1 = sX(2:end); % shuffled time series of X at time step t+1
    sYt = sY(1:end-1); % shuffled time series of Y at time step t
    sYt1 = sY(2:end); % shuffled time series of Y at time step t+1
    
    sI = compute_entropy([sXt1 sXt sYt Zt])-compute_entropy([sXt sYt Zt])+compute_entropy([sYt1 sXt sYt Zt])-compute_entropy([sXt1 sYt1 sXt sYt Zt]);
    perms = [perms; sI];
    
end
percentile = prctile(perms,95); % compute the 95th percentile
prcnt = [1-transpose([1:0.0001:100])/100 prctile(perms,[1:0.0001:100],'all')];
p_value = prcnt(find(abs(prcnt(:,2)-I)==min(abs(prcnt(:,2)-I)),1,'first'),1); % compute the p-value


end