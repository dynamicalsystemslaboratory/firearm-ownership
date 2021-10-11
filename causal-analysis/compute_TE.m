function [TE,distributions,percentiles_95,p_values] = compute_TE(XYZX1Y1Z1,num_iterations)

[TE_X_YZ,TE_X_ZY,TE_Y_XZ,TE_Y_ZX,TE_Z_XY,TE_Z_YX] = conditional_TE(XYZX1Y1Z1);
TE = [TE_X_YZ,TE_X_ZY,TE_Y_XZ,TE_Y_ZX,TE_Z_XY,TE_Z_YX];

distributions = array2table(nan(num_iterations,6),'VariableNames',{'TE_X_YZ','TE_X_ZY','TE_Y_XZ','TE_Y_ZX','TE_Z_XY','TE_Z_YX'}); % empty table to be filled with surrogate distributions
for iter = 1:num_iterations
    
    % shuffle X with respect to Y and Z (compute TE X->Y|Z and X->Z|Y)
    shuffled_ts = XYZX1Y1Z1(:,1:3);
    for y = [0 1]
        for z = [0 1]
            ind = intersect(find(shuffled_ts.Y == y),find(shuffled_ts.Z == z));
            shuffled_ts.X(ind) = shuffled_ts.X(ind(randperm(length(ind))));            
        end 
    end
    [TE_X_YZ,TE_X_ZY,~,~,~,~] = conditional_TE(array2table([shuffled_ts{1:end-1,:} shuffled_ts{2:end,:}],'VariableNames',{'X','Y','Z','X1','Y1','Z1'})); % time series for t and t+1

    % shuffle Y with respect to X and Z (compute TE Y->X|Z and TE Y->Z|X)
    shuffled_ts = XYZX1Y1Z1(:,1:3);
    for x = [0 1]
        for z = [0 1]
            ind = intersect(find(shuffled_ts.X == x),find(shuffled_ts.Z == z));
            shuffled_ts.Y(ind) = shuffled_ts.Y(ind(randperm(length(ind))));            
        end 
    end
    [~,~,TE_Y_XZ,TE_Y_ZX,~,~] = conditional_TE(array2table([shuffled_ts{1:end-1,:} shuffled_ts{2:end,:}],'VariableNames',{'X','Y','Z','X1','Y1','Z1'})); % time series for t and t+1

    % shuffle Z with respect to X and Y (compute TE Z->X|Y and Z->Y|X)
    shuffled_ts = XYZX1Y1Z1(:,1:3);
    for x = [0 1]
        for y = [0 1]
            ind = intersect(find(shuffled_ts.X == x),find(shuffled_ts.Y == y));
            shuffled_ts.Z(ind) = shuffled_ts.Z(ind(randperm(length(ind))));            
        end 
    end
    [~,~,~,~,TE_Z_XY,TE_Z_YX] = conditional_TE(array2table([shuffled_ts{1:end-1,:} shuffled_ts{2:end,:}],'VariableNames',{'X','Y','Z','X1','Y1','Z1'})); % time series for t and t+1

    distributions(iter,:) = array2table([TE_X_YZ,TE_X_ZY,TE_Y_XZ,TE_Y_ZX,TE_Z_XY,TE_Z_YX]);
end

percentiles_95 = [quantile(distributions.TE_X_YZ,0.95) quantile(distributions.TE_X_ZY,0.95) quantile(distributions.TE_Y_XZ,0.95) quantile(distributions.TE_Y_ZX,0.95) quantile(distributions.TE_Z_XY,0.95) quantile(distributions.TE_Z_YX,0.95)];

p_values = [];
for a = 1:6
    prcnt = [1-transpose([1:0.001:100])/100 prctile(distributions{:,a},[1:0.001:100],'all')];
    p_values = [p_values prcnt(find(abs(prcnt(:,2)-TE(a))==min(abs(prcnt(:,2)-TE(a))),1,'first'),1)];
end


end