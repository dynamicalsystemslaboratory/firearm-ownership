% this function computes conditional transfer entropy between three time
% series X, Y, and Z of length N, given a matrix [X Y Z].

% Variable names X, Y, and Z denote time series X_t, Y_t, and Z_t
% Variable names X1, Y1, and Z1 denote time series X_t+1, Y_t+1, and Z_t+1

function [TE_X_YZ,TE_X_ZY,TE_Y_XZ,TE_Y_ZX,TE_Z_XY,TE_Z_YX] = conditional_TE(XYZX1Y1Z1)

% compute joint probabilities

[comb_y1yz, ~, ind_comb_y1yz] = unique([XYZX1Y1Z1.Y1 XYZX1Y1Z1.Y XYZX1Y1Z1.Z],'rows');
P_Y1YZ = [comb_y1yz accumarray(ind_comb_y1yz, 1)/size(XYZX1Y1Z1,1)];

[comb_z1yz, ~, ind_comb_z1yz] = unique([XYZX1Y1Z1.Z1 XYZX1Y1Z1.Y XYZX1Y1Z1.Z],'rows');
P_Z1YZ = [comb_z1yz accumarray(ind_comb_z1yz, 1)/size(XYZX1Y1Z1,1)];

[comb_x1xz, ~, ind_comb_x1xz] = unique([XYZX1Y1Z1.X1 XYZX1Y1Z1.X XYZX1Y1Z1.Z],'rows');
P_X1XZ = [comb_x1xz accumarray(ind_comb_x1xz, 1)/size(XYZX1Y1Z1,1)];

[comb_z1xz, ~, ind_comb_z1xz] = unique([XYZX1Y1Z1.Z1 XYZX1Y1Z1.X XYZX1Y1Z1.Z],'rows');
P_Z1XZ = [comb_z1xz accumarray(ind_comb_z1xz, 1)/size(XYZX1Y1Z1,1)];

[comb_x1xy, ~, ind_comb_x1xy] = unique([XYZX1Y1Z1.X1 XYZX1Y1Z1.X XYZX1Y1Z1.Y],'rows');
P_X1XY = [comb_x1xy accumarray(ind_comb_x1xy, 1)/size(XYZX1Y1Z1,1)];

[comb_y1xy, ~, ind_comb_y1xy] = unique([XYZX1Y1Z1.Y1 XYZX1Y1Z1.X XYZX1Y1Z1.Y],'rows');
P_Y1XY = [comb_y1xy accumarray(ind_comb_y1xy, 1)/size(XYZX1Y1Z1,1)];

[comb_x1xyz, ~, ind_comb_x1xyz] = unique([XYZX1Y1Z1.X1 XYZX1Y1Z1.X XYZX1Y1Z1.Y XYZX1Y1Z1.Z],'rows');
P_X1XYZ = [comb_x1xyz accumarray(ind_comb_x1xyz, 1)/size(XYZX1Y1Z1,1)];

[comb_y1xyz, ~, ind_comb_y1xyz] = unique([XYZX1Y1Z1.Y1 XYZX1Y1Z1.X XYZX1Y1Z1.Y XYZX1Y1Z1.Z],'rows');
P_Y1XYZ = [comb_y1xyz accumarray(ind_comb_y1xyz, 1)/size(XYZX1Y1Z1,1)];

[comb_z1xyz, ~, ind_comb_z1xyz] = unique([XYZX1Y1Z1.Z1 XYZX1Y1Z1.X XYZX1Y1Z1.Y XYZX1Y1Z1.Z],'rows');
P_Z1XYZ = [comb_z1xyz accumarray(ind_comb_z1xyz, 1)/size(XYZX1Y1Z1,1)];


% compute conditional probabilities

[comb_y1_yz, ~, ind_comb_y1_yz] = unique([XYZX1Y1Z1.Y1 XYZX1Y1Z1.Y XYZX1Y1Z1.Z],'rows');
[~, ~, ind_comb_yz] = unique(comb_y1_yz(:,2:end),'rows');
denominator = accumarray(ind_comb_yz,accumarray(ind_comb_y1_yz, 1));
P_Y1_YZ = [comb_y1_yz accumarray(ind_comb_y1_yz, 1)./denominator(ind_comb_yz)];

[comb_z1_yz, ~, ind_comb_z1_yz] = unique([XYZX1Y1Z1.Z1 XYZX1Y1Z1.Y XYZX1Y1Z1.Z],'rows');
[~, ~, ind_comb_yz] = unique(comb_z1_yz(:,2:end),'rows');
denominator = accumarray(ind_comb_yz,accumarray(ind_comb_z1_yz, 1));
P_Z1_YZ = [comb_z1_yz accumarray(ind_comb_z1_yz, 1)./denominator(ind_comb_yz)];

[comb_x1_xz, ~, ind_comb_x1_xz] = unique([XYZX1Y1Z1.X1 XYZX1Y1Z1.X XYZX1Y1Z1.Z],'rows');
[~, ~, ind_comb_xz] = unique(comb_x1_xz(:,2:end),'rows');
denominator = accumarray(ind_comb_xz,accumarray(ind_comb_x1_xz, 1));
P_X1_XZ = [comb_x1_xz accumarray(ind_comb_x1_xz, 1)./denominator(ind_comb_xz)];

[comb_z1_xz, ~, ind_comb_z1_xz] = unique([XYZX1Y1Z1.Z1 XYZX1Y1Z1.X XYZX1Y1Z1.Z],'rows');
[~, ~, ind_comb_xz] = unique(comb_z1_xz(:,2:end),'rows');
denominator = accumarray(ind_comb_xz,accumarray(ind_comb_z1_xz, 1));
P_Z1_XZ = [comb_z1_xz accumarray(ind_comb_z1_xz, 1)./denominator(ind_comb_xz)];

[comb_x1_xy, ~, ind_comb_x1_xy] = unique([XYZX1Y1Z1.X1 XYZX1Y1Z1.X XYZX1Y1Z1.Y],'rows');
[~, ~, ind_comb_xy] = unique(comb_x1_xy(:,2:end),'rows');
denominator = accumarray(ind_comb_xy,accumarray(ind_comb_x1_xy, 1));
P_X1_XY = [comb_x1_xy accumarray(ind_comb_x1_xy, 1)./denominator(ind_comb_xy)];

[comb_y1_xy, ~, ind_comb_y1_xy] = unique([XYZX1Y1Z1.Y1 XYZX1Y1Z1.X XYZX1Y1Z1.Y],'rows');
[~, ~, ind_comb_xy] = unique(comb_y1_xy(:,2:end),'rows');
denominator = accumarray(ind_comb_xy,accumarray(ind_comb_y1_xy, 1));
P_Y1_XY = [comb_y1_xy accumarray(ind_comb_y1_xy, 1)./denominator(ind_comb_xy)];

[comb_x1_xyz, ~, ind_comb_x1_xyz] = unique([XYZX1Y1Z1.X1 XYZX1Y1Z1.X XYZX1Y1Z1.Y XYZX1Y1Z1.Z],'rows');
[~, ~, ind_comb_xyz] = unique(comb_x1_xyz(:,2:end),'rows');
denominator = accumarray(ind_comb_xyz,accumarray(ind_comb_x1_xyz, 1));
P_X1_XYZ = [comb_x1_xyz accumarray(ind_comb_x1_xyz, 1)./denominator(ind_comb_xyz)];

[comb_y1_xyz, ~, ind_comb_y1_xyz] = unique([XYZX1Y1Z1.Y1 XYZX1Y1Z1.X XYZX1Y1Z1.Y XYZX1Y1Z1.Z],'rows');
[~, ~, ind_comb_xyz] = unique(comb_y1_xyz(:,2:end),'rows');
denominator = accumarray(ind_comb_xyz,accumarray(ind_comb_y1_xyz, 1));
P_Y1_XYZ = [comb_y1_xyz accumarray(ind_comb_y1_xyz, 1)./denominator(ind_comb_xyz)];

[comb_z1_xyz, ~, ind_comb_z1_xyz] = unique([XYZX1Y1Z1.Z1 XYZX1Y1Z1.X XYZX1Y1Z1.Y XYZX1Y1Z1.Z],'rows');
[~, ~, ind_comb_xyz] = unique(comb_z1_xyz(:,2:end),'rows');
denominator = accumarray(ind_comb_xyz,accumarray(ind_comb_z1_xyz, 1));
P_Z1_XYZ = [comb_z1_xyz accumarray(ind_comb_z1_xyz, 1)./denominator(ind_comb_xyz)];


% % compute entropies
H_Y1_YZ = 0;
for r = 1:size(P_Y1YZ,1)
    if ~isnan(P_Y1YZ(r,end)*log2(P_Y1_YZ(all(P_Y1YZ(r,1:end-1)==P_Y1_YZ(:,1:end-1),2),end)))
        H_Y1_YZ = H_Y1_YZ - P_Y1YZ(r,end)*log2(P_Y1_YZ(all(P_Y1YZ(r,1:end-1)==P_Y1_YZ(:,1:end-1),2),end));
    else
        continue
    end
end

H_Z1_YZ = 0;
for r = 1:size(P_Z1YZ,1)
    if ~isnan(P_Z1YZ(r,end)*log2(P_Z1_YZ(all(P_Z1YZ(r,1:end-1)==P_Z1_YZ(:,1:end-1),2),end)))    
        H_Z1_YZ = H_Z1_YZ - P_Z1YZ(r,end)*log2(P_Z1_YZ(all(P_Z1YZ(r,1:end-1)==P_Z1_YZ(:,1:end-1),2),end));
    else
        continue
    end
end

H_X1_XZ = 0;
for r = 1:size(P_X1XZ,1)
    if ~isnan(P_X1XZ(r,end)*log2(P_X1_XZ(all(P_X1XZ(r,1:end-1)==P_X1_XZ(:,1:end-1),2),end)))
        H_X1_XZ = H_X1_XZ - P_X1XZ(r,end)*log2(P_X1_XZ(all(P_X1XZ(r,1:end-1)==P_X1_XZ(:,1:end-1),2),end));
    else
        continue
    end
end

H_Z1_XZ = 0;
for r = 1:size(P_Z1XZ,1)
    if ~isnan(P_Z1XZ(r,end)*log2(P_Z1_XZ(all(P_Z1XZ(r,1:end-1)==P_Z1_XZ(:,1:end-1),2),end)))
        H_Z1_XZ = H_Z1_XZ - P_Z1XZ(r,end)*log2(P_Z1_XZ(all(P_Z1XZ(r,1:end-1)==P_Z1_XZ(:,1:end-1),2),end));
    else
        continue
    end
end

H_X1_XY = 0;
for r = 1:size(P_X1XY,1)
    if ~isnan(P_X1XY(r,end)*log2(P_X1_XY(all(P_X1XY(r,1:end-1)==P_X1_XY(:,1:end-1),2),end)))
        H_X1_XY = H_X1_XY - P_X1XY(r,end)*log2(P_X1_XY(all(P_X1XY(r,1:end-1)==P_X1_XY(:,1:end-1),2),end));
    else
        continue
    end
end

H_Y1_XY = 0;
for r = 1:size(P_Y1XY,1)
    if ~isnan(P_Y1XY(r,end)*log2(P_Y1_XY(all(P_Y1XY(r,1:end-1)==P_Y1_XY(:,1:end-1),2),end)))
        H_Y1_XY = H_Y1_XY - P_Y1XY(r,end)*log2(P_Y1_XY(all(P_Y1XY(r,1:end-1)==P_Y1_XY(:,1:end-1),2),end));
    else
        continue
    end
end

H_X1_XYZ = 0;
for r = 1:size(P_X1XYZ,1)
    if ~isnan(P_X1XYZ(r,end)*log2(P_X1_XYZ(all(P_X1XYZ(r,1:end-1)==P_X1_XYZ(:,1:end-1),2),end)))
        H_X1_XYZ = H_X1_XYZ - P_X1XYZ(r,end)*log2(P_X1_XYZ(all(P_X1XYZ(r,1:end-1)==P_X1_XYZ(:,1:end-1),2),end));
    else
        continue
    end
end

H_Y1_XYZ = 0;
for r = 1:size(P_Y1XYZ,1)
    if ~isnan(P_Y1XYZ(r,end)*log2(P_Y1_XYZ(all(P_Y1XYZ(r,1:end-1)==P_Y1_XYZ(:,1:end-1),2),end)))
        H_Y1_XYZ = H_Y1_XYZ - P_Y1XYZ(r,end)*log2(P_Y1_XYZ(all(P_Y1XYZ(r,1:end-1)==P_Y1_XYZ(:,1:end-1),2),end));
    else
        continue
    end
end

H_Z1_XYZ = 0;
for r = 1:size(P_Z1XYZ,1)
    if ~isnan(P_Z1XYZ(r,end)*log2(P_Z1_XYZ(all(P_Z1XYZ(r,1:end-1)==P_Z1_XYZ(:,1:end-1),2),end)))
        H_Z1_XYZ = H_Z1_XYZ - P_Z1XYZ(r,end)*log2(P_Z1_XYZ(all(P_Z1XYZ(r,1:end-1)==P_Z1_XYZ(:,1:end-1),2),end));
    else
        continue
    end
end

% compute transfer entropy
TE_X_YZ = H_Y1_YZ-H_Y1_XYZ; % should be 0.0184 / 0.0146
TE_X_ZY = H_Z1_YZ-H_Z1_XYZ; % should be 0.0029 / 0.0267
TE_Y_XZ = H_X1_XZ-H_X1_XYZ; % should be 0.0134 / 0.0111
TE_Y_ZX = H_Z1_XZ-H_Z1_XYZ; % should be 0.0086 / 0.0345 
TE_Z_XY = H_X1_XY-H_X1_XYZ; % should be 0.0389 / 0.0018
TE_Z_YX = H_Y1_XY-H_Y1_XYZ; % should be 0.0114 / 0.0082


end