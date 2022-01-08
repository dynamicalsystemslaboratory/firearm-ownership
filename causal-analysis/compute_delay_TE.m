function [TE_delay] = compute_delay_TE(XYZ,d1,d2)

if d1 > d2

    % X->Y|Z
    XYZX1Y1Z1 = [XYZ.X(1:(end-1-d1)) XYZ.Y((1+d1):(end-1)) XYZ.Z((1+d1-d2):(end-1-d2)) XYZ.X((2+d1):end) XYZ.Y((2+d1):end) XYZ.Z((2+d1):end)];
    XYZX1Y1Z1 = array2table(XYZX1Y1Z1,'VariableNames',{'X','Y','Z','X1','Y1','Z1'});
    [TE_X_YZ,~,~,~,~,~] = conditional_TE(XYZX1Y1Z1);

    % X->Z|Y
    XYZX1Y1Z1 = [XYZ.X(1:(end-1-d1)) XYZ.Y((1+d1-d2):(end-1-d2)) XYZ.Z((1+d1):(end-1)) XYZ.X((2+d1):end) XYZ.Y((2+d1):end) XYZ.Z((2+d1):end)];
    XYZX1Y1Z1 = array2table(XYZX1Y1Z1,'VariableNames',{'X','Y','Z','X1','Y1','Z1'});
    [~,TE_X_ZY,~,~,~,~] = conditional_TE(XYZX1Y1Z1);

    % Y->X|Z
    XYZX1Y1Z1 = [XYZ.X((1+d1):(end-1)) XYZ.Y(1:(end-1-d1)) XYZ.Z((1+d1-d2):(end-1-d2)) XYZ.X((2+d1):end) XYZ.Y((2+d1):end) XYZ.Z((2+d1):end)];
    XYZX1Y1Z1 = array2table(XYZX1Y1Z1,'VariableNames',{'X','Y','Z','X1','Y1','Z1'});
    [~,~,TE_Y_XZ,~,~,~] = conditional_TE(XYZX1Y1Z1);
    
    % Y->Z|X
    XYZX1Y1Z1 = [XYZ.X((1+d1-d2):(end-1-d2)) XYZ.Y(1:(end-1-d1)) XYZ.Z((1+d1):(end-1)) XYZ.X((2+d1):end) XYZ.Y((2+d1):end) XYZ.Z((2+d1):end)];
    XYZX1Y1Z1 = array2table(XYZX1Y1Z1,'VariableNames',{'X','Y','Z','X1','Y1','Z1'});
    [~,~,~,TE_Y_ZX,~,~] = conditional_TE(XYZX1Y1Z1);

    % Z->X|Y
    XYZX1Y1Z1 = [XYZ.X((1+d1):(end-1)) XYZ.Y((1+d1-d2):(end-1-d2)) XYZ.Z(1:(end-1-d1)) XYZ.X((2+d1):end) XYZ.Y((2+d1):end) XYZ.Z((2+d1):end)];
    XYZX1Y1Z1 = array2table(XYZX1Y1Z1,'VariableNames',{'X','Y','Z','X1','Y1','Z1'});
    [~,~,~,~,TE_Z_XY,~] = conditional_TE(XYZX1Y1Z1);
    
    % Z->Y|X
    XYZX1Y1Z1 = [XYZ.X((1+d1-d2):(end-1-d2)) XYZ.Y((1+d1):(end-1)) XYZ.Z(1:(end-1-d1)) XYZ.X((2+d1):end) XYZ.Y((2+d1):end) XYZ.Z((2+d1):end)];
    XYZX1Y1Z1 = array2table(XYZX1Y1Z1,'VariableNames',{'X','Y','Z','X1','Y1','Z1'});
    [~,~,~,~,~,TE_Z_YX] = conditional_TE(XYZX1Y1Z1);


elseif d2 > d1
    
    % X->Y|Z
    XYZX1Y1Z1 = [XYZ.X((1+d2-d1):(end-1-d1)) XYZ.Y((1+d2):(end-1)) XYZ.Z(1:(end-1-d2)) XYZ.X((2+d2):end) XYZ.Y((2+d2):end) XYZ.Z((2+d2):end)];
    XYZX1Y1Z1 = array2table(XYZX1Y1Z1,'VariableNames',{'X','Y','Z','X1','Y1','Z1'});
    [TE_X_YZ,~,~,~,~,~] = conditional_TE(XYZX1Y1Z1);

    % X->Z|Y
    XYZX1Y1Z1 = [XYZ.X((1+d2-d1):(end-1-d1)) XYZ.Y(1:(end-1-d2)) XYZ.Z((1+d2):(end-1)) XYZ.X((2+d2):end) XYZ.Y((2+d2):end) XYZ.Z((2+d2):end)];
    XYZX1Y1Z1 = array2table(XYZX1Y1Z1,'VariableNames',{'X','Y','Z','X1','Y1','Z1'});
    [~,TE_X_ZY,~,~,~,~] = conditional_TE(XYZX1Y1Z1);

    % Y->X|Z
    XYZX1Y1Z1 = [XYZ.X((1+d2):(end-1)) XYZ.Y((1+d2-d1):(end-1-d1)) XYZ.Z(1:(end-1-d2)) XYZ.X((2+d2):end) XYZ.Y((2+d2):end) XYZ.Z((2+d2):end)];
    XYZX1Y1Z1 = array2table(XYZX1Y1Z1,'VariableNames',{'X','Y','Z','X1','Y1','Z1'});
    [~,~,TE_Y_XZ,~,~,~] = conditional_TE(XYZX1Y1Z1);
    
    % Y->Z|X
    XYZX1Y1Z1 = [XYZ.X(1:(end-1-d2)) XYZ.Y((1+d2-d1):(end-1-d1)) XYZ.Z((1+d2):(end-1)) XYZ.X((2+d2):end) XYZ.Y((2+d2):end) XYZ.Z((2+d2):end)];
    XYZX1Y1Z1 = array2table(XYZX1Y1Z1,'VariableNames',{'X','Y','Z','X1','Y1','Z1'});
    [~,~,~,TE_Y_ZX,~,~] = conditional_TE(XYZX1Y1Z1);

    % Z->X|Y
    XYZX1Y1Z1 = [XYZ.X((1+d2):(end-1)) XYZ.Y(1:(end-1-d2)) XYZ.Z((1+d2-d1):(end-1-d1)) XYZ.X((2+d2):end) XYZ.Y((2+d2):end) XYZ.Z((2+d2):end)];
    XYZX1Y1Z1 = array2table(XYZX1Y1Z1,'VariableNames',{'X','Y','Z','X1','Y1','Z1'});
    [~,~,~,~,TE_Z_XY,~] = conditional_TE(XYZX1Y1Z1);
    
    % Z->Y|X
    XYZX1Y1Z1 = [XYZ.X(1:(end-1-d2)) XYZ.Y((1+d2):(end-1)) XYZ.Z((1+d2-d1):(end-1-d1)) XYZ.X((2+d2):end) XYZ.Y((2+d2):end) XYZ.Z((2+d2):end)];
    XYZX1Y1Z1 = array2table(XYZX1Y1Z1,'VariableNames',{'X','Y','Z','X1','Y1','Z1'});
    [~,~,~,~,~,TE_Z_YX] = conditional_TE(XYZX1Y1Z1);
    
else % d1=d2
    
    % X->Y|Z
    XYZX1Y1Z1 = [XYZ.X(1:(end-1-d1)) XYZ.Y((1+d1):(end-1)) XYZ.Z(1:(end-1-d2)) XYZ.X((2+d1):(end)) XYZ.Y((2+d1):(end)) XYZ.Z((2+d1):(end))];
    XYZX1Y1Z1 = array2table(XYZX1Y1Z1,'VariableNames',{'X','Y','Z','X1','Y1','Z1'});
    [TE_X_YZ,~,~,~,~,~] = conditional_TE(XYZX1Y1Z1);

    % X->Z|Y
    XYZX1Y1Z1 = [XYZ.X(1:(end-1-d1)) XYZ.Y(1:(end-1-d2)) XYZ.Z((1+d1):(end-1)) XYZ.X((2+d1):(end)) XYZ.Y((2+d1):(end)) XYZ.Z((2+d1):(end))];
    XYZX1Y1Z1 = array2table(XYZX1Y1Z1,'VariableNames',{'X','Y','Z','X1','Y1','Z1'});
    [~,TE_X_ZY,~,~,~,~] = conditional_TE(XYZX1Y1Z1);

    % Y->X|Z
    XYZX1Y1Z1 = [XYZ.X((1+d1):(end-1)) XYZ.Y(1:(end-1-d1)) XYZ.Z(1:(end-1-d2)) XYZ.X((2+d1):(end)) XYZ.Y((2+d1):(end)) XYZ.Z((2+d1):(end))];
    XYZX1Y1Z1 = array2table(XYZX1Y1Z1,'VariableNames',{'X','Y','Z','X1','Y1','Z1'});
    [~,~,TE_Y_XZ,~,~,~] = conditional_TE(XYZX1Y1Z1);
    
    % Y->Z|X
    XYZX1Y1Z1 = [XYZ.X(1:(end-1-d2)) XYZ.Y(1:(end-1-d1)) XYZ.Z((1+d1):(end-1)) XYZ.X((2+d1):(end)) XYZ.Y((2+d1):(end)) XYZ.Z((2+d1):(end))];
    XYZX1Y1Z1 = array2table(XYZX1Y1Z1,'VariableNames',{'X','Y','Z','X1','Y1','Z1'});
    [~,~,~,TE_Y_ZX,~,~] = conditional_TE(XYZX1Y1Z1);

    % Z->X|Y
    XYZX1Y1Z1 = [XYZ.X((1+d1):(end-1)) XYZ.Y(1:(end-1-d2)) XYZ.Z(1:(end-1-d1)) XYZ.X((2+d1):(end)) XYZ.Y((2+d1):(end)) XYZ.Z((2+d1):(end))];
    XYZX1Y1Z1 = array2table(XYZX1Y1Z1,'VariableNames',{'X','Y','Z','X1','Y1','Z1'});
    [~,~,~,~,TE_Z_XY,~] = conditional_TE(XYZX1Y1Z1);
    
    % Z->Y|X
    XYZX1Y1Z1 = [XYZ.X(1:(end-1-d2)) XYZ.Y((1+d1):(end-1)) XYZ.Z(1:(end-1-d1)) XYZ.X((2+d1):(end)) XYZ.Y((2+d1):(end)) XYZ.Z((2+d1):(end))];
    XYZX1Y1Z1 = array2table(XYZX1Y1Z1,'VariableNames',{'X','Y','Z','X1','Y1','Z1'});
    [~,~,~,~,~,TE_Z_YX] = conditional_TE(XYZX1Y1Z1);
      
end


TE_delay = [TE_X_YZ,TE_X_ZY,TE_Y_XZ,TE_Y_ZX,TE_Z_XY,TE_Z_YX];


end 