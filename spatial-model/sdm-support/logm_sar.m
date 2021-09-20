function logm = logm_sar(y,xo,W,lflag)
rmin = -0.9999;
rmax = 0.9999;
mcorder = 30;
mciter = 50;
nobs = length(y);

x=xo;
xpx = x'*x;
lndetx_sar = log(det(xpx));
dof = (nobs -1)/2;
D = (1 - 1/rmin); % from uniform prior on rho
logC_sar = -log(D) + gammaln(dof) - dof*log(2*pi) -0.5*lndetx_sar;
Wy = sparse(W)*y;
bo = (xpx)\(x'*y);
bd = (xpx)\(x'*Wy);
eo = y - x*bo;
ed = Wy - x*bd;
epeo = eo'*eo;
eped = ed'*ed;
epeod = ed'*eo;
incr = 0.001;
xxp=rmin:incr:rmax;
xx = xxp';
ngrid = length(xx);
iotan = ones(ngrid,1);
logm_sar_profile = -dof*log(epeo*iotan - 2*xx*epeod + (xx.*xx)*eped) ...
+ lndetmc_v1(mcorder,mciter,W,xx);
[adj,mind] = max(logm_sar_profile);
yy = exp(logm_sar_profile -adj);
% trapezoid rule integration
isum = trapz(xx,yy);
isum = isum + adj;
logm = isum + logC_sar; % we put back the scale adjustment here


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




