function nearest = isnearest(uc, nuc)
%%% jusify whether uc is the nearest unit cell to nuc
% (uc)by(nuc)
[xa,xb] = meshgrid(nuc(:,1),uc(:,1));
%xa ->nuc (#uc)by(#nuc)
%xb ->uc  (#uc)by(#nuc)
delx = xa-xb; %(#uc)by(#nuc), x_nuc - x_uc
[ya,yb] = meshgrid(nuc(:,2),uc(:,2));
dely = ya - yb;
tol = 10^-4;
% nearest.id = abs(abs(delx)-0.5)<tol&abs(abs(dely))- 1/2/sqrt(3)<tol;
nearest.id = abs(sqrt(delx.^2+dely.^2) - 1/sqrt(3))<tol;%0 or 1
nearest.xa = xa;
nearest.xb = xb;
nearest.ya = ya;
nearest.yb = yb;
end