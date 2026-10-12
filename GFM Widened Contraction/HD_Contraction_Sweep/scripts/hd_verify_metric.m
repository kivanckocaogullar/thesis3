function v=hd_verify_metric(C,grid,lam,cen,sc,p,n,mode)
% Verbatim copy of the local function 'verify' from
% GFM Contraction Analsys/limiter region/widen_nu_contraction.m (itself copied from
% radial_inactive_metric_rate_trial.m). Only function names were changed.
% H and D enter exclusively through p.mass (=2H) and p.D.
 x=cen+grid.*sc;d=hd_coupled(x(:,1),1,p,n,mode);v=struct('minM',zeros(size(grid,1),1),'maxM',zeros(size(grid,1),1),'maxQ',zeros(size(grid,1),1),'active',d.active);
 for k=1:size(grid,1)
  nu=x(k,2);f=[p.omega*(nu-1),((1-d.P(k))/nu-p.D*(nu-1))/p.mass];J=[0 p.omega;-d.totalPg(k)/(p.mass*nu),(-(1-d.P(k))/nu^2-p.D)/p.mass];
  [M,Mdot]=hd_metric_poly(C,grid(k,:),f./sc);Q=J'*M+M*J+Mdot+2*lam*M;
  v.minM(k)=min(eig(M));v.maxM(k)=max(eig(M));v.maxQ(k)=max(eig(Q));
 end
end
