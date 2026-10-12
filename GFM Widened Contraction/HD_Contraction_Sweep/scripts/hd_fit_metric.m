function [C,tm]=hd_fit_metric(grid,lam,ep,up,cen,sc,p,n,mode)
% Verbatim copy of the local function 'fit' from
% GFM Contraction Analsys/limiter region/widen_nu_contraction.m (itself copied from
% radial_inactive_metric_rate_trial.m). Only function names were changed.
% H and D enter exclusively through p.mass (=2H) and p.D.
 setlmis([]);ids=zeros(1,6);ids(1)=lmivar(3,[1 2;2 -1]);for j=2:6,ids(j)=lmivar(1,[2 1]);end
 base=[0 0;0 1];
 for k=1:size(grid,1)
  z=grid(k,:);x=cen+sc.*z;d=hd_coupled(x(1),1,p,n,mode);v=x(2);f=[p.omega*(v-1),((1-d.P)/v-p.D*(v-1))/p.mass];J=[0 p.omega;-d.totalPg/(p.mass*v),(-(1-d.P)/v^2-p.D)/p.mass];
  b=[1 z z(1)^2 z(1)*z(2) z(2)^2];dz=f./sc;bd=[0 dz 2*z(1)*dz(1) z(2)*dz(1)+z(1)*dz(2) 2*z(2)*dz(2)];
  il=newlmi;iu=newlmi;iq=newlmi;
  lmiterm([il 1 1 0],ep*eye(2)-base);lmiterm([iu 1 1 0],base-up*eye(2));lmiterm([iq 1 1 0],J'*base+base*J+2*lam*base+ep*eye(2));
  for j=1:6
   lmiterm([il 1 1 ids(j)],-b(j),1);lmiterm([iu 1 1 ids(j)],b(j),1);
   lmiterm([iq 1 1 ids(j)],b(j),J,'s');lmiterm([iq 1 1 ids(j)],bd(j)+2*lam*b(j),1);
  end
 end
 sys=getlmis;[tm,xx]=feasp(sys,[1e-8 200 1e4 0 1]);C=zeros(2,2,6);
 if ~isempty(xx),for j=1:6,C(:,:,j)=dec2mat(sys,xx,ids(j));end;C(:,:,1)=C(:,:,1)+base;end
end
