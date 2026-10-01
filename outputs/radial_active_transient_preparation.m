%% Physical-time, fixed-Pstar active-transient preparation; no metric search.
root=fileparts(mfilename('fullpath'));out=fullfile(root,['radial_active_transient_' char(datetime('now','Format','yyyyMMdd_HHmmss'))]);mkdir(out);diary(fullfile(out,'console.txt'));
a=load(fullfile(root,'current_limiter_SCR5_XR10','current_limiter_SCR5_XR10_results.mat'),'results');p=a.results.parameters;n0=a.results.network;p.omega=2*pi*50;p.Pstar=1;
xadd=fzero(@(z)imag(rebuild(z,p,n0).Zeq)-.350,[0 .6]);n=rebuild(xadd,p,n0);
gnom=fzero(@(g)coupled(deg2rad(g),1,p,n,0).P-1,[15 25]);nom=coupled(deg2rad(gnom),1,p,n,0);xd=[deg2rad(gnom);1];
% Primary two-state reduction: Q-V filter replaced by instantaneous DC relation.
% PLL frequency fixed to nominal. Electrical impedances evaluated at 50 Hz.
% This is NOT a time rescaling of the earlier normalized-angle model.
retained=[.8 .7 .6];duration=.1;selected=false;
for sag=retained
 pf=p;pf.U=p.U*sag;nf=rebuild(xadd,pf,n0);
 trial=experiment(xd,duration,p,n,pf,nf);
 if any(trial.dFault.active),selected=true;break;end
end
assert(selected,'No activation obtained in defined sag trials');
fprintf('Nominal: gamma=%.12f deg, Eref=%.12f, I=%.12f, u=%.12f, Vc=%.12f, PCC=%.12f\n',gnom,nom.E,abs(nom.ic),abs(nom.u),abs(nom.vc),abs(nom.vg));
fprintf('Selected retained voltage=%.3f (fault U=%.12f); duration=%.3f s.\n',sag,pf.U,duration);
% Locate fixed-Pstar roots in each electrical mode, not an equilibrium family.
rootsByMode=cell(2,1);
for mode=1:2
 pp=p;nn=n;if mode==2,pp=pf;nn=nf;end
 gg=linspace(0,180,36001)';dd=coupled(deg2rad(gg),1,pp,nn,0);rr=find_roots(gg,dd.P-1,@(g)coupled(deg2rad(g),1,pp,nn,0).P-1);vals=cell(numel(rr),1);
 for j=1:numel(rr),z=coupled(deg2rad(rr(j)),1,pp,nn,0);J=[0 p.omega;-z.totalPg/p.mass -p.D/p.mass];vals{j}=struct('gamma',rr(j),'active',z.active,'u',abs(z.u),'eig',eig(J),'Pg',z.totalPg);end
 rootsByMode{mode}=vals;fprintf('Mode %d fixed-Pstar roots (deg): %s\n',mode,mat2str(rr',12));
end
runs=cell(27,1);rows=cell(27,12);k=0;cloud=[];faultCloud=[];postCloud=[];
for dg=[-.1 0 .1]
 for dv=[-1e-4 0 1e-4]
  for tf=[.08 .10 .12]
   k=k+1;r=experiment(xd+[deg2rad(dg);dv],tf,p,n,pf,nf);runs{k}=r;
   faultCloud=[faultCloud;r.xFault(r.dFault.active,:)];postCloud=[postCloud;r.xPost(r.dPost.active,:)];
   if dg==0&&dv==0&&abs(tf-.1)<1e-12,primary=r;end
   rows(k,:)={dg,dv,tf,mat2str(r.currentEventsFault',8),mat2str(r.currentEventsPost',8),rad2deg(r.clearState(1)),r.clearState(2),rad2deg(r.finalState(1)),r.finalState(2),r.recovery,max([abs(r.dFault.u);abs(r.dPost.u)]),r.valid};
  end
 end
end
T=cell2table(rows,'VariableNames',{'InitialAngleOffset_deg','InitialNuOffset','ClearingTime_s','FaultCurrentEvents_s','PostCurrentEvents_s','ClearingGamma_deg','ClearingNu','FinalGamma_deg','FinalNu','Recovery','MaxRequiredCommand','Valid'});disp(T);
% Separate candidate boxes for different vector fields: no common-domain claim.
boxes=cell(2,1);
for mode=1:2
 pp=pf;nn=nf;cc=faultCloud;if mode==2,pp=p;nn=n;cc=postCloud;end
 if isempty(cc),fprintf('No active cloud in mode %d\n',mode);continue;end
 bounds=[min(cc);max(cc)];pad=[deg2rad(.05) 2e-5];bounds(1,:)=bounds(1,:)-pad;bounds(2,:)=bounds(2,:)+pad;
 [G,V]=ndgrid(linspace(bounds(1,1),bounds(2,1),501),linspace(bounds(1,2),bounds(2,2),101));z=coupled(G,1,pp,nn,0);
 A=-z.totalPg./(p.mass*V);B=(-(1-z.P)./V.^2-p.D)/p.mass;
 lambdaMaxReal=real((B+sqrt(complex(B.^2+4*p.omega*A)))/2);
 box=struct('bounds',bounds,'allActive',all(z.active,'all'),'rawMargin',min(z.rawI-p.Imax,[],'all'),'voltageMargin',p.Vmax-max(abs(z.u),[],'all'),'currentMax',max(abs(z.ic),[],'all'),'PgRange',[min(z.totalPg,[],'all') max(z.totalPg,[],'all')],'J21range',[min(A,[],'all') max(A,[],'all')],'J22range',[min(B,[],'all') max(B,[],'all')],'maxRealEig',[min(lambdaMaxReal,[],'all') max(lambdaMaxReal,[],'all')]);
 boxes{mode}=box;fprintf('MODE %d raw box gamma=[%.9f %.9f] deg, nu=[%.12f %.12f], active=%d, rawMargin=%.9g, voltageMargin=%.9g, Pg=[%.9g %.9g]\n',mode,rad2deg(bounds(:,1)),bounds(:,2),box.allActive,box.rawMargin,box.voltageMargin,box.PgRange);
end
% Export actual time series, with both sides of clearing preserved by Stage.
writetable(T,fullfile(out,'ensemble_summary.csv'));
for k=1:2
 tt=primary.tFault;xx=primary.xFault;dd=primary.dFault;label="fault";if k==2,tt=primary.tPost;xx=primary.xPost;dd=primary.dPost;label="postfault";end
 S=table(tt,rad2deg(xx(:,1)),xx(:,2),dd.P,dd.E,dd.rawI,abs(dd.ic),abs(dd.u),dd.active,'VariableNames',{'t_s','gamma_deg','nu','Pe','Eref','RawCurrent','LimitedCurrent','RequiredCommand','Active'});writetable(S,fullfile(out,label+"_trajectory.csv"));
end
save(fullfile(out,'results.mat'),'p','n','pf','nf','nom','gnom','sag','duration','primary','runs','T','boxes','rootsByMode','faultCloud','postCloud','xadd');diary off;
function r=experiment(x0,tf,p,n,pf,nf)
 opts=odeset('RelTol',1e-9,'AbsTol',[1e-11 1e-12],'MaxStep',.002,'Events',@(t,x)events(t,x,pf,nf));
 [t1,x1,te1,~,ie1]=ode45(@(t,x)rhs(x,pf,nf),[0 tf],x0,opts);
 valid=abs(t1(end)-tf)<1e-10;
 if valid
  opts=odeset(opts,'Events',@(t,x)events(t,x,p,n),'MaxStep',.005);
  z=coupled(x1(end,1),1,p,n,0);
  if abs(z.u)>p.Vmax,valid=false;t2=tf;x2=x1(end,:);te2=[];ie2=[];
  else,[t2,x2,te2,~,ie2]=ode45(@(t,x)rhs(x,p,n),[tf 8],x1(end,:)',opts);valid=t2(end)>=8-1e-10;end
 else,t2=t1(end);x2=x1(end,:);te2=[];ie2=[];end
 d1=coupled(x1(:,1),1,pf,nf,0);d2=coupled(x2(:,1),1,p,n,0);
 % Tail recovery to original nominal, NOT to an arbitrary 360-degree branch.
 gn=fzero(@(g)coupled(g,1,p,n,0).P-1,deg2rad([15 25]));tail=t2>=7;
 rec="not yet recovered";
 if ~valid,rec="voltage/model validity stop";
 elseif any(tail)&&max(abs(rad2deg(x2(tail,1)-gn)))<.01&&max(abs(x2(tail,2)-1))<1e-5&&abs(d2.P(end)-1)<1e-4,rec="recovered original equilibrium";
 elseif abs(x2(end,2)-1)<1e-5&&abs(d2.P(end)-1)<1e-4,rec="other branch or unresolved tail";end
 r=struct('tFault',t1,'xFault',x1,'dFault',d1,'tPost',t2,'xPost',x2,'dPost',d2,'currentEventsFault',te1(ie1==1),'currentEventsPost',te2(ie2==1),'clearState',x1(end,:),'finalState',x2(end,:),'recovery',rec,'valid',valid);
end
function f=rhs(x,p,n)
 z=coupled(x(1),1,p,n,0);f=[p.omega*(x(2)-1);((1-z.P)/x(2)-p.D*(x(2)-1))/p.mass];
end
function [v,term,dir]=events(~,x,p,n)
 z=coupled(x(1),1,p,n,0);v=[z.rawI-p.Imax;p.Vmax-abs(z.u);x(2)-.5];term=[0;1;1];dir=[0;-1;-1];
end
function d=coupled(g,vs,p,n,mode)
 E=vs+zeros(size(g));
 for it=1:60
  d=raw(g,E,p,n,mode);h=E+.1*d.Q-vs;he=1+.1*d.QE;
  if max(abs(h),[],'all')<2e-14,break;end
  E=E-max(-.1,min(.1,h./he));
 end
 assert(all(E>0,'all')&&max(abs(h),[],'all')<1e-11&&all(abs(he)>1e-6,'all'));
 d.E=E;d.hE=he;d.Eg=-.1*d.Qg./he;d.totalPg=d.Pg+d.PE.*d.Eg;
end
function d=raw(g,E,p,n,mode)
 rot=exp(1i*(g+angle(n.Vth)));v=E.*rot;q=(v-n.Vth)/n.Zeq;qg=1i*v/n.Zeq;qE=rot/n.Zeq;r=abs(q);
 active=r>p.Imax;if mode==1,active=false(size(q));elseif mode==2,active=true(size(q));end
 ic=q;icg=qg;icE=qE;
 if any(active(:))
  ii=active;qa=q(ii);ra=r(ii);
  ic(ii)=p.Imax*qa./ra;
  icg(ii)=p.Imax*(qg(ii)./ra-qa.*real(conj(qa).*qg(ii))./ra.^3);
  icE(ii)=p.Imax*(qE(ii)./ra-qa.*real(conj(qa).*qE(ii))./ra.^3);
 end
 vc=n.Vth+n.Zeq*ic;vg=n.vgE*vc+n.vgU*p.U;ig=n.igE*vc+n.igU*p.U;
 vcg=n.Zeq*icg;vcE=n.Zeq*icE;
 S=vg.*conj(ig);Sg=n.vgE*vcg.*conj(ig)+vg.*conj(n.igE*vcg);SE=n.vgE*vcE.*conj(ig)+vg.*conj(n.igE*vcE);
 d=struct('P',real(S),'Q',imag(S),'Pg',real(Sg),'Qg',imag(Sg),'PE',real(SE),'QE',imag(SE),'vc',vc,'vg',vg,'ig',ig,'ic',ic,'u',vc+(n.Zconverter+p.AD)*ic,'rawI',r,'active',active);
end
function n=rebuild(x,p,n0)
n=n0;z1=n0.Zt1+1i*x;z2=n0.Zt2;zg=n0.Zgrid;zb=n0.Zbreaker;
y=[1/z1+1/z2+n0.Ym,-1/z2;-1/z2,1/z2+1/(zg+zb)];co=y\diag([1/z1,1/(zg+zb)]);
iv=co(2,1)/(zg+zb);iu=(co(2,2)-1)/(zg+zb);av=co(2,1)-zb*iv;au=co(2,2)-zb*iu;
icv=(1-co(1,1))/z1+n0.Yc;icu=-co(1,2)/z1;
n.Zt1=z1;n.nodalY=y;n.nodeCoefficients=co;n.vgE=av;n.vgU=au;n.igE=iv;n.igU=iu;n.icE=icv;n.icU=icu;n.Zeq=1/icv;n.Vth=-icu/icv*p.U;
end

function r=find_roots(g,y,f)
ix=find(y(1:end-1).*y(2:end)<0);r=g(abs(y)<1e-12);for j=ix',r(end+1,1)=fzero(f,g(j+[0 1]),optimset('TolX',1e-12));end;r=uniquetol(r,1e-9,'DataScale',1);
end




