function R=widen_nu_contraction(cfg)
% Refit the SAME quadratic metric LMI (same lambda, eps, M<upper*I) on rectangles
% with a wider nu range. The gamma range is kept from the certified rectangles
% unless cfg.gamma is given. Bisection between the old nu range and cfg.target.
% Result is a sampled (301x301) check only; run export_widened_metric_json and
% widened_interval_checker.py for the continuous-domain certificate.
% The certified result is in metrics/ (target {[.98 1.02],[.98 1.02]}, gamma {[10 24.3234],[24.4234 43]}).
if nargin<1,cfg=struct();end
root=fileparts(mfilename('fullpath'));a=load(fullfile(root,'inputs','model_inputs.mat'),'p','n');p=a.p;n=a.n;
src={fullfile('reference_old_metrics','inactive'),fullfile('reference_old_metrics','active')};mode=[1 0];labels={'inactive','active'};
tgt={[.990 1.010],[.990 1.015]};if isfield(cfg,'target'),tgt=cfg.target;end
nbis=8;if isfield(cfg,'nbis'),nbis=cfg.nbis;end
out=fullfile(root,'search_runs',['widened_nu_' char(datetime('now','Format','yyyyMMdd_HHmmss'))]);mkdir(out);diary(fullfile(out,'console.txt'));
R=struct();
for j=1:2
 m=load(fullfile(root,src{j},'metric.mat'),'bounds','lambda','epsilon','upper');b0=m.bounds;if isfield(cfg,'lambda'),m.lambda=cfg.lambda(j);end
 if isfield(cfg,'gamma'),b0(:,1)=deg2rad(cfg.gamma{j}(:));end
 bt=b0;bt(:,2)=tgt{j}(:);
 fprintf('\n=== %s | lambda=%g | gamma [%.4f %.4f] deg | old nu [%.6f %.6f] -> target [%.6f %.6f]\n',labels{j},m.lambda,rad2deg(b0(:,1)),b0(:,2),bt(:,2));
 [ok,C,dense]=attempt(bt,m,mode(j),p,n);best=[];
 if ok,best=struct('t',1,'bounds',bt,'C',C,'dense',dense);
 else
  lo=0;hi=1;
  for it=1:nbis
   t=(lo+hi)/2;bb=b0+t*(bt-b0);[ok,C,dense]=attempt(bb,m,mode(j),p,n);
   fprintf('  bisect t=%.4f nu [%.6f %.6f] -> %d\n',t,bb(:,2),ok);
   if ok,lo=t;best=struct('t',t,'bounds',bb,'C',C,'dense',dense);else,hi=t;end
  end
 end
 if isempty(best),fprintf('  no widening found\n');R.(labels{j})=[];continue;end
 fprintf('  RESULT nu [%.6f %.6f] (t=%.4f) minM %.3g maxM %.3g worstQ %.3g\n',best.bounds(:,2),best.t,min(best.dense.minM),max(best.dense.maxM),max(best.dense.maxQ));
 C=best.C;bounds=best.bounds;center=mean(bounds,1);scale=diff(bounds,1)/2;lambda=m.lambda;epsilon=m.epsilon;upper=m.upper;dense=best.dense;
 save(fullfile(out,[labels{j} '_metric.mat']),'C','lambda','epsilon','upper','bounds','center','scale','dense');
 R.(labels{j})=best;
end
diary off;
end

function [ok,C,dense]=attempt(bounds,m,mode,p,n)
 center=mean(bounds,1);scale=diff(bounds,1)/2;lam=m.lambda;ep=m.epsilon;up=m.upper;
 [g,v]=ndgrid(linspace(-1,1,9));grid=[g(:),v(:)];ok=false;C=[];dense=[];
 for it=1:8
  [C,tm]=fit(grid,lam,ep,up,center,scale,p,n,mode);if ~(tm<0),return;end
  [g,v]=ndgrid(linspace(-1,1,101));test=[g(:),v(:)];z=verify(C,test,lam,center,scale,p,n,mode);
  scores=max([ep-z.minM,z.maxM-up,z.maxQ+ep],[],2);bad=find(scores>0);if isempty(bad),break;end
  [~,ix]=sort(scores,'descend');grid=unique([grid;test(ix(1:min(30,numel(bad))),:)],'rows');
 end
 if ~isempty(bad),return;end
 [g,v]=ndgrid(linspace(-1,1,301));dense=verify(C,[g(:),v(:)],lam,center,scale,p,n,mode);
 ok=min(dense.minM)>ep&&max(dense.maxM)<up&&max(dense.maxQ)<-ep;
end

function [C,tm]=fit(grid,lam,ep,up,cen,sc,p,n,mode)
 setlmis([]);ids=zeros(1,6);ids(1)=lmivar(3,[1 2;2 -1]);for j=2:6,ids(j)=lmivar(1,[2 1]);end
 base=[0 0;0 1];
 for k=1:size(grid,1)
  z=grid(k,:);x=cen+sc.*z;d=coupled(x(1),1,p,n,mode);v=x(2);f=[p.omega*(v-1),((1-d.P)/v-p.D*(v-1))/p.mass];J=[0 p.omega;-d.totalPg/(p.mass*v),(-(1-d.P)/v^2-p.D)/p.mass];
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
function v=verify(C,grid,lam,cen,sc,p,n,mode)
 x=cen+grid.*sc;d=coupled(x(:,1),1,p,n,mode);v=struct('minM',zeros(size(grid,1),1),'maxM',zeros(size(grid,1),1),'maxQ',zeros(size(grid,1),1),'active',d.active);
 for k=1:size(grid,1)
  nu=x(k,2);f=[p.omega*(nu-1),((1-d.P(k))/nu-p.D*(nu-1))/p.mass];J=[0 p.omega;-d.totalPg(k)/(p.mass*nu),(-(1-d.P(k))/nu^2-p.D)/p.mass];
  [M,Mdot]=metric(C,grid(k,:),f./sc);Q=J'*M+M*J+Mdot+2*lam*M;
  v.minM(k)=min(eig(M));v.maxM(k)=max(eig(M));v.maxQ(k)=max(eig(Q));
 end
end
function [M,Mdot]=metric(C,z,dz)
 b=[1 z z(1)^2 z(1)*z(2) z(2)^2];bd=[0 dz 2*z(1)*dz(1) z(2)*dz(1)+z(1)*dz(2) 2*z(2)*dz(2)];M=zeros(2);Mdot=M;for j=1:6,M=M+b(j)*C(:,:,j);Mdot=Mdot+bd(j)*C(:,:,j);end
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
