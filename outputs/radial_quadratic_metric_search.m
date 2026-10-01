%% One common quadratic metric: NOMINAL voltage, Pstar=1, physical seconds.
root=fileparts(mfilename('fullpath'));out=fullfile(root,['radial_quadratic_metric_' char(datetime('now','Format','yyyyMMdd_HHmmss'))]);mkdir(out);diary(fullfile(out,'console.txt'));
a=load(fullfile(root,'radial_active_transient_20261001_152918','results.mat'),'p','n','primary');p=a.p;n=a.n;
bounds=[deg2rad(24.423430419712) .996174514620;deg2rad(33.858410380369) 1.005178107392];center=mean(bounds);scale=diff(bounds)/2;
sw=fzero(@(g)coupled(deg2rad(g),1,p,n,0).rawI-1.2,[24 25]);saddle=fzero(@(g)coupled(deg2rad(g),1,p,n,0).P-1,[60 70]);assert(sw<rad2deg(bounds(1,1))&&saddle>rad2deg(bounds(2,1)));
fprintf('Rechecked switch=%.12f deg, saddle=%.12f deg\n',sw,saddle);
[gg,vv]=ndgrid(linspace(-1,1,9));grid=[gg(:),vv(:)];epsilon=1e-6;upper=10;
% Fixed-lambda search is LMI in all 17 independent coefficients; trace M0=1.
lo=0;hi=.01;history=[];best=[];
for j=1:12
 [C,tm]=fit(grid,hi,epsilon,upper,center,scale,p,n);history=[history;hi tm];fprintf('lambda %.8g tmin %.9g\n',hi,tm);
 if tm>=-1e-9,break;end
 best=C;lo=hi;hi=2*hi;
end
if isempty(best),error('Initial full-domain degree-2 LMI failed: inspect smaller subsets before conclusion');end
for j=1:14
 mid=(lo+hi)/2;[C,tm]=fit(grid,mid,epsilon,upper,center,scale,p,n);history=[history;mid tm];
 if tm < -1e-9,lo=mid;best=C;else,hi=mid;end
end
lambda=.8*lo;fprintf('Grid rate bracket [%.12g %.12g], selected safety rate %.12g\n',lo,hi,lambda);
for it=1:6
 [C,tm]=fit(grid,lambda,epsilon,upper,center,scale,p,n);
 assert(tm<0);[zg,zv]=ndgrid(linspace(-1,1,101));test=[zg(:),zv(:)];v=verify(C,test,lambda,center,scale,p,n);
 fprintf('Refinement %d: minM %.12g maxM %.12g worstQ %.12g\n',it,min(v.minM),max(v.maxM),max(v.maxQ));
 bad=find(v.minM<epsilon|v.maxM>upper|v.maxQ>-epsilon);
 if isempty(bad),break;end
 scores=max([epsilon-v.minM,v.maxM-upper,v.maxQ+epsilon],[],2);[~,ix]=sort(scores,'descend');grid=unique([grid;test(ix(1:min(30,numel(bad))),:)],'rows');
end
assert(isempty(bad),'Adaptive refinement did not obtain feasible grid');
[zg,zv]=ndgrid(linspace(-1,1,301));dense=verify(C,[zg(:),zv(:)],lambda,center,scale,p,n);
assert(min(dense.minM)>epsilon&&max(dense.maxQ)<-epsilon&&max(dense.maxM)<upper);
fprintf('Dense minM %.15g, maxM %.15g, worstQ %.15g, rate %.15g\n',min(dense.minM),max(dense.maxM),max(dense.maxQ),lambda);
format long g;for j=1:6,fprintf('M%d\n',j-1);disp(C(:,:,j));end
fprintf('Center=[%.17g %.17g], scale=[%.17g %.17g]\n',center,scale);
fprintf('Instantaneous max-real-eigenvalue range [%.12g %.12g]; raw branch margin %.12g; voltage margin %.12g\n',min(dense.eigJ),max(dense.eigJ),min(dense.rawI)-1.2,1.3-max(dense.u));
% Check Mdot by directional differentiation at off-equilibrium points.
err=0;for z1=linspace(-.9,.9,7),for z2=linspace(-.9,.9,7)
 z=[z1 z2];x=center+scale.*z;d=coupled(x(1),1,p,n,0);f=[p.omega*(x(2)-1),((1-d.P)/x(2)-25*(x(2)-1))/6.8];
 [M,Mdot]=metric(C,z,f./scale);h=1e-7;Mp=metric(C,z+h*f./scale,[0 0]);Mm=metric(C,z-h*f./scale,[0 0]);err=max(err,max(abs((Mp-Mm)/(2*h)-Mdot),[],'all'));
end,end;assert(err<1e-6);fprintf('Mdot directional derivative error %.12g\n',err);
pr=a.primary;active=pr.dPost.active;inside=pr.xPost(:,1)>=bounds(1,1)&pr.xPost(:,1)<=bounds(2,1)&pr.xPost(:,2)>=bounds(1,2)&pr.xPost(:,2)<=bounds(2,2);fprintf('Primary sampled active postfault points inside domain: %d/%d; excluded only near release? %d\n',nnz(active&inside),nnz(active),all(pr.xPost(active&~inside,1)<bounds(1,1)));
containedInterval=[pr.tPost(find(active&inside,1)),pr.tPost(find(active&inside,1,'last'))];
save(fullfile(out,'metric.mat'),'C','lambda','bounds','center','scale','epsilon','upper','history','dense','grid','sw','saddle','containedInterval','err');
% JSON for independent interval checking; coefficients in PHYSICAL x metric.
names={'Zeq','Vth','vgE','vgU','igE','igU','Zconverter'};net=struct();for j=1:numel(names),z=n.(names{j});net.(names{j})=[real(z),imag(z)];end
coef=cell(6,1);for j=1:6,coef{j}=C(:,:,j);end
export=struct('coefficients',{coef},'lambda',lambda,'bounds',bounds,'center',center,'scale',scale,'epsilon',epsilon,'upper',upper,'network',net,'U',p.U,'AD',p.AD,'omega',p.omega,'mass',p.mass,'D',p.D,'Imax',p.Imax,'umax',p.Vmax);
fid=fopen(fullfile(out,'metric.json'),'w');fwrite(fid,jsonencode(export),'char');fclose(fid);diary off;
function [C,tm]=fit(grid,lam,ep,up,cen,sc,p,n)
 setlmis([]);ids=zeros(1,6);ids(1)=lmivar(3,[1 2;2 -1]);for j=2:6,ids(j)=lmivar(1,[2 1]);end
 base=[0 0;0 1];
 for k=1:size(grid,1)
  z=grid(k,:);x=cen+sc.*z;d=coupled(x(1),1,p,n,0);v=x(2);f=[p.omega*(v-1),((1-d.P)/v-p.D*(v-1))/p.mass];J=[0 p.omega;-d.totalPg/(p.mass*v),(-(1-d.P)/v^2-p.D)/p.mass];
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
function v=verify(C,grid,lam,cen,sc,p,n)
 x=cen+grid.*sc;d=coupled(x(:,1),1,p,n,0);v=struct('minM',zeros(size(grid,1),1),'maxM',zeros(size(grid,1),1),'maxQ',zeros(size(grid,1),1),'eigJ',zeros(size(grid,1),1),'rawI',d.rawI,'u',abs(d.u));
 for k=1:size(grid,1)
  nu=x(k,2);f=[p.omega*(nu-1),((1-d.P(k))/nu-p.D*(nu-1))/p.mass];J=[0 p.omega;-d.totalPg(k)/(p.mass*nu),(-(1-d.P(k))/nu^2-p.D)/p.mass];
  [M,Mdot]=metric(C,grid(k,:),f./sc);Q=J'*M+M*J+Mdot+2*lam*M;
  v.minM(k)=min(eig(M));v.maxM(k)=max(eig(M));v.maxQ(k)=max(eig(Q));v.eigJ(k)=max(real(eig(J)));
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




