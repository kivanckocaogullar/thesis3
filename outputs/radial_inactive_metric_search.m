root=fileparts(mfilename('fullpath'));a=load(fullfile(root,'radial_active_transient_20261001_152918','results.mat'),'p','n','primary');p=a.p;n=a.n;
out=fullfile(root,['radial_inactive_metric_' char(datetime('now','Format','yyyyMMdd_HHmmss'))]);mkdir(out);diary(fullfile(out,'console.txt'));
gsw=24.373430419712*pi/180; xx=a.primary.xPost; inactive=~a.primary.dPost.active; xr=xx(inactive,:);
bounds=[min(xr(:,1))-.05*pi/180,min(xr(:,2))-2e-5;gsw-.05*pi/180,max(xr(:,2))+2e-5];center=mean(bounds,1);scale=diff(bounds,1)/2;
[g,v]=ndgrid(linspace(-1,1,9));grid=[g(:),v(:)];
inside=all(xr>=bounds(1,:) & xr<=bounds(2,:),2);fprintf('Bounds deg %.12f %.12f nu %.12f %.12f; inactive samples %d/%d\n',bounds(:,1)*180/pi,bounds(:,2),sum(inside),numel(inside));
save(fullfile(out,'containment.mat'),'inside','xr','bounds');
lambda=.1;epsilon=1e-5;upper=10;
for it=1:6
 [C,tm]=fit(grid,lambda,epsilon,upper,center,scale,p,n);assert(tm<0);
 [g,v]=ndgrid(linspace(-1,1,101));test=[g(:),v(:)];z=verify(C,test,lambda,center,scale,p,n);
 fprintf('Round %d tmin=%.12g minM=%.12g maxM=%.12g worstQ=%.12g\n',it,tm,min(z.minM),max(z.maxM),max(z.maxQ));
 scores=max([epsilon-z.minM,z.maxM-upper,z.maxQ+epsilon],[],2);bad=find(scores>0);if isempty(bad),break;end
 [~,ix]=sort(scores,'descend');grid=unique([grid;test(ix(1:min(30,numel(bad))),:)],'rows');
end
assert(isempty(bad));[g,v]=ndgrid(linspace(-1,1,301));dense=verify(C,[g(:),v(:)],lambda,center,scale,p,n);assert(min(dense.minM)>epsilon&&max(dense.maxQ)<-epsilon);
format long g;for j=1:6,fprintf('M%d\n',j-1);disp(C(:,:,j));end
fprintf('Dense minM %.15g maxM %.15g worstQ %.15g\n',min(dense.minM),max(dense.maxM),max(dense.maxQ));
coef=cell(6,1);for j=1:6,coef{j}=C(:,:,j);end
names={'Zeq','Vth','vgE','vgU','igE','igU','Zconverter'};net=struct();for j=1:numel(names),zz=n.(names{j});net.(names{j})=[real(zz),imag(zz)];end
export=struct('coefficients',{coef},'lambda',lambda,'bounds',bounds,'center',center,'scale',scale,'epsilon',epsilon,'upper',upper,'network',net,'U',p.U,'AD',p.AD,'omega',p.omega,'mass',p.mass,'D',p.D,'Imax',p.Imax,'umax',p.Vmax);
fid=fopen(fullfile(out,'metric.json'),'w');fwrite(fid,jsonencode(export),'char');fclose(fid);save(fullfile(out,'metric.mat'),'C','lambda','epsilon','upper','bounds','center','scale','dense','grid');diary off;
function [C,tm]=fit(grid,lam,ep,up,cen,sc,p,n)
 setlmis([]);ids=zeros(1,6);ids(1)=lmivar(3,[1 2;2 -1]);for j=2:6,ids(j)=lmivar(1,[2 1]);end
 base=[0 0;0 1];
 for k=1:size(grid,1)
  z=grid(k,:);x=cen+sc.*z;d=coupled(x(1),1,p,n,1);v=x(2);f=[p.omega*(v-1),((1-d.P)/v-p.D*(v-1))/p.mass];J=[0 p.omega;-d.totalPg/(p.mass*v),(-(1-d.P)/v^2-p.D)/p.mass];
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
 x=cen+grid.*sc;d=coupled(x(:,1),1,p,n,1);v=struct('minM',zeros(size(grid,1),1),'maxM',zeros(size(grid,1),1),'maxQ',zeros(size(grid,1),1),'eigJ',zeros(size(grid,1),1),'rawI',d.rawI,'u',abs(d.u));
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






