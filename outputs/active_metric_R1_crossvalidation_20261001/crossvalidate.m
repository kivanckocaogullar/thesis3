% Read-only cross-validation of the full-precision certified active metric.
out=fileparts(mfilename('fullpath')); root=fileparts(out);
source=fullfile(root,'radial_quadratic_robust_20261001_161903','metric.json');
m=jsondecode(fileread(source)); C=permute(m.coefficients,[2 3 1]);
a=load(fullfile(root,'radial_active_transient_20261001_152918','results.mat'),'p','n');p=a.p;n=a.n;
assert(abs(p.U-m.U)<1e-14 && abs(n.Zeq-(m.network.Zeq(1)+1i*m.network.Zeq(2)))<1e-14);
lo=[deg2rad(16.426345885756),0.99652713910517932];hi=[deg2rad(24.373430419712),1.0025];
N=1001; [g,v]=ndgrid(linspace(lo(1),hi(1),N),linspace(lo(2),hi(2),N));
d=coupled(g(:,1),1,p,n,1);
S=screen(g,v,d,C,m); S.domain=[lo;hi];S.gridSize=[N N];S.metricSource=source;
S.rawCurrentRange=[min(d.rawI),max(d.rawI)];S.maxVoltageCommand=max(abs(d.u));
gs=fzero(@(g)switcherror(g,p,n),[.42 .43],optimset('TolX',eps));
vb=linspace(lo(2),hi(2),10001);di=coupled(gs,1,p,n,1);da=coupled(gs,1,p,n,2);
S.switch=struct('gamma_deg',rad2deg(gs),'Pe_inactive',di.P,'Pe_active',da.P,'Pg_inactive',di.totalPg,'Pg_active',da.totalPg,'inactive',screen(gs+zeros(size(vb)),vb,di,C,m),'active',screen(gs+zeros(size(vb)),vb,da,C,m));
S.originalCertificate=jsondecode(fileread(fullfile(fileparts(source),'interval_validation.json')));
S.uniformR1Verification='Not attempted: numerical counterexamples invalidate the requested full-domain certificate.';
S.checkMdot=0;S.checkJacobian=0;
for x=[lo; (lo+hi)/2; hi]'
 x=x';dd=coupled(x(1),1,p,n,1);f=[p.omega*(x(2)-1),((1-dd.P)/x(2)-25*(x(2)-1))/6.8];
 [M,MD]=metricpoint(x,f,C,m);h=1e-7;
 MP=metricpoint(x+h*f,[0 0],C,m);MM=metricpoint(x-h*f,[0 0],C,m);
 S.checkMdot=max(S.checkMdot,max(abs((MP-MM)/(2*h)-MD),[],'all'));
 J=[0 p.omega;-dd.totalPg/(6.8*x(2)),(-(1-dd.P)/x(2)^2-25)/6.8];Jn=zeros(2);
 for j=1:2,dx=[0 0];dx(j)=h;Jn(:,j)=(field(x+dx,p,n)-field(x-dx,p,n))'/(2*h);end
 S.checkJacobian=max(S.checkJacobian,max(abs(J-Jn),[],'all'));
end
assert(S.checkMdot<1e-5 && S.checkJacobian<1e-5);
fid=fopen(fullfile(out,'results.json'),'w');fwrite(fid,jsonencode(S,PrettyPrint=true));fclose(fid);
save(fullfile(out,'results.mat'),'S','m','C','p','n');disp(jsonencode(S,PrettyPrint=true));
function s=screen(g,v,d,C,m)
 z1=(g-m.center(1))/m.scale(1);z2=(v-m.center(2))/m.scale(2);
 f1=m.omega*(v-1);f2=((1-d.P)./v-m.D*(v-1))/m.mass;
 dz1=f1/m.scale(1);dz2=f2/m.scale(2);
 b={ones(size(g)),z1,z2,z1.^2,z1.*z2,z2.^2};bd={zeros(size(g)),dz1,dz2,2*z1.*dz1,z2.*dz1+z1.*dz2,2*z2.*dz2};
 A=zeros(size(g));B=A;D=A;Ad=A;Bd=A;Dd=A;
 for k=1:6,A=A+C(1,1,k)*b{k};B=B+C(1,2,k)*b{k};D=D+C(2,2,k)*b{k};Ad=Ad+C(1,1,k)*bd{k};Bd=Bd+C(1,2,k)*bd{k};Dd=Dd+C(2,2,k)*bd{k};end
 emin=(A+D)/2-hypot((A-D)/2,B);[s.minEigM,ii]=min(emin,[],'all','linear');s.minM_state=[rad2deg(g(ii)),v(ii)];s.nonPositiveM_count=nnz(emin<=0);
 j21=-d.totalPg./(m.mass*v);j22=(-(1-d.P)./v.^2-m.D)/m.mass;
 for lambda=[.1 0]
  qa=2*j21.*B+Ad+2*lambda*A;qb=m.omega*A+j22.*B+j21.*D+Bd+2*lambda*B;qd=2*m.omega*B+2*j22.*D+Dd+2*lambda*D;
  emax=(qa+qd)/2+hypot((qa-qd)/2,qb);[worst,ii]=max(emax,[],'all','linear');
  r=struct('lambda',lambda,'maxEigQc',worst,'worst_state',[rad2deg(g(ii)),v(ii)],'nonNegativeQc_count',nnz(emax>=0),'M_atWorst',[A(ii),B(ii);B(ii),D(ii)],'Qc_atWorst',[qa(ii),qb(ii);qb(ii),qd(ii)]);
  if lambda>0,s.rate01=r;else,s.rate0=r;end
 end
end
function [M,MD]=metricpoint(x,f,C,m)
 z=(x-m.center')./m.scale';dz=f./m.scale';b=[1 z z(1)^2 prod(z) z(2)^2];bd=[0 dz 2*z(1)*dz(1) z(2)*dz(1)+z(1)*dz(2) 2*z(2)*dz(2)];M=zeros(2);MD=M;
 for k=1:6,M=M+C(:,:,k)*b(k);MD=MD+C(:,:,k)*bd(k);end
end
function f=field(x,p,n)
 d=coupled(x(1),1,p,n,1);f=[p.omega*(x(2)-1),((1-d.P)/x(2)-25*(x(2)-1))/6.8];
end
function e=switcherror(g,p,n)
 d=coupled(g,1,p,n,1);e=d.rawI-p.Imax;
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





