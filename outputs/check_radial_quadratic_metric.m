root=fileparts(mfilename('fullpath'));out=fullfile(root,'radial_quadratic_robust_20261001_161903');load(fullfile(out,'metric.mat'));a=load(fullfile(root,'radial_active_transient_20261001_152918','results.mat'),'p','n','primary');p=a.p;n=a.n;
err=0;for z1=linspace(-.9,.9,7),for z2=linspace(-.9,.9,7)
 z=[z1 z2];x=center+scale.*z;d=coupled(x(1),1,p,n,0);f=[p.omega*(x(2)-1),((1-d.P)/x(2)-25*(x(2)-1))/6.8];[M,MD]=metric(C,z,f./scale);h=1e-7;MP=metric(C,z+h*f./scale,[0 0]);MM=metric(C,z-h*f./scale,[0 0]);err=max(err,max(abs((MP-MM)/(2*h)-MD),[],'all'));
end,end;assert(err<1e-6);
r=a.primary;ix=find(r.xPost(1:end-1,1)>=bounds(1,1)&r.xPost(2:end,1)<bounds(1,1),1);
exitTime=fzero(@(t)interp1(r.tPost,r.xPost(:,1),t,'pchip')-bounds(1,1),r.tPost(ix+[0 1]));
active=r.dPost.active;inside=all(r.xPost>=bounds(1,:)&r.xPost<=bounds(2,:),2);assert(all(r.xPost(active&~inside,1)<bounds(1,1)));
fprintf('Mdot error=%.12g, center trace=%.17g, active contained samples %d/%d, interpolated lower-bound exit %.12f s, limiter release %.12f s\n',err,trace(C(:,:,1)),nnz(active&inside),nnz(active),exitTime,r.currentEventsPost(end));
S=struct('MdotError',err,'traceCenter',trace(C(:,:,1)),'insideActiveSamples',nnz(active&inside),'totalActiveSamples',nnz(active),'exitTimeInterpolated_s',exitTime,'releaseTime_s',r.currentEventsPost(end),'startTime_s',r.tPost(1));fid=fopen(fullfile(out,'independent_checks.json'),'w');fwrite(fid,jsonencode(S),'char');fclose(fid);
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





