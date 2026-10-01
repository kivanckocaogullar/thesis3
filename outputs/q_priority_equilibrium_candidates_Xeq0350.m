%% Five equilibrium candidates: same physical network and q-priority proxy.
% Read-only inputs; no Simulink load/edit, metric optimization or file export.
% q(gamma)=(E exp(j theta)-Vth)/Zeq is an IDEAL raw-current proxy, NOT the
% actual voltage PI reference (which includes PI memory/feedforward).
% theta=gamma+arg(Vth); r_dq=exp(-j theta)*q; apply q priority to r_dq;
% Iconv=exp(j theta)*s_dq; Vc=Vth+Zeq*Iconv. Reconstruct grid power afterward.
% At EACH row nu0=1 and Pstar=Pe(gamma0); these are different power setpoints.
root=fileparts(mfilename('fullpath'));
ref=load(fullfile(root,'current_limiter_SCR5_XR10','current_limiter_SCR5_XR10_results.mat'),'results');
p=ref.results.parameters;n0=ref.results.network;
assert(p.H==3.4&&p.K==20&&p.Kd==5&&p.D==25&&p.Imax==1.2&&p.Vmax==1.3&&p.AD==.2);
extraX=fzero(@(z)imag(rebuild(z,p,n0).Zeq)-.350,[0 .6],optimset('TolX',1e-13));n=rebuild(extraX,p,n0);
assert(abs(imag(n.Zeq)-.350)<1e-12&&n.Zgrid==n0.Zgrid);
angles=[30;35;40;45;47];summary=zeros(5,10);states=cell(5,1);
validation=struct('frameToDQ',0,'frameToNetwork',0,'networkKCL',0,'power',0,'PeDerivative',0,'Jacobian',0,'equilibriumResidual',0);
for k=1:5
 g=deg2rad(angles(k));theta=g+angle(n.Vth);d=map(g,p,n,3,0);ps=d.Pe;
 % Independent clip implementation checks the rotating-frame priority map.
 rawGlobal=(p.E*exp(1i*theta)-n.Vth)/n.Zeq;r=rawGlobal*exp(-1i*theta);
 sq=max(-p.Imax,min(p.Imax,imag(r)));cap=sqrt(max(0,p.Imax^2-sq^2));sd=max(-cap,min(cap,real(r)));s=sd+1i*sq;
 validation.frameToDQ=max(validation.frameToDQ,abs(r-d.raw));validation.frameToNetwork=max(validation.frameToNetwork,abs(exp(1i*theta)*s-d.ic));
 nodes=n.nodeCoefficients*[d.vc;p.U];validation.networkKCL=max(validation.networkKCL,abs((d.vc-nodes(1))/n.Zt1+n.Yc*d.vc-d.ic));
 validation.power=max(validation.power,abs(d.Pe-real(d.Vgrid*conj(d.Igrid))));
 J=[0 1;-d.Pg/p.mass -p.D/p.mass];e=sort(real(eig(J)),'descend');assert(all(abs(imag(eig(J)))<1e-12));
 h=1e-6;dp=map(g+h,p,n,3,0);dm=map(g-h,p,n,3,0);validation.PeDerivative=max(validation.PeDerivative,abs((dp.Pe-dm.Pe)/(2*h)-d.Pg));
 xe=[g;1];Jn=zeros(2);for j=1:2,hh=zeros(2,1);hh(j)=h;Jn(:,j)=(rhs(xe+hh,ps,p,n)-rhs(xe-hh,ps,p,n))/(2*h);end
 validation.Jacobian=max(validation.Jacobian,norm(J-Jn,inf));validation.equilibriumResidual=max(validation.equilibriumResidual,norm(rhs(xe,ps,p,n),inf));
 currentOK=abs(d.ic)<=p.Imax+1e-12;voltageOK=abs(d.u)<=p.Vmax;candidate=currentOK&&voltageOK&&d.Pg>0&&all(e<0)&&d.active;
 summary(k,:)=[angles(k),ps,abs(d.ic),abs(d.u),p.Vmax-abs(d.u),d.Pg,e',currentOK&&voltageOK,candidate];
 states{k}=struct('electrical',d,'rawProxyVSM',r,'limitedReferenceVSM',s,'Pstar',ps,'J',J,'eigenvalues',e,'candidate',candidate);
end
assert(validation.frameToDQ<1e-12&&validation.frameToNetwork<1e-12&&validation.networkKCL<1e-11&&validation.power<1e-12);
assert(validation.PeDerivative<1e-7&&validation.Jacobian<1e-7&&validation.equilibriumResidual<1e-12);
fprintf(' gamma0    Pe=Pstar    |Iconv|    |u_req|    1.3-|u|    Pe_gamma      lambda1       lambda2    I/u OK  Candidate\n');
fprintf('  (deg)      (pu)       (pu)      (pu)       (pu)      (pu/rad)      normalized analytical time^-1      (0/1)    (0/1)\n');
for k=1:5,fprintf(' %5.0f   %.9f  %.6f  %.6f  %.6f  %+.9f  %+.9f  %+.9f     %d       %d\n',summary(k,:));end

function dx=rhs(x,ps,p,n)
d=map(x(1),p,n,3,0);dx=[x(2)-1;((ps-d.Pe)/x(2)-p.D*(x(2)-1))/p.mass];
end

function n=rebuild(x,p,n0)
n=n0;z1=n0.Zt1+1i*x;z2=n0.Zt2;zg=n0.Zgrid;zb=n0.Zbreaker;
y=[1/z1+1/z2+n0.Ym,-1/z2;-1/z2,1/z2+1/(zg+zb)];co=y\diag([1/z1,1/(zg+zb)]);
iv=co(2,1)/(zg+zb);iu=(co(2,2)-1)/(zg+zb);av=co(2,1)-zb*iv;au=co(2,2)-zb*iu;icv=(1-co(1,1))/z1+n0.Yc;icu=-co(1,2)/z1;
n.Zt1=z1;n.nodalY=y;n.nodeCoefficients=co;n.vgE=av;n.vgU=au;n.igE=iv;n.igU=iu;n.icE=icv;n.icU=icu;n.Zeq=1/icv;n.Vth=-icu/icv*p.U;
end
function d=map(g,p,n,strategy,mode)
theta=g+angle(n.Vth);rot=exp(1i*theta);v=p.E*rot;q=(v-n.Vth)/n.Zeq;
% r(gamma)=(E-|Vth| exp(-j gamma))/Zeq is the raw-reference proxy
% expressed in the VSM frame, not in the global network phasor frame.
r=(p.E-abs(n.Vth)*exp(-1i*g))/n.Zeq;rg=1i*abs(n.Vth)*exp(-1i*g)/n.Zeq;rgg=abs(n.Vth)*exp(-1i*g)/n.Zeq;
active=abs(r)>p.Imax;if mode==1,active=false(size(g));elseif mode==2,active=true(size(g));end
s=r;sg=rg;sgg=rgg;
if strategy==1
 a=active;if any(a(:)),ra=abs(r(a));qa=r(a);ga=rg(a);gga=rgg(a);ag=real(conj(qa).*ga)./ra;agg=(abs(ga).^2+real(conj(qa).*gga))./ra-ag.^2./ra;
 s(a)=p.Imax*qa./ra;sg(a)=p.Imax*(ga./ra-qa.*ag./ra.^2);sgg(a)=p.Imax*(gga./ra-2*ga.*ag./ra.^2-qa.*agg./ra.^2+2*qa.*ag.^2./ra.^3);end
else
 a=real(r);ag=real(rg);agg=real(rgg);b=imag(r);if strategy==3,a=imag(r);ag=imag(rg);agg=imag(rgg);b=real(r);end
 full=active & abs(a)>=p.Imax;partial=active & ~full;
 sa=a;sag=ag;sagg=agg;sb=b;sbg=zeros(size(g));sbgg=sbg;
 sa(full)=sign(a(full))*p.Imax;sag(full)=0;sagg(full)=0;sb(full)=0;
 f=sqrt(max(0,p.Imax^2-a(partial).^2));sig=sign(b(partial));
 sb(partial)=sig.*f;sbg(partial)=-sig.*a(partial).*ag(partial)./f;
 sbgg(partial)=-sig.*((ag(partial).^2+a(partial).*agg(partial))./f+a(partial).^2.*ag(partial).^2./f.^3);
 if strategy==2,s(active)=sa(active)+1i*sb(active);sg(active)=sag(active)+1i*sbg(active);sgg(active)=sagg(active)+1i*sbgg(active);
 else,s(active)=sb(active)+1i*sa(active);sg(active)=sbg(active)+1i*sag(active);sgg(active)=sbgg(active)+1i*sagg(active);end
end
ic=rot.*s;icg=rot.*(sg+1i*s);icgg=rot.*(sgg+2i*sg-s);vc=n.Vth+n.Zeq*ic;vcg=n.Zeq*icg;vcgg=n.Zeq*icgg;
vg=n.vgE*vc+n.vgU*p.U;ig=n.igE*vc+n.igU*p.U;vg1=n.vgE*vcg;ig1=n.igE*vcg;vg2=n.vgE*vcgg;ig2=n.igE*vcgg;
d=struct('raw',r,'q',q,'ic',ic,'vc',vc,'Vgrid',vg,'Igrid',ig,'active',active,'Pe',real(vg.*conj(ig)), ...
 'Pg',real(vg1.*conj(ig)+vg.*conj(ig1)),'Pgg',real(vg2.*conj(ig)+2*vg1.*conj(ig1)+vg.*conj(ig2)), ...
 'u',vc+(n.Zconverter+p.AD)*ic,'ug',vcg+(n.Zconverter+p.AD)*icg);
end
