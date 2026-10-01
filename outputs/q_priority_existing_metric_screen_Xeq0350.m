%% Existing Pourazad-style construction on the unchanged q-priority model.
% Same definitions as gfm_current_limiter_SCR5_XR10.m, no optimization:
% beta=-Pe_U/(2H), c=-Pe_gammaU/(2H), Theta=[-beta 0;-c -beta].
% M=Theta'*Theta; F=Theta*J/Theta. Coefficients use THIS model's derivatives.
% Each equilibrium has its own reference Theta, frozen for the local test.
% This is NOT one numerical M valid for the whole interval, and NOT a
% state-dependent-metric proof: no Theta_dot term or invariant region claimed.
% q(gamma) is the ideal raw PI-reference proxy; rotate into the VSM dq frame,
% apply q-priority clipping, rotate back, and reconstruct all network nodes.
% No changes to saved model, old scripts/results, controls, or network tuning.
root=fileparts(mfilename('fullpath'));
a=load(fullfile(root,'current_limiter_SCR5_XR10','current_limiter_SCR5_XR10_results.mat'),'results');p=a.results.parameters;n0=a.results.network;
assert(p.H==3.4&&p.K==20&&p.Kd==5&&p.D==25&&p.mass==6.8&&p.Imax==1.2&&p.Vmax==1.3&&p.AD==.2&&p.tol==1e-9);
extraX=fzero(@(z)imag(rebuild(z,p,n0).Zeq)-.350,[0 .6],optimset('TolX',1e-13));n=rebuild(extraX,p,n0);
% Refine the supplied rounded interval endpoints using the same model.
lo=fzero(@(z)electric(deg2rad(z),p.U,p,n).Pg,[29 30],optimset('TolX',1e-12));
hi=fzero(@(z)abs(electric(deg2rad(z),p.U,p,n).u)-p.Vmax,[47 48],optimset('TolX',1e-12));
gd=linspace(lo,hi,400001)';s=screen(gd,p,n);assert(all(s.currentOK)&&all(s.active)&&all(s.partial));
passing=intervals(gd,s.margin,@(g)screen(g,p,n).margin);
% Validate refined boundaries from both sides, not just the original grid.
for k=1:size(passing,1)
 epsDeg=min(1e-5,(passing(k,2)-passing(k,1))/100);
 assert(screen(passing(k,1)+epsDeg,p,n).margin>0&&screen(passing(k,2)-epsDeg,p,n).margin>0);
 if passing(k,1)>lo+epsDeg,assert(screen(passing(k,1)-epsDeg,p,n).margin<0);end
 assert(screen(passing(k,2)+epsDeg,p,n).margin<0);
end
% Exact zero-tolerance boundaries distinguish mathematical sign from the
% retained numerical tolerances; the original 1e-9 test is the main result.
z=p;z.tol=0;zs=screen(gd,z,n);zeroPassing=intervals(gd,zs.margin,@(g)screen(g,z,n).margin);
singularRoots=roots_grid(gd,s.beta,@(g)screen(g,p,n).beta);
% Refine all sampled interior minima of lambda_max(SymF), including endpoints.
ii=find(s.lambda(2:end-1)<=s.lambda(1:end-2)&s.lambda(2:end-1)<=s.lambda(3:end))+1;
minG=[lo;hi];for j=ii',minG(end+1,1)=fminbnd(@(g)screen(g,p,n).lambda,gd(j-1),gd(j+1),optimset('TolX',1e-11));end
mins=screen(minG,p,n);[bestLambda,ib]=min(mins.lambda);bestAngle=minG(ib);
validation=struct('PeU',0,'PeGammaU_viaPg',0,'PeGammaU_viaPU',0,'matrix',0,'KCL',0,'equilibriumResidual',0);
for gdCheck=linspace(lo+.01,hi-.01,101)
 g=deg2rad(gdCheck);h=1e-6;d=electric(g,p.U,p,n);up=electric(g,p.U+h,p,n);um=electric(g,p.U-h,p,n);gp=electric(g+h,p.U,p,n);gm=electric(g-h,p.U,p,n);
 validation.PeU=max(validation.PeU,abs((up.Pe-um.Pe)/(2*h)-d.PU));
 validation.PeGammaU_viaPg=max(validation.PeGammaU_viaPg,abs((up.Pg-um.Pg)/(2*h)-d.PgU));
 validation.PeGammaU_viaPU=max(validation.PeGammaU_viaPU,abs((gp.PU-gm.PU)/(2*h)-d.PgU));
 t=screen(gdCheck,p,n);T=[-t.beta 0;-t.c -t.beta];J=[0 1;-d.Pg/p.mass -p.D/p.mass];F=T*J/T;M=T'*T;
 validation.matrix=max(validation.matrix,max(abs([max(eig((F+F')/2))-t.lambda,min(eig(M))-t.minM])));
 nodes=n.nodeCoefficients*[d.vc;p.U];validation.KCL=max(validation.KCL,abs((d.vc-nodes(1))/n.Zt1+n.Yc*d.vc-d.ic));
 validation.equilibriumResidual=max(validation.equilibriumResidual,abs((d.Pe-d.Pe)/p.mass));
end
assert(max([validation.PeU validation.PeGammaU_viaPg validation.PeGammaU_viaPU])<1e-6&&validation.matrix<1e-8&&validation.KCL<1e-11);
points=screen([30;35;40;45;47],p,n);
fprintf(' gamma0(deg)  Pe=Pstar(pu)   lambda_min(M)   lambda_max(SymF)   max Re eig(J)  |Iconv|(pu)  |u_req|(pu)  Pass\n');
for k=1:numel(points.g)
 fprintf(' %8.3f    %.9f    %+.9e    %+.9e    %+.9e   %.6f     %.6f      %d\n',points.g(k),points.Pe(k),points.minM(k),points.lambda(k),points.maxEigJ(k),points.I(k),points.u(k),points.margin(k)>0&&points.currentOK(k));
end
fprintf('\nScreened open domain: (%.12f, %.12f) deg. Existing tolerance %.1e.\n',lo,hi,p.tol);
fprintf('Passing equilibrium intervals (tol=1e-9): %s\n',fmt(passing));
fprintf('Passing equilibrium intervals (zero tolerance): %s\n',fmt(zeroPassing));
fprintf('Minimum lambda_max(SymF) = %+.12g at gamma0=%.12f deg.\n',bestLambda,bestAngle);
fprintf('Metric min eigenvalue range: [%.12g, %.12g]; singular Theta roots: %s\n',min(s.minM),max(s.minM),mat2str(singularRoots,12));
fprintf('Current magnitude range: [%.12g,%.12g] <=1.2; voltage-command range: [%.12g,%.12g] <=1.3.\n',min(s.I),max(s.I),min(s.u),max(s.u));
boundaryTable=table();
for j=1:size(passing,1)
 for side=1:2
  g=passing(j,side);v=screen(g,p,n);causes=["metric positive-definiteness tolerance","Sym(F) negative-definiteness tolerance","local-stability tolerance","voltage-command limit","limiter-activation boundary"];
  margins=[v.minM-p.tol,-v.lambda-p.tol,-v.maxEigJ-p.tol,p.Vmax-v.u,v.rawMagnitude-p.Imax];[~,which]=min(abs(margins));cause=causes(which);
  if abs(g-lo)<1e-8,cause="input domain lower endpoint: Pe_gamma=0";elseif abs(g-hi)<1e-8,cause="voltage-command limit |u|=1.3";end
  boundaryTable=[boundaryTable;table(j,side,g,cause,'VariableNames',{'interval','side','gamma_deg','condition'})];
  fprintf('Interval %d boundary %.12f deg: %s\n',j,g,cause);
 end
end
results=struct('parameters',p,'network',n,'domain_deg',[lo hi],'passing_deg',passing,'zeroTolerancePassing_deg',zeroPassing,'boundaries',boundaryTable,'samplePoints',points,'minimumLambda',bestLambda,'minimumLambdaAngle_deg',bestAngle,'validation',validation);
% results remains in the MATLAB workspace. No output files are overwritten.

function d=electric(g,U,p,n)
% Fixed network; Vth(U)=t*U. The VSM offset arg(t) is constant for U>0.
t=n.Vth/p.U;theta=g+angle(t);rot=exp(1i*theta);e=exp(-1i*g);
r=(p.E-abs(t)*U.*e)/n.Zeq;rg=1i*abs(t)*U.*e/n.Zeq;rU=-abs(t)*e/n.Zeq;rgU=1i*abs(t)*e/n.Zeq;
active=abs(r)>p.Imax;primary=imag(r);partial=active&abs(primary)<p.Imax;full=active&~partial;
s=r;sg=rg;sU=rU;sgU=rgU;
if any(partial(:))
 a=primary(partial);ag=imag(rg(partial));aU=imag(rU(partial));agU=imag(rgU(partial));signD=sign(real(r(partial)));f=sqrt(p.Imax^2-a.^2);
 fp=-a./f;fpp=-p.Imax^2./f.^3;
 s(partial)=signD.*f+1i*a;sg(partial)=signD.*fp.*ag+1i*ag;sU(partial)=signD.*fp.*aU+1i*aU;
 sgU(partial)=signD.*(fpp.*ag.*aU+fp.*agU)+1i*agU;
end
if any(full(:)),s(full)=1i*sign(primary(full))*p.Imax;sg(full)=0;sU(full)=0;sgU(full)=0;end
ic=rot.*s;icg=rot.*(sg+1i*s);icU=rot.*sU;icgU=rot.*(sgU+1i*sU);
vc=t*U+n.Zeq*ic;vcg=n.Zeq*icg;vcU=t+n.Zeq*icU;vcgU=n.Zeq*icgU;
vg=n.vgE*vc+n.vgU*U;ig=n.igE*vc+n.igU*U;
vgg=n.vgE*vcg;igg=n.igE*vcg;vgU=n.vgE*vcU+n.vgU;igU=n.igE*vcU+n.igU;vgm=n.vgE*vcgU;igm=n.igE*vcgU;
d=struct('Pe',real(vg.*conj(ig)),'Pg',real(vgg.*conj(ig)+vg.*conj(igg)), ...
 'PU',real(vgU.*conj(ig)+vg.*conj(igU)), ...
 'PgU',real(vgm.*conj(ig)+vgg.*conj(igU)+vgU.*conj(igg)+vg.*conj(igm)), ...
 'vc',vc,'ic',ic,'u',vc+(n.Zconverter+p.AD)*ic,'active',active,'partial',partial,'rawMagnitude',abs(r));
end
function s=screen(gd,p,n)
g=deg2rad(gd);d=electric(g,p.U,p,n);beta=-d.PU/p.mass;c=-d.PgU/p.mass;r=c./beta;b=p.D/p.mass;a=-d.Pg/p.mass;
off=(1+a+b*r-r.^2)/2;lm=-b/2+hypot(b/2-r,off);
m11=beta.^2+c.^2;m12=beta.*c;m22=beta.^2;mm=(m11+m22+hypot(m11-m22,2*m12))/2;mi=beta.^4./mm;
e1=(-b+sqrt(complex(b^2+4*a)))/2;e2=(-b-sqrt(complex(b^2+4*a)))/2;local=max(real(e1),real(e2));
% Equality |I|=Imax is feasible on the entire active branch, NOT a strict
% positive margin to insert into the metric's min-margin test.
currentOK=abs(d.ic)<=p.Imax+1e-12;
margin=min([mi(:)-p.tol,-lm(:)-p.tol,-local(:)-p.tol,p.Vmax-abs(d.u(:)),d.rawMagnitude(:)-p.Imax],[],2);
margin(abs(beta(:))<=1e-12|~currentOK(:))=-Inf;
s=struct('g',gd,'Pe',d.Pe,'Pg',d.Pg,'PU',d.PU,'PgU',d.PgU,'beta',beta,'c',c,'minM',mi,'lambda',lm,'maxEigJ',local,'I',abs(d.ic),'u',abs(d.u),'currentOK',currentOK,'active',d.active,'partial',d.partial,'rawMagnitude',d.rawMagnitude,'margin',margin);
end
function I=intervals(g,y,f)
ok=y>0;beg=find(diff([false;ok])==1);fin=find(diff([ok;false])==-1);I=zeros(numel(beg),2);
for j=1:numel(beg),lo=g(beg(j));hi=g(fin(j));if beg(j)>1,lo=fzero(f,g(beg(j)+[-1 0]),optimset('TolX',1e-12));end;if fin(j)<numel(g),hi=fzero(f,g(fin(j)+[0 1]),optimset('TolX',1e-12));end;I(j,:)=[lo hi];end
end
function r=roots_grid(g,y,f)
ix=find(y(1:end-1).*y(2:end)<0);r=g(abs(y)<1e-12);for k=ix',r(end+1,1)=fzero(f,g(k+[0 1]),optimset('TolX',1e-12));end;r=uniquetol(r,1e-9,'DataScale',1);
end
function t=fmt(I)
if isempty(I),t="none";return;end;s=strings(size(I,1),1);for j=1:size(I,1),s(j)=sprintf('(%.12f, %.12f)',I(j,:));end;t=strjoin(s,' U ');
end

function n=rebuild(x,p,n0)
n=n0;z1=n0.Zt1+1i*x;z2=n0.Zt2;zg=n0.Zgrid;zb=n0.Zbreaker;
y=[1/z1+1/z2+n0.Ym,-1/z2;-1/z2,1/z2+1/(zg+zb)];co=y\diag([1/z1,1/(zg+zb)]);
iv=co(2,1)/(zg+zb);iu=(co(2,2)-1)/(zg+zb);av=co(2,1)-zb*iv;au=co(2,2)-zb*iu;icv=(1-co(1,1))/z1+n0.Yc;icu=-co(1,2)/z1;
n.Zt1=z1;n.nodalY=y;n.nodeCoefficients=co;n.vgE=av;n.vgU=au;n.igE=iv;n.igU=iu;n.icE=icv;n.icU=icu;n.Zeq=1/icv;n.Vth=-icu/icv*p.U;
end
