%% Common constant metric, frozen q-priority network. New outputs only.
% Numerical enclosure evidence, NOT interval-arithmetic/SOS certification.
% q(gamma) is the ideal raw PI-reference proxy, rotated into VSM dq before
% q-priority limiting. The physical PI/current-loop dynamics are not included.
root=fileparts(mfilename('fullpath'));
out=fullfile(root,['q_priority_common_metric_Xeq0350_' char(datetime('now','Format','yyyyMMdd_HHmmss'))]);
assert(~isfolder(out),'Refusing to overwrite previous results');mkdir(out);
diary(fullfile(out,'console.txt'));
a=load(fullfile(root,'current_limiter_SCR5_XR10','current_limiter_SCR5_XR10_results.mat'),'results');p=a.results.parameters;n0=a.results.network;
assert(p.H==3.4&&p.D==25&&p.Imax==1.2&&p.Vmax==1.3);
xadd=fzero(@(x)imag(rebuild(x,p,n0).Zeq)-.350,[0 .6]);n=rebuild(xadd,p,n0);
sw=fzero(@(g)electric(deg2rad(g),p.U,p,n).rawMagnitude-p.Imax,[24 25]);
turn=fzero(@(g)electric(deg2rad(g),p.U,p,n).Pg,[29 30]);
vbound=fzero(@(g)abs(electric(deg2rad(g),p.U,p,n).u)-p.Vmax,[47 48]);
before=[fzero(@(g)screen(g,p,n).lambda,[20 21]),sw;fzero(@(g)screen(g,p,n).lambda,[45 46]),vbound];
assert(max(abs(before(:)-[20.167968;45.459623897403;24.421813;47.521506]))<2e-6);
fprintf('Original metric equilibrium screening (not nonlinear domains):\n');disp(before);
fprintf('Active stable feasible interval: (%.12f, %.12f) deg\n',turn,vbound);
% Explicit search design: lambda=1e-5 / normalized time, nu halfwidth=.01.
% Refine lower endpoints of common boxes, maintaining .02 deg boundary guards.
% Every selected Pstar is fixed; vertex enclosure covers ALL independent
% Pstar in Pe(family), gamma in box and nu in box (no Pstar=Pe(gamma) substitution).
rate=1e-5;nuBounds=[.99 1.01];guard=.02;
upper=[sw-guard,vbound-guard];lowerSeed=[1,30];searchFloor=[0,turn+1e-5];
resCell=cell(1,2);
for branch=1:2
 hi=lowerSeed(branch);lo=searchFloor(branch);
 test=fitbox([hi upper(branch)],nuBounds,rate,p,n);
 assert(test.margin>1e-10,'Seed box failed; quadratic metric required');
 floorTest=fitbox([lo upper(branch)],nuBounds,rate,p,n);
 if floorTest.margin>1e-10,hi=lo;else
  for j=1:24
   mid=(lo+hi)/2;trial=fitbox([mid upper(branch)],nuBounds,rate,p,n);
   if trial.margin>1e-10,hi=mid;else,lo=mid;end
  end
 end
 % Back away from numerical feasibility edge to obtain a robust slack.
 selectedLower=hi+.02;
 r=fitbox([selectedLower upper(branch)],nuBounds,rate,p,n);
 r.searchLower=hi;r.family=[selectedLower+guard upper(branch)-guard];
 r.branch=branch;r.rate=rate;r.nuBounds=nuBounds;
 % Independent much denser verification of electrical extrema and LMI.
 dd=electric(deg2rad(linspace(r.box(1),r.box(2),100001)),p.U,p,n);
 assert(all(dd.Pg>0));assert(all(abs(dd.u)<p.Vmax));
 if branch==1,assert(all(~dd.active));else,assert(all(dd.partial));end
 r.currentMargin=p.Imax-max(abs(dd.ic));r.voltageMargin=p.Vmax-max(abs(dd.u));
 r.minM=min(eig(r.M));r.worstQ=-r.margin;
 r.aDense=[min(-dd.Pg/(p.mass*nuBounds(1))) max(-dd.Pg/(p.mass*nuBounds(2)))];
 assert(r.aDense(1)>=r.a(1)&&r.aDense(2)<=r.a(2));
 % Invariant ellipsoid contained in box: e'Me <= rho. Half the squared
 % distance to each rectangle face leaves room for all family members.
 invM=inv(r.M);r.rho=.5*min([deg2rad(guard)^2/invM(1,1),.01^2/invM(2,2)]);
 r.axisGamma_deg=rad2deg(sqrt(r.rho*invM(1,1)));r.axisNu=sqrt(r.rho*invM(2,2));
  r.PstarFamily=electric(deg2rad(r.family),p.U,p,n).Pe;
 R=chol(r.M);r.numericalRateBound=Inf;
 for av=r.a,for bv=r.b
  J=[0 1;av bv];Q=J'*r.M+r.M*J;
  r.numericalRateBound=min(r.numericalRateBound,-max(eig(R'\Q/R))/2);
 end,end
 % Independently verify the off-equilibrium Jacobian at FIXED Pstar.
 r.jacobianError=0;
 for gg=linspace(r.box(1),r.box(2),11)
  for vv=[.99 1 1.01]
   for ps=r.PstarFamily
    xx=[deg2rad(gg);vv];ddx=electric(xx(1),p.U,p,n);
    J=[0 1;-ddx.Pg/(p.mass*vv),(-(ps-ddx.Pe)/vv^2-p.D)/p.mass];
    JJ=zeros(2);h=1e-6;
    for col=1:2,ee=zeros(2,1);ee(col)=h;JJ(:,col)=(dynamics(xx+ee,ps,p,n)-dynamics(xx-ee,ps,p,n))/(2*h);end
    r.jacobianError=max(r.jacobianError,max(abs(J-JJ),[],'all'));
   end
  end
 end
 assert(r.jacobianError<1e-7&&r.numericalRateBound>rate);
 resCell{branch}=r;
 fprintf('Box=[%.12f,%.12f] deg; family=[%.12f,%.12f] deg; search lower=%.12f deg\n',r.box,r.family,r.searchLower);
 fprintf('minM=%.12g; worstQ=%.12g; numerical rate bound=%.12g; I margin=%.12g; voltage margin=%.12g\n',r.minM,r.worstQ,r.numericalRateBound,r.currentMargin,r.voltageMargin);
 fprintf('Invariant ellipsoid rho=%.12g; coordinate extents=(%.12g deg, %.12g nu); fixed-P Jacobian error=%.3g\n',r.rho,r.axisGamma_deg,r.axisNu,r.jacobianError);
 fprintf('Branch %d common M (trace=1):\n',branch);disp(r.M);
end
res=[resCell{:}];
T=table(["Inactive";"Active"],before(:,1),before(:,2),vertcat(res.family), 'VariableNames',{'Branch','OriginalLower_deg','OriginalUpper_deg','NewFamily_deg'});
T.GammaBox_deg=vertcat(res.box);T.NuBox=repmat(nuBounds,2,1);T.minEigM=[res.minM]';T.WorstQEigenvalue=[res.worstQ]';T.Lambda=repmat(rate,2,1);T.CurrentMargin=[res.currentMargin]';T.VoltageMargin=[res.voltageMargin]';
T.NumericalRateBound=[res.numericalRateBound]';
disp(T);writetable(T,fullfile(out,'comparison.csv'));
save(fullfile(out,'results.mat'),'p','n','xadd','before','sw','turn','vbound','res','T');
fprintf('Degree 0 succeeded: Mdot=0 exactly; quadratic optimization unnecessary for this enlargement.\n');
fprintf('Numerical common-Jacobian vertex test; electrical extrema use dense samples and refined extrema, not directed rounding.\n');
fprintf('Therefore numerical nonlinear-neighbourhood evidence, NOT a mathematically uniform certificate.\n');
fprintf('Separate common M per smooth branch; no switched-system certificate or bridging.\n');
fprintf('Invariant ellipsoid construction is conditional on the numerical domain inequalities holding everywhere.\n');
diary off;

function f=dynamics(x,ps,p,n)
 d=electric(x(1),p.U,p,n);
 f=[x(2)-1;((ps-d.Pe)/x(2)-p.D*(x(2)-1))/p.mass];
end
function r=fitbox(box,nb,lam,p,n)
 g=linspace(box(1),box(2),4001);d=electric(deg2rad(g),p.U,p,n);
 % Refine all detected stationary extrema; retain conservative floating pads.
 pg=extrema(g,d.Pg,@(z)electric(deg2rad(z),p.U,p,n).Pg);
 pe=extrema(g,d.Pe,@(z)electric(deg2rad(z),p.U,p,n).Pe);
 aa=[-pg(2)/(p.mass*nb(1)) -pg(1)/(p.mass*nb(2))]+[-1 1]*1e-10;
 delta=pe(2)-pe(1);bb=[(-p.D-delta/nb(1)^2)/p.mass,(-p.D+delta/nb(1)^2)/p.mass]+[-1 1]*1e-10;
 Js=zeros(2,2,4);k=0;for a=aa,for b=bb,k=k+1;Js(:,:,k)=[0 1;a b];end,end
 opt=optimoptions('fmincon','Display','off','Algorithm','sqp','ConstraintTolerance',1e-12,'OptimalityTolerance',1e-12,'MaxIterations',400);
 best=-Inf;
 for sh=[.15 .25]
  M0=[1 sh;sh .1];M0=M0/trace(M0);z0=[M0(1,1);M0(1,2);-1e-4];
  [z,~,flag]=fmincon(@(z)-z(3),z0,[],[],[],[],[1e-6;-1;-1],[1-1e-6;1;1],@(z)con(z,Js,lam),opt);
  M=[z(1) z(2);z(2) 1-z(1)];margin=Inf;
  for k=1:4,J=Js(:,:,k);margin=min(margin,-max(eig(J'*M+M*J+2*lam*M)));end
  if min(eig(M))>1e-6&&margin>best,best=margin;chosen=M;exitflag=flag;end
 end
 if ~exist('chosen','var'),chosen=eye(2)/2;best=-Inf;exitflag=-99;end
 r=struct('box',box,'M',chosen,'margin',best,'a',aa,'b',bb,'PeRange',pe,'PgRange',pg,'exitflag',exitflag);
end
function [c,ceq]=con(z,Js,lam)
 M=[z(1) z(2);z(2) 1-z(1)];c=1e-6-min(eig(M));
 for k=1:4,J=Js(:,:,k);c(end+1)=max(eig(J'*M+M*J+2*lam*M))+z(3);end
 ceq=[];
end
function bounds=extrema(g,y,f)
 vals=[y(1) y(end)];
 for signF=[1 -1]
  inds=find(signF*y(2:end-1)<=signF*y(1:end-2)&signF*y(2:end-1)<=signF*y(3:end))+1;
  for k=inds,[~,v]=fminbnd(@(z)signF*f(z),g(k-1),g(k+1),optimset('TolX',1e-11));vals(end+1)=signF*v;end
 end
 bounds=[min([y vals])-1e-9 max([y vals])+1e-9];
end
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



