%% Coupled steady Q-V radial limiter benchmark: new outputs only.
% Eref is the target, NOT the physical |Vc| during saturation.
% gamma is VSM reference angle minus arg(Vth), not saturated capacitor angle.
root=fileparts(mfilename('fullpath'));
out=fullfile(root,['radial_QV_consistency_' char(datetime('now','Format','yyyyMMdd_HHmmss'))]);mkdir(out);diary(fullfile(out,'console.txt'));
a=load(fullfile(root,'current_limiter_SCR5_XR10','current_limiter_SCR5_XR10_results.mat'),'results');p=a.results.parameters;n0=a.results.network;
assert(p.H==3.4&&p.K==20&&p.Kd==5&&p.D==25&&p.Imax==1.2&&p.Vmax==1.3);
x=fzero(@(z)imag(rebuild(z,p,n0).Zeq)-.350,[0 .6]);n=rebuild(x,p,n0);
Eold=[.875 .890];oldIntervals=[24.994574206146 28.493771;25.040119361334 26.628001];rows=cell(2,14);results=cell(2,1);
for k=1:2
 e0=Eold(k);gn=fzero(@(g)raw(deg2rad(g),e0,p,n,0).P-1,[20 25]);nom=raw(deg2rad(gn),e0,p,n,0);vs=e0+.1*nom.Q;
 gg=linspace(0,180,18001)';d=coupled(deg2rad(gg),vs,p,n,0);
 sw=fzero(@(g)coupled(deg2rad(g),vs,p,n,0).rawI-p.Imax,[20 30]);
 rootsSlope=find_roots(gg,d.totalPg,@(g)coupled(deg2rad(g),vs,p,n,0).totalPg);
 rootsU=find_roots(gg,abs(d.u)-p.Vmax,@(g)abs(coupled(deg2rad(g),vs,p,n,0).u)-p.Vmax);
 rootsP=find_roots(gg,d.P,@(g)coupled(deg2rad(g),vs,p,n,0).P);
 % Remove roots returned solely across the derivative discontinuity.
 rootsSlope=rootsSlope(abs(rootsSlope-sw)>1e-5);
 edges=sort(unique([0;sw;rootsSlope;rootsU;rootsP;180]));usable=[];
 for j=1:numel(edges)-1
  z=coupled(deg2rad(mean(edges(j:j+1))),vs,p,n,0);
  if z.active&&z.totalPg>0&&z.P>0&&abs(z.u)<=p.Vmax&&abs(z.ic)<=p.Imax+1e-12,usable(end+1,:)=edges(j:j+1)';end
 end
 assert(size(usable,1)==1,'Inspect multiple/absent intervals explicitly');
 activeStart=sw;activeEnd=180; % Verified active throughout the remainder.
 assert(all(d.active(gg>sw))&&all(~d.active(gg<sw)));
 nomNew=coupled(deg2rad(gn),vs,p,n,0);assert(abs(nomNew.P-1)<1e-11&&abs(nomNew.E-e0)<1e-11&&~nomNew.active);
 ss=coupled(deg2rad(linspace(usable(1)+1e-7,usable(2)-1e-7,10001)),vs,p,n,0);
 lm=(-p.D/p.mass+sqrt(complex((p.D/p.mass)^2-4*ss.totalPg/p.mass)))/2;
 lf=(-p.D/p.mass-sqrt(complex((p.D/p.mass)^2-4*ss.totalPg/p.mass)))/2;
 assert(all(real(lm)<0)&&all(real(lf)<0));
 dl=coupled(deg2rad(sw),vs,p,n,1);dr=coupled(deg2rad(sw),vs,p,n,2);
 vm=coupled(deg2rad(mean(usable)),vs,p,n,0);J=[0 1;-vm.totalPg/p.mass -p.D/p.mass];ev=eig(J);
 err=0;
 for g=linspace(usable(1)+.001,usable(2)-.001,31)
  h=1e-6;z=coupled(deg2rad(g),vs,p,n,0);zp=coupled(deg2rad(g)+h,vs,p,n,0);zm=coupled(deg2rad(g)-h,vs,p,n,0);
  err=max(err,abs((zp.P-zm.P)/(2*h)-z.totalPg));
  nodes=n.nodeCoefficients*[z.vc;p.U];assert(abs((z.vc-nodes(1))/n.Zt1+n.Yc*z.vc-z.ic)<1e-11);
 end
 assert(err<1e-7);
 % Multistart scalar root scan checks for additional algebraic roots over
 % Eref=.05..2 pu at 181 angles, in addition to continuation/Newton solution.
 maxRootCount=0;
 for g=0:180
  egrid=linspace(.05,2,1001)';r=raw(deg2rad(g)+zeros(size(egrid)),egrid,p,n,0);hh=egrid+.1*r.Q-vs;
  rr=find_roots(egrid,hh,@(e)e+.1*raw(deg2rad(g),e,p,n,0).Q-vs);
  maxRootCount=max(maxRootCount,numel(rr));assert(numel(rr)==1);
 end
 fmt=@(z)sprintf('[%.9f,%.9f]',min(z),max(z));
 rows(k,:)={e0,vs,sw,sprintf('(%.9f,%.9f)',usable),dr.P,dl.totalPg,dr.totalPg,fmt(ss.E),fmt(abs(ss.vc)),fmt(abs(ss.vg)),fmt(abs(ss.u)),fmt(ss.P),min(d.hE),fmt(real(lm))};
 results{k}=struct('Eold',e0,'Vset',vs,'oldInterval',oldIntervals(k,:),'interval',usable,'fullActiveInterval',[activeStart activeEnd],'switchLeft',dl,'switchRight',dr,'nominal',nomNew,'gammaNominal_deg',gn,'midpoint_deg',mean(usable),'midpoint',vm,'Jmidpoint',J,'eigMidpoint',ev,'gamma_deg',gg,'scan',d,'usableScan',ss,'slowEigenvalues',lm,'fastEigenvalues',lf,'derivativeError',err,'maxAlgebraicRootCount',maxRootCount);
 fprintf('Candidate Eold %.3f; Vset=%.12f; nominal gamma=%.12f, I=%.12f, u=%.12f, Vc=%.12f, PCC=%.12f\n',e0,vs,gn,abs(nomNew.ic),abs(nomNew.u),abs(nomNew.vc),abs(nomNew.vg));
 fprintf('Midpoint gamma=%.12f; Pe=%.12f; fixed-E partial=%.12f, Eprime=%.12f, total Pg=%.12f; eigJ=(%.12g,%.12g)\n',mean(usable),vm.P,vm.Pg,vm.Eg,vm.totalPg,ev);
 fprintf('Branch algebraic denominator min %.12f, max residual %.3g, derivative error %.3g\n',min(d.hE),max(abs(d.E+.1*d.Q-vs)),err);
end
T=cell2table(rows,'VariableNames',{'OldE','FixedVset','Switch_deg','StableFeasibleActive_deg','Pe_switch','TotalPg_left','TotalPg_right','ErefRange','VcMagnitudeRange','PCCMagnitudeRange','uRequiredRange','PstarRange','MinImplicitDenominator','SlowEigenvalueRange'});disp(T);
writetable(T,fullfile(out,'comparison.csv'));save(fullfile(out,'results.mat'),'results','T','p','n','x');diary off;
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



