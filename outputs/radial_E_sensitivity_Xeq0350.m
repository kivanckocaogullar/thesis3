%% Radial limiter E sensitivity. New study only; no frozen file writes.
root=fileparts(mfilename('fullpath'));
out=fullfile(root,['radial_E_sensitivity_' char(datetime('now','Format','yyyyMMdd_HHmmss'))]);mkdir(out);diary(fullfile(out,'console.txt'));
a=load(fullfile(root,'current_limiter_SCR5_XR10','current_limiter_SCR5_XR10_results.mat'),'results');p=a.results.parameters;n0=a.results.network;
assert(p.Imax==1.2&&p.Vmax==1.3&&p.D==25&&p.H==3.4&&p.K==20&&p.Kd==5);
x=fzero(@(z)imag(rebuild(z,p,n0).Zeq)-.350,[0 .6]);n=rebuild(x,p,n0);
Es=[.85 .86 .875 .89 .90 .905 .95 p.E 1 1.05 1.10 1.15 1.20 1.30];rows=cell(numel(Es),14);details=cell(size(Es));
for k=1:numel(Es)
 p.E=Es(k);co=(p.E^2+abs(n.Vth)^2-(p.Imax*abs(n.Zeq))^2)/(2*p.E*abs(n.Vth));
 sw=acosd(co);assert(isreal(sw));g=linspace(sw+1e-8,180,100001)';d=electrical(deg2rad(g),p,n,2);
 rootsPg=find_roots(g,d.Pg,@(z)electrical(deg2rad(z),p,n,2).Pg);
 rootsP=find_roots(g,d.Pe,@(z)electrical(deg2rad(z),p,n,2).Pe);
 rootsU=find_roots(g,abs(d.u)-p.Vmax,@(z)abs(electrical(deg2rad(z),p,n,2).u)-p.Vmax);
 cuts=sort(unique([sw;rootsPg;rootsP;rootsU;180]));usable=[];
 for j=1:numel(cuts)-1
  mid=electrical(deg2rad(mean(cuts(j:j+1))),p,n,2);
  if mid.Pg>0&&mid.Pe>0&&abs(mid.u)<=p.Vmax&&abs(mid.ic)<=p.Imax+1e-12,usable(end+1,:)=cuts(j:j+1)';end
 end
 left=electrical(deg2rad(sw),p,n,1);right=electrical(deg2rad(sw),p,n,2);
 gg=linspace(0,180,100001)';dd=electrical(deg2rad(gg),p,n,0);
 nominalRoots=find_roots(gg,dd.Pe-1,@(z)electrical(deg2rad(z),p,n,0).Pe-1);
 nominal=NaN;In=NaN;Un=NaN;Vset=NaN;nomOK=false;
 for z=nominalRoots'
  t=electrical(deg2rad(z),p,n,0);
  if t.Pg>0
   nominal=z;In=abs(t.ic);Un=abs(t.u);Vset=p.E+.1*t.Q;nomOK=~t.active&&In<=p.Imax&&Un<=p.Vmax;break;
  end
 end
 textInterval="none";psRange=[NaN NaN];vsRange=[NaN NaN];maxEig=NaN;voltageMargin=NaN;
 if ~isempty(usable)
  strs=strings(size(usable,1),1);pss=[];vss=[];ees=[];us=[];
  for j=1:size(usable,1)
   strs(j)=sprintf('(%.6f,%.6f)',usable(j,:));
   ts=electrical(deg2rad(linspace(usable(j,1)+1e-7,usable(j,2)-1e-7,10001)),p,n,2);
   pss=[pss ts.Pe];vss=[vss p.E+.1*ts.Q];us=[us abs(ts.u)];
   eigSlow=(-p.D/p.mass+sqrt(complex((p.D/p.mass)^2-4*ts.Pg/p.mass)))/2;ees=[ees real(eigSlow)];
  end
  psRange=[min(pss) max(pss)];vsRange=[min(vss) max(vss)];maxEig=max(ees);voltageMargin=p.Vmax-max(us);textInterval=strjoin(strs,' U ');
 end
 for z=linspace(sw+.01,179,21)
  h=1e-6;gr=deg2rad(z);v=electrical(gr,p,n,2);vp=electrical(gr+h,p,n,2);vm=electrical(gr-h,p,n,2);
  assert(abs((vp.Pe-vm.Pe)/(2*h)-v.Pg)<1e-7);
  nodes=n.nodeCoefficients*[v.vc;p.U];assert(abs((v.vc-nodes(1))/n.Zt1+n.Yc*v.vc-v.ic)<1e-11);
 end
 rows(k,:)={p.E,sw,right.Pe,left.Pg,right.Pg,textInterval,nominal,In,Un,nomOK,Vset,mat2str(psRange,6),mat2str(vsRange,6),abs(right.u)};
 details{k}=struct('E',p.E,'usable',usable,'PstarRange',psRange,'requiredVsetRange',vsRange,'maxRealEig',maxEig,'voltageMargin',voltageMargin,'positiveSlopeRoots',rootsPg,'voltageRoots',rootsU);
end
T=cell2table(rows,'VariableNames',{'E','gamma_switch_deg','Pe_switch','Pg_left','Pg_right','FeasibleStableActive_deg','gamma_nom_deg','I_nom','u_nom','NominalUnsaturatedOK','RequiredVset_nom_Qstar0','ActivePstarRange','ActiveRequiredVsetRange','u_switch'});disp(T);
% Refine the low-E window in which an unsaturated nominal and a nonempty
% positive-slope radial active interval coexist in this ideal benchmark.
Elo=fzero(@(e)switchValue(e,p,n,1),[.85 .9]);
Ehi=fzero(@(e)switchValue(e,p,n,2),[.9 .95]);
opt=optimoptions('fsolve','Display','off','FunctionTolerance',1e-12,'StepTolerance',1e-12);
[nomQV,~,flag]=fsolve(@(y)qvNominal(y,p,n),[1;deg2rad(21)],opt);assert(flag>0&&norm(qvNominal(nomQV,p,n))<1e-9);
fprintf('Low-E candidate window: %.12f < E < %.12f (nominal-unsaturated lower boundary, active-slope upper boundary).\n',Elo,Ehi);
fprintf('Unchanged Qstar=0,Vset=1 quasi-steady nominal: E=%.12f, gamma=%.12f deg.\n',nomQV(1),rad2deg(nomQV(2)));
writetable(T,fullfile(out,'comparison.csv'));save(fullfile(out,'results.mat'),'T','details','n','x','Es','Elo','Ehi','nomQV');diary off;
function y=switchValue(E,p,n,which)
 p.E=E;g=acos((E^2+abs(n.Vth)^2-(p.Imax*abs(n.Zeq))^2)/(2*E*abs(n.Vth)));
 d=electrical(g,p,n,2);if which==1,y=d.Pe-1;else,y=d.Pg;end
end
function y=qvNominal(x,p,n)
 p.E=x(1);d=electrical(x(2),p,n,0);y=[d.Pe-1;p.E-1+.1*d.Q];
end
function n=rebuild(x,p,n0)
n=n0;z1=n0.Zt1+1i*x;z2=n0.Zt2;zg=n0.Zgrid;zb=n0.Zbreaker;
y=[1/z1+1/z2+n0.Ym,-1/z2;-1/z2,1/z2+1/(zg+zb)];co=y\diag([1/z1,1/(zg+zb)]);
iv=co(2,1)/(zg+zb);iu=(co(2,2)-1)/(zg+zb);av=co(2,1)-zb*iv;au=co(2,2)-zb*iu;
icv=(1-co(1,1))/z1+n0.Yc;icu=-co(1,2)/z1;
n.Zt1=z1;n.nodalY=y;n.nodeCoefficients=co;n.vgE=av;n.vgU=au;n.igE=iv;n.igU=iu;n.icE=icv;n.icU=icu;n.Zeq=1/icv;n.Vth=-icu/icv*p.U;
end
function d=electrical(g,p,n,mode)
v=p.E*exp(1i*(g+angle(n.Vth)));q=(v-n.Vth)/n.Zeq;r=abs(q);qg=1i*v/n.Zeq;qgg=-v/n.Zeq;
active=r>p.Imax;if mode==1,active=false(size(g));elseif mode==2,active=true(size(g));end
ic=q;ig1=qg;ig2=qgg;
if any(active(:))
 a=active;qa=q(a);ra=r(a);ga=qg(a);gga=qgg(a);rg=real(conj(qa).*ga)./ra;rgg=(abs(ga).^2+real(conj(qa).*gga))./ra-rg.^2./ra;
 ic(a)=p.Imax*qa./ra;ig1(a)=p.Imax*(ga./ra-qa.*rg./ra.^2);
 ig2(a)=p.Imax*(gga./ra-2*ga.*rg./ra.^2-qa.*rgg./ra.^2+2*qa.*rg.^2./ra.^3);
end
vc=n.Vth+n.Zeq*ic;vg1=n.Zeq*ig1;vg2=n.Zeq*ig2;vc(~active)=v(~active);
vg=n.vgE*vc+n.vgU*p.U;ig=n.igE*vc+n.igU*p.U;a=n.vgE*vg1;b=n.igE*vg1;aa=n.vgE*vg2;bb=n.igE*vg2;
d=struct('q',q,'ic',ic,'vc',vc,'active',active,'Pe',real(vg.*conj(ig)), ...
 'Q',imag(vg.*conj(ig)),'Pg',real(a.*conj(ig)+vg.*conj(b)),'Pgg',real(aa.*conj(ig)+2*a.*conj(b)+vg.*conj(bb)), ...
 'u',vc+(n.Zconverter+p.AD)*ic);
end
function r=find_roots(g,y,f)
ix=find(y(1:end-1).*y(2:end)<0);r=g(abs(y)<1e-12);for j=ix',r(end+1,1)=fzero(f,g(j+[0 1]),optimset('TolX',1e-12));end;r=uniquetol(r,1e-9,'DataScale',1);
end


