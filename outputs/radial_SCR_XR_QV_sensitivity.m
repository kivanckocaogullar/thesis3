%% External SCR-X/R only; fixed physical added reactor and frozen controls.
root=fileparts(mfilename('fullpath'));out=fullfile(root,['radial_SCR_XR_QV_' char(datetime('now','Format','yyyyMMdd_HHmmss'))]);mkdir(out);diary(fullfile(out,'console.txt'));
a=load(fullfile(root,'current_limiter_SCR5_XR10','current_limiter_SCR5_XR10_results.mat'),'results');p=a.results.parameters;n0=a.results.network;
assert(p.S==2e6&&p.H==3.4&&p.K==20&&p.Kd==5&&p.D==25&&p.Imax==1.2&&p.Vmax==1.3&&p.AD==.2);
% Calibrate ONCE against saved physical SCR5/XR10 reference; never retune.
xadd=fzero(@(z)imag(rebuild(z,p,n0).Zeq)-.350,[0 .6]);ref=rebuild(xadd,p,n0);
Ladd=xadd*(p.Vlv^2/p.S)/(2*pi*p.f);
Slist=[3 5 7 10];XRlist=[3 5 10 15];rows=cell(16,9);allResults=cell(16,1);nomRows={};detailRows={};netRows=cell(16,10);k=0;
for scr=Slist
 for xr=XRlist
  k=k+1;zmag=(p.U*p.Vhv)^2/(p.S*scr);Rg=zmag/sqrt(1+xr^2);Xg=xr*Rg;Lg=Xg/(2*pi*p.f);
  nb=n0;nb.Zgrid=(Rg+1i*Xg)/(p.Vhv^2/p.S);n=rebuild(xadd,p,nb);
  assert(n.Zt1==ref.Zt1&&n.Zt2==ref.Zt2&&n.Yc==ref.Yc&&n.Ym==ref.Ym&&n.Zconverter==ref.Zconverter);
  if scr==5&&xr==10
   assert(abs(n.Zeq-ref.Zeq)<1e-8&&abs(n.Vth-ref.Vth)<1e-8);
   gr=fzero(@(g)coupled(deg2rad(g),1,p,n,0).P-1,[15 25]);dd=coupled(deg2rad(gr),1,p,n,0);
   assert(abs(dd.E-1.000698736668)<1e-8&&abs(gr-20.633904873783)<1e-6);
   fprintf('SCR5/XR10 reproduction: dZ=%.3g,dVth=%.3g, E_nom=%.12f, gamma_nom=%.12f; fixed added L=%.12g H\n',abs(n.Zeq-ref.Zeq),abs(n.Vth-ref.Vth),dd.E,gr,Ladd);
  end
  gg=linspace(0,180,36001)';d=coupled(deg2rad(gg),1,p,n,0);
  switches=find_roots(gg,d.rawI-p.Imax,@(g)coupled(deg2rad(g),1,p,n,0).rawI-p.Imax);
  assert(numel(switches)==1,'Handle multiple activation intervals explicitly');sw=switches(1);
  % Separate smooth branches before finding derivative zeros.
  cuts=[0;sw;180];slopeRoots=[];pRoots=[];uRoots=[];usable=[];
  for j=1:2
   g=linspace(cuts(j)+1e-7,cuts(j+1)-1e-7,20001)';v=coupled(deg2rad(g),1,p,n,0);
   slopeRoots=[slopeRoots;find_roots(g,v.totalPg,@(z)coupled(deg2rad(z),1,p,n,0).totalPg)];
   pRoots=[pRoots;find_roots(g,v.P,@(z)coupled(deg2rad(z),1,p,n,0).P)];
   uRoots=[uRoots;find_roots(g,abs(v.u)-p.Vmax,@(z)abs(coupled(deg2rad(z),1,p,n,0).u)-p.Vmax)];
  end
  edges=uniquetol(sort([cuts;slopeRoots;pRoots;uRoots]),1e-10,'DataScale',1);
  for j=1:numel(edges)-1
   v=coupled(deg2rad(mean(edges(j:j+1))),1,p,n,0);
   if v.active&&v.P>0&&v.totalPg>0&&abs(v.ic)<=p.Imax+1e-12&&abs(v.u)<=p.Vmax
    usable(end+1,:)=edges(j:j+1)';
   end
  end
  left=coupled(deg2rad(sw),1,p,n,1);right=coupled(deg2rad(sw),1,p,n,2);
  assert(abs(left.P-right.P)<1e-10&&abs(left.E-right.E)<1e-10);
  nominalRoots=find_roots(gg,d.P-1,@(z)coupled(deg2rad(z),1,p,n,0).P-1);
  nomStatus="none";nominals=cell(numel(nominalRoots),1);
  for jj=1:numel(nominalRoots)
   z=coupled(deg2rad(nominalRoots(jj)),1,p,n,0);stable=z.totalPg>0;feasible=abs(z.ic)<=p.Imax+1e-12&&abs(z.u)<=p.Vmax;
   if stable&&feasible
    if z.active,nomStatus="yes/stable/active";else,nomStatus="yes/stable/inactive";end
   end
   nominals{jj}=z;nomRows(end+1,:)={scr,xr,nominalRoots(jj),z.E,abs(z.vc),abs(z.vg),abs(z.ic),abs(z.u),z.totalPg,stable,feasible,z.active};
  end
  if nomStatus=="none"&&~isempty(nominalRoots),nomStatus="roots exist; none stable+feasible";end
  intervalText="none";marginText="--";width=0;pieces=cell(size(usable,1),1);margins=[];
  if ~isempty(usable)
   itxt=strings(size(usable,1),1);
   for jj=1:size(usable,1)
    ints=usable(jj,:);g=linspace(ints(1)+1e-8,ints(2)-1e-8,10001);z=coupled(deg2rad(g),1,p,n,0);
    range=@(y)[min(y) max(y)];
    % Refine every sampled local voltage maximum and include endpoints.
    umax=max(abs(z.u));iv=find(abs(z.u(2:end-1))>=abs(z.u(1:end-2))&abs(z.u(2:end-1))>=abs(z.u(3:end)))+1;
    for ix=iv,[~,val]=fminbnd(@(a)-abs(coupled(deg2rad(a),1,p,n,0).u),g(ix-1),g(ix+1));umax=max(umax,-val);end
    ends=coupled(deg2rad(ints),1,p,n,0);umax=max([umax abs(ends.u)]);margin=p.Vmax-umax;margins(end+1)=margin;
    eig1=(-p.D/p.mass+sqrt(complex((p.D/p.mass)^2-4*z.totalPg/p.mass)))/2;eig2=-p.D/p.mass-eig1;assert(all(real(eig1)<0)&&all(real(eig2)<0));
    lowerCause=cause(ints(1),sw,p,n);upperCause=cause(ints(2),sw,p,n);
    pieces{jj}=struct('gamma',g,'data',z,'eig1',eig1,'eig2',eig2,'lowerCause',lowerCause,'upperCause',upperCause,'margin',margin);
    detailRows(end+1,:)={scr,xr,ints(1),ints(2),lowerCause,upperCause,range(z.E),range(z.P),range(z.totalPg),range(abs(z.vc)),range(abs(z.vg)),range(abs(z.ic)),range(abs(z.u)),range(real(eig1)),range(real(eig2)),margin};
    itxt(jj)=sprintf('(%.9f,%.9f)',ints);
   end
   intervalText=strjoin(itxt,' U ');marginText=strjoin(compose('%.9f',margins),'; ');width=sum(usable(:,2)-usable(:,1));
  end
  err=0;kcl=0;
  for g=[sw*.2 sw*.8 sw+.01 linspace(sw+1,179,17)]
   h=1e-6;z=coupled(deg2rad(g),1,p,n,0);zp=coupled(deg2rad(g)+h,1,p,n,0);zm=coupled(deg2rad(g)-h,1,p,n,0);
   err=max(err,abs((zp.P-zm.P)/(2*h)-z.totalPg));nodes=n.nodeCoefficients*[z.vc;p.U];kcl=max(kcl,abs((z.vc-nodes(1))/n.Zt1+n.Yc*z.vc-z.ic));
  end
    % Separate one-sided derivative checks at the switching equilibrium.
  for mode=1:2
   h=1e-6;sgn=2*mode-3;g=deg2rad(sw);z=coupled(g,1,p,n,mode);z1=coupled(g+sgn*h,1,p,n,mode);z2=coupled(g+sgn*2*h,1,p,n,mode);
   err=max(err,abs(sgn*(-3*z.P+4*z1.P-z2.P)/(2*h)-z.totalPg));
  end
  assert(err<1e-7&&kcl<1e-11&&min(d.hE)>.1);
  % Numerical root multiplicity check over E=.05..2 pu on 37 angles.
  for g=0:5:180
   ee=linspace(.05,2,1001)';zz=raw(deg2rad(g)+zeros(size(ee)),ee,p,n,0);rr=find_roots(ee,ee+.1*zz.Q-1,@(e)e+.1*raw(deg2rad(g),e,p,n,0).Q-1);assert(numel(rr)==1);
  end
  rows(k,:)={scr,xr,imag(n.Zeq),sw,right.totalPg,intervalText,width,marginText,nomStatus};
  netRows(k,:)={scr,xr,Rg,Lg,real(n.Zeq),imag(n.Zeq),real(n.Vth),imag(n.Vth),xadd,Ladd};
  allResults{k}=struct('SCR',scr,'XR',xr,'network',n,'Rgrid',Rg,'Lgrid',Lg,'switch',sw,'left',left,'right',right,'intervals',usable,'pieces',{pieces},'gamma_deg',gg,'scan',d,'nominalRoots',nominalRoots,'nominals',{nominals},'derivativeError',err,'KCLerror',kcl,'slopeRoots',slopeRoots,'powerRoots',pRoots,'voltageRoots',uRoots);
 end
end
T=cell2table(rows,'VariableNames',{'SCR','XR','Xeq','gammaSwitch_deg','TotalPg_afterSwitch','FeasibleActive_deg','Width_deg','MinVoltageMargin','NominalStatus'});disp(T);
N=cell2table(nomRows,'VariableNames',{'SCR','XR','gamma_deg','Eref','Vc','PCC','Iconv','uRequired','TotalPg','Stable','Feasible','Active'});
W=cell2table(netRows,'VariableNames',{'SCR','XR','Rgrid_ohm','Lgrid_H','Req','Xeq','VthRe','VthIm','FixedAddedX_LV','FixedAddedL_H'});
writetable(T,fullfile(out,'comparison.csv'));writetable(N,fullfile(out,'nominal_points.csv'));writetable(W,fullfile(out,'networks.csv'));
if ~isempty(detailRows)
 D=cell2table(detailRows,'VariableNames',{'SCR','XR','Lower_deg','Upper_deg','LowerCondition','UpperCondition','ErefRange','PstarRange','TotalPgRange','VcRange','PCCRange','IconvRange','uRange','eig1Range','eig2Range','VoltageMargin'});disp(D);writetable(D,fullfile(out,'feasible_intervals.csv'));
else,D=table();fprintf('No feasible positive-width limiter-active intervals in all 16 cases.\n');end
% Additional numerical mechanism and physical-voltage summary, no new sweep.
mech=cell(16,13);
for k=1:16
 r=allResults{k};n=r.network;d=r.scan;g=r.gamma_deg;mask=d.active&d.P>0;
 pgmax=max([r.right.totalPg;d.totalPg(mask)]);
 % Refine all sampled interior local maxima of total slope in exporting active branch.
 inds=find(mask(2:end-1)&d.totalPg(2:end-1)>=d.totalPg(1:end-2)&d.totalPg(2:end-1)>=d.totalPg(3:end))+1;
 for ii=inds'
  [~,val]=fminbnd(@(a)-coupled(deg2rad(a),1,p,n,0).totalPg,g(ii-1),g(ii+1),optimset('TolX',1e-10));pgmax=max(pgmax,-val);
 end
 for pp=r.powerRoots'
  if pp>r.switch,pgmax=max(pgmax,coupled(deg2rad(pp),1,p,n,0).totalPg);end
 end
 nr=[];for jj=1:numel(r.nominals),v=r.nominals{jj};if v.totalPg>0&&abs(v.u)<=p.Vmax,nr=v;break;end,end
 mech(k,:)={r.SCR,r.XR,pgmax,r.left.totalPg,r.right.totalPg,p.Vmax-abs(r.right.u),abs(nr.vc),abs(nr.vg),abs(nr.ic),abs(nr.u),r.derivativeError,r.KCLerror,min(d.hE)};
end
B=cell2table(mech,'VariableNames',{'SCR','XR','MaxTotalPg_active_export','Pg_switch_left','Pg_switch_right','VoltageMargin_switch','NominalVc','NominalPCC','NominalI','NominalU','DerivativeError','KCLerror','MinImplicitDenominator'});disp(B);writetable(B,fullfile(out,'mechanism_and_voltages.csv'));
save(fullfile(out,'results.mat'),'allResults','T','N','W','D','p','xadd','Ladd','ref','B');diary off;
function c=cause(g,sw,p,n)
 if abs(g-sw)<1e-6,c="current activation";return;end
 z=coupled(deg2rad(g),1,p,n,0);[~,ii]=min(abs([z.totalPg,z.P,abs(z.u)-p.Vmax]));names=["total Pe_gamma=0","Pe=0","voltage limit"];c=names(ii);
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





