%% Reevaluate five saved FCT cases against the widened, bounded metrics.
% No ODE simulation, CCT search, or metric optimization is performed.
root=fileparts(mfilename('fullpath'));project=fileparts(root);out=fullfile(root,'results');addpath(root);
a=load(fullfile(root,'inputs','saved_fct_analysis.mat'));
cctRoot=fullfile(project,'CCT_analysis_20261011_181410');regionRoot=fullfile(project,'contraction region');
b=load(fullfile(cctRoot,'baseline_results.mat'));c=load(fullfile(cctRoot,'model_inputs.mat'),'p','n');w=load(fullfile(regionRoot,'inputs','model_inputs.mat'),'p','n');
p=a.p;n=a.n;c.p.Pstar=1;c.p.omega=100*pi;w.p.Pstar=1;w.p.omega=100*pi;
assert(isequaln(p,b.p)&&isequaln(n,b.n)&&isequaln(p,c.p)&&isequaln(n,c.n)&&isequaln(p,w.p)&&isequaln(n,w.n),'Physical model mismatch: do not reuse cache');
fields={'tEnd','RelTol','AbsTol','faultMaxStep','postMaxStep'};for k=1:numel(fields),assert(isequaln(a.cfg.(fields{k}),b.cfg.(fields{k})));end
assert(abs(a.cfg.gammaEq-b.cfg.gammaEq)<1e-12);assert(a.Uret==.4&&isequal(a.FCT,[.05 .10 .14 .15 .25]));
assert(norm(a.x0-[a.cfg.gammaEq;1],inf)<1e-12);assert(a.pf.U==a.p.U*.4);assert(isequaln(rebuild(0,a.pf,n),a.nf));
runs=a.runs;results=a.results;cfg=a.cfg;FCT=a.FCT;x0=a.x0;Uret=a.Uret;metrics=cell(2,1);oldMetrics=metrics;certificates=metrics;
names={'inactive','active'};boundsTable=table();
for j=1:2
 dirMetric=fullfile(regionRoot,'metrics',names{j});m=load(fullfile(dirMetric,'metric.mat'));cert=jsondecode(fileread(fullfile(dirMetric,'interval_verification.json')));jm=jsondecode(fileread(fullfile(dirMetric,'metric.json')));
 assert(strcmp(cert.status,'VERIFIED')&&cert.unresolved_boxes==0&&cert.counterexample_boxes==0);
 assert(max(abs(m.C-permute(jm.coefficients,[2 3 1])),[],'all')<1e-12&&max(abs(m.bounds-jm.bounds),[],'all')<1e-12);
 assert(max(abs(rad2deg(m.bounds(:,1))-cert.gamma_bounds_deg(:)))<1e-10&&max(abs(m.bounds(:,2)-cert.nu_bounds(:)))<1e-12);
 assert(jm.U==p.U&&jm.mass==p.mass&&jm.D==p.D&&jm.Imax==p.Imax&&jm.umax==p.Vmax&&jm.AD==p.AD&&abs(jm.omega-p.omega)<1e-12);
 nflds=fieldnames(jm.network);for q=1:numel(nflds),z=jm.network.(nflds{q});assert(abs(n.(nflds{q})-complex(z(1),z(2)))<1e-12);end
 metrics{j}=m;certificates{j}=cert;oldMetrics{j}=load(fullfile(regionRoot,'reference_old_metrics',names{j},'metric.mat'));
 boundsTable=[boundsTable;table(string(names{j}),rad2deg(m.bounds(1,1)),rad2deg(m.bounds(2,1)),m.bounds(1,2),m.bounds(2,2),m.lambda,m.upper,cert.epsilon_margin_holds_in_all_boxes,'VariableNames',{'Mode','GammaMin_deg','GammaMax_deg','NuMin_pu','NuMax_pu','Rate_per_s','MetricUpper','UniformEpsilonMarginCertified'})];
end
comparison=table();
for k=1:5
 x=[deg2rad(results.gamma_clear_deg(k)),results.nu_clear(k)];r=runs{k};assert(max(abs(r.postfault.x(1,:)-x))<1e-12&&abs(r.postfault.t(1)-FCT(k))<1e-12);
 d=coupled(x(1),1,p,n,0);j=1+double(d.active);m=metrics{j};old=oldMetrics{j};inside=all(x>=m.bounds(1,:)&x<=m.bounds(2,:));oldInside=all(x>=old.bounds(1,:)&x<=old.bounds(2,:));
 [mn,mx,ms]=metric_at(x,m,p,d);voltage=abs(d.u)<=p.Vmax;rate=mn>0&&mx<m.upper&&ms<-m.lambda&&voltage;
 omega=p.omega*x(2);fHz=50*x(2);
 comparison=[comparison;table(FCT(k),x(1),rad2deg(x(1)),x(2),omega,fHz,string(names{j}),oldInside,inside,mn,mx,mx/mn,ms,rate,voltage,results.classification(k),'VariableNames',{'FCT_s','GammaClear_rad','GammaClear_deg','NuClear_pu','OmegaClear_rad_s','FrequencyClear_Hz','PostfaultMode','InsideOldRectangle','InsideWidenedRectangle','MinEigM','MaxEigM','ConditionM','MaxEigSymF','PointwiseBoundedRatePass','VoltageValid','NumericalRecovery'})];
end
cct=readtable(fullfile(cctRoot,'results_refined_20261011_184900_867','cct_results.csv'));cct=cct(abs(cct.Uret-.4)<1e-12,:);assert(height(cct)==1);writetable(cct,fullfile(out,'cct_reference.csv'));
writetable(comparison,fullfile(out,'clearing_state_comparison.csv'));writetable(boundsTable,fullfile(out,'region_bounds.csv'));writetable(a.clearingTable,fullfile(out,'clearing_states.csv'));
writetable(results,fullfile(out,'assessment_details.csv'));writetable(table(true,true,true,true,'VariableNames',{'CachedModelMatchesCCT','CachedModelMatchesRegions','SolverSettingsMatch','FaultScenarioMatches'}),fullfile(out,'cache_validation.csv'));
save(fullfile(out,'widened_fct_analysis.mat'),'runs','results','cfg','FCT','x0','Uret','p','n','metrics','oldMetrics','certificates','comparison','cct','-v7.3');
disp(comparison);disp(cct);
plot_saved_widened_FCT(out);plot_widened_containment(out);

function [minM,maxM,maxS]=metric_at(x,m,p,d)
z=(x-m.center)./m.scale;f=[p.omega*(x(2)-1),((1-d.P)/x(2)-p.D*(x(2)-1))/p.mass];dz=f./m.scale;
b=[1 z z(1)^2 z(1)*z(2) z(2)^2];bd=[0 dz 2*z(1)*dz(1) z(2)*dz(1)+z(1)*dz(2) 2*z(2)*dz(2)];
M=zeros(2);MD=M;for j=1:6,M=M+b(j)*m.C(:,:,j);MD=MD+bd(j)*m.C(:,:,j);end
ev=eig((M+M')/2);minM=min(ev);maxM=max(ev);maxS=NaN;
[R,flag]=chol(M);if flag==0
 J=[0 p.omega;-d.totalPg/(p.mass*x(2)),(-(1-d.P)/x(2)^2-p.D)/p.mass];
 S=(R'\(J'*M+M*J+MD)/R)/2;maxS=max(eig((S+S')/2));
end
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



















