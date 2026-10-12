%% CCT per retained voltage: unchanged saved GFM physics; status-aware CCT visualization.
% Every execution creates a fresh directory; results never overwrite a prior run.
% NOT_RECOVERED means failure of the finite-horizon recovery test, NOT proven instability.
% Optional plotting-only mode: set cctPlotResultsDir to an existing run directory.
% Example: cctPlotResultsDir=fullfile(pwd,'results_...'); run this script.
if exist('cctRefineResultsDir','var') && strlength(string(cctRefineResultsDir))>0
 refine_endpoint_case(char(cctRefineResultsDir));
 clear cctRefineResultsDir;return;
end
if exist('cctPlotResultsDir','var') && strlength(string(cctPlotResultsDir))>0
 loadedCCT=load(fullfile(cctPlotResultsDir,'cct_results.mat'),'results');
 plot_cct_results(loadedCCT.results,char(cctPlotResultsDir));
 clear cctPlotResultsDir loadedCCT;
 return;
end
clearvars;clc;root=fileparts(mfilename('fullpath'));
runID=char(datetime('now','Format','yyyyMMdd_HHmmss_SSS'));out=fullfile(root,['results_' runID]);assert(~isfolder(out));mkdir(out);
a=load(fullfile(root,'model_inputs.mat'),'p','n');base=load(fullfile(root,'baseline_results.mat'),'p','n','cfg','x0','summary');
p=a.p;n=a.n;p.Pstar=1;p.omega=100*pi;assert(isequaln(p,base.p)&&isequaln(n,base.n),'Saved baseline/input mismatch');
assert(p.Imax==1.2&&p.Vmax==1.3&&p.mass==6.8&&p.D==25);
cfg=base.cfg;cfg.gammaEq=deg2rad(20.633904875622);x0=[cfg.gammaEq;1];assert(norm(rhs(x0,p,n),inf)<1e-10);
search=struct('tClearMin',.001,'tClearMax',1,'cctTolerance',1e-3,'coarseStep',.05,'extendedEnd',16,'tightRelTol',1e-11,'tightAbsTol',[1e-13 1e-14]);
assert(cfg.tEnd==8&&cfg.RelTol==1e-9&&isequal(cfg.AbsTol,[1e-11 1e-12])&&cfg.faultMaxStep==.002&&cfg.postMaxStep==.005);
Uret=[.90 .80 .70 .65 .60 .55 .50 .45 .40 .35 .30 .25];
coarse=unique([search.tClearMin:.05:search.tClearMax,search.tClearMax,.1,.25]);
criteria="Full valid integration; all stored postfault samples at t>=0.9*tEnd have abs(gamma-gammaEq)<0.1 deg and abs(nu-1)<1e-4; final abs(Pe-1)<1e-4. Absolute tEnd=8 s; endpoint horizon check tEnd=16 s. Unwrapped gamma, no modulo.";
trials=table();results=table();baselineCheck=table();endpointRuns=cell(12,1);
fprintf('Saved baseline clearing time = %.6f s. New results: %s\n',base.cfg.tClear,out);
for j=1:numel(Uret)
 u=Uret(j);pf=p;pf.U=p.U*u;nf=rebuild(0,pf,n);
 [b,~]=assess(base.cfg.tClear,cfg,false,x0,p,n,pf,nf);old=base.summary(abs(base.summary.Uret-u)<1e-10,:);
 expected=old.classification;if expected=="NON_RECOVERED_WITHIN_HORIZON",expected="NOT_RECOVERED";end
 delta=max(abs([b.gamma_clear_deg-old.gamma_clear_deg,b.nu_clear-old.nu_clear]));
 assert(b.classification==expected&&delta<1e-6,'Baseline regression failed');
 baselineCheck=[baselineCheck;table(u,base.cfg.tClear,delta,b.classification,'VariableNames',{'Uret','baseline_tClear_s','max_clearing_error','classification'})];
 fprintf('Uret %.2f: baseline verified; coarse scan (%d times)\n',u,numel(coarse));
 rows=table();
 for k=1:numel(coarse)
  [q,~]=assess(coarse(k),cfg,false,x0,p,n,pf,nf);q.Uret=u;q.phase="coarse";rows=[rows;struct2table(q)];
 end
 trials=[trials;rows];classes=rows.classification;isR=classes=="RECOVERED";isN=classes=="NOT_RECOVERED";known=isR|isN;
 changes=find(classes(1:end-1)~=classes(2:end));lo=NaN;hi=NaN;cct=NaN;status="UNRESOLVED";detail="";verified=false;
 if ~all(known)
  status="UNRESOLVED_COARSE";detail="Validity stops or numerical failures prevent a single assessed transition.";
  if any(isR),lo=max(rows.tClear_s(isR));end
  cand=rows.tClear_s(isN & rows.tClear_s>lo);if ~isempty(cand),hi=min(cand);end
 elseif all(isR)
  lo=search.tClearMax;status="NO_TRANSITION_UP_TO_MAX";detail="All sampled clearing times recovered; no upper CCT bound.";
 elseif all(isN)
  hi=search.tClearMin;status="NO_RECOVERY_IN_SAMPLED_RANGE";detail="No recovered lower endpoint; finite-horizon result only.";
 elseif numel(changes)~=1||~isR(1)||~isN(end)
  status="MULTIPLE_OR_REVERSED_TRANSITIONS";detail="Sampled recovery is not a single recovered-to-not-recovered transition; no bisection.";
 else
  k=changes(1);lo=rows.tClear_s(k);hi=rows.tClear_s(k+1);status="BRACKETED";
  while hi-lo>search.cctTolerance
   mid=(lo+hi)/2;[q,~]=assess(mid,cfg,false,x0,p,n,pf,nf);q.Uret=u;q.phase="bisection";trials=[trials;struct2table(q)];
   if q.classification=="RECOVERED",lo=mid;elseif q.classification=="NOT_RECOVERED",hi=mid;else,status="UNRESOLVED_BISECTION";detail="Ambiguous midpoint; preserved last assessed bracket.";break;end
  end
  if status=="BRACKETED"
   extended=cfg;extended.tEnd=search.extendedEnd;ok=true;storeRuns=cell(2,2);
   for e=1:2
    tc=[lo hi];want=["RECOVERED" "NOT_RECOVERED"];
    for h=1:2
     testcfg=cfg;if h==2,testcfg=extended;end
     [q,runData]=assess(tc(e),testcfg,true,x0,p,n,pf,nf);q.Uret=u;q.phase="endpoint_tight_"+string(testcfg.tEnd)+"s";trials=[trials;struct2table(q)];storeRuns{e,h}=runData;ok=ok&&(q.classification==want(e));
    end
   end
   endpointRuns{j}=storeRuns;verified=ok;
   if ok,cct=(lo+hi)/2;status="NUMERICAL_CCT";detail="Single transition on coarse grid; <=1 ms bracket; tight 8 s and 16 s endpoint classifications agree. Not a proof of asymptotic stability/instability or global monotonicity.";
   else,status="UNRESOLVED_ENDPOINT_CHECK";detail="Tighter integration or extended horizon changed endpoint classification; no CCT estimate reported.";end
  end
 end
 results=[results;table(u,100*(1-u),cct,status,lo,hi,verified,string(detail),'VariableNames',{'Uret','Sag_Depth_percent','CCT_s','Status','CCT_lower_s','CCT_upper_s','Endpoints_verified','Detail'})];
 writetable(results,fullfile(out,'cct_results.csv'));writetable(trials,fullfile(out,'cct_trials.csv'));writetable(baselineCheck,fullfile(out,'baseline_regression.csv'));
 save(fullfile(out,'cct_results.mat'),'results','trials','baselineCheck','endpointRuns','p','n','cfg','search','x0','Uret','coarse','criteria','-v7.3');
 fprintf('Uret %.2f: %s, estimate %.6f, bounds [%.6f, %.6f]\n',u,status,cct,lo,hi);
end
assert(height(results)==12&&numel(unique(results.Uret))==12);
disp(results(:,{'Uret','Sag_Depth_percent','CCT_s','Status','CCT_lower_s','CCT_upper_s'}));
fprintf('%s\nNOT_RECOVERED is finite-horizon non-recovery, not proven instability.\nSaved to %s\n',criteria,out);

%% Visualization of the already calculated and saved CCT results.
plot_cct_results(results,out);

function [q,r]=assess(tc,cfg,tighter,x0,p,n,pf,nf)
q=struct('tClear_s',tc,'horizon_s',cfg.tEnd,'tighter',tighter,'classification',"UNRESOLVED",'reason',"",'gamma_clear_deg',NaN,'nu_clear',NaN,'final_gamma_deg',NaN,'final_nu',NaN,'final_power',NaN,'max_voltage_pu',NaN,'max_limited_current_pu',NaN,'tail_gamma_error_deg',NaN,'tail_nu_error',NaN);
r=struct();rel=cfg.RelTol;ab=cfg.AbsTol;mf=cfg.faultMaxStep;mp=cfg.postMaxStep;
if tighter,rel=1e-11;ab=[1e-13 1e-14];mf=mf/2;mp=mp/2;end
try
 sf=segment([0 tc],x0,pf,nf,rel,ab,mf);r.fault=sf;
 [state,why]=segment_status(sf,tc,pf,nf);if state~="OK",q.classification=state;q.reason="fault: "+why;return;end
 xc=sf.x(end,:)';q.gamma_clear_deg=rad2deg(xc(1));q.nu_clear=xc(2);
 sp=segment([tc cfg.tEnd],xc,p,n,rel,ab,mp);r.postfault=sp;
 [state,why]=segment_status(sp,cfg.tEnd,p,n);if state~="OK",q.classification=state;q.reason="postfault: "+why;return;end
 df=coupled(sf.x(:,1),1,pf,nf,0);dp=coupled(sp.x(:,1),1,p,n,0);
 q.max_voltage_pu=max([abs(df.u);abs(dp.u)]);q.max_limited_current_pu=max([abs(df.ic);abs(dp.ic)]);
 assert(q.max_limited_current_pu<=p.Imax+1e-9,'Current limiter bound violated');
 tail=sp.t>=.9*cfg.tEnd;q.tail_gamma_error_deg=max(abs(rad2deg(sp.x(tail,1)-cfg.gammaEq)));q.tail_nu_error=max(abs(sp.x(tail,2)-1));
 q.final_gamma_deg=rad2deg(sp.x(end,1));q.final_nu=sp.x(end,2);q.final_power=dp.P(end);
 rec=any(tail)&&q.tail_gamma_error_deg<.1&&q.tail_nu_error<1e-4&&abs(q.final_power-1)<1e-4;
 if rec,q.classification="RECOVERED";q.reason="Existing recovery tolerances satisfied";else,q.classification="NOT_RECOVERED";q.reason="Valid full run, recovery test failed within finite horizon; not proof of instability";end
catch err
 q.classification="UNRESOLVED";q.reason=string(err.identifier)+": "+string(err.message);
end
end
function [state,why]=segment_status(s,target,p,n)
state="OK";why="";d=coupled(s.x(:,1),1,p,n,0);
if any(abs(d.u)>p.Vmax+1e-9)||any(s.x(:,2)<=.5)||any(ismember(s.sol.ie,[2 3]))
 state="NOT_ASSESSED_VALIDITY_STOP";why="Converter voltage limit or nu<=0.5 event";
elseif abs(s.t(end)-target)>1e-9
 if numel(s.t)==1 && (abs(d.u)>p.Vmax||s.x(1,2)<=.5),state="NOT_ASSESSED_VALIDITY_STOP";why="Invalid segment initial state";
 else,state="UNRESOLVED";why="Integrator ended before requested endpoint without validity event";end
elseif any(~isfinite(s.x),'all'),state="UNRESOLVED";why="Nonfinite state";
end
end
function seg=segment(span,x0,p,n,rel,ab,maxstep)
 d=coupled(x0(1),1,p,n,0);
 if abs(d.u)>p.Vmax||x0(2)<=.5
  seg=struct('t',span(1),'x',x0','sol',struct('xe',[],'ye',zeros(2,0),'ie',[]),'activeDuration',0);return;
 end
 opt=odeset('RelTol',rel,'AbsTol',ab,'MaxStep',maxstep,'Events',@(t,x)events(t,x,p,n));
 sol=ode45(@(t,x)rhs(x,p,n),span,x0,opt);
 count=6001;if span(1)==0,count=501;end
 t=unique(sort([linspace(span(1),sol.x(end),count) sol.x sol.xe]));x=deval(sol,t)';
 cut=unique(sort([span(1),sol.xe(sol.ie==1),sol.x(end)]));duration=0;
 for j=1:numel(cut)-1,mid=mean(cut(j:j+1));xx=deval(sol,mid);dd=coupled(xx(1),1,p,n,0);if dd.rawI>1.2,duration=duration+diff(cut(j:j+1));end;end
 seg=struct('t',t','x',x,'sol',sol,'activeDuration',duration);
end
function f=rhs(x,p,n)
 d=coupled(x(1),1,p,n,0);f=rhs_power(x,p,d.P);
end
function f=rhs_power(x,p,P)
 f=[p.omega*(x(2)-1);((1-P)/x(2)-p.D*(x(2)-1))/p.mass];
end
function [v,term,dir]=events(~,x,p,n)
 d=coupled(x(1),1,p,n,0);v=[d.rawI-p.Imax;p.Vmax-abs(d.u);x(2)-.5];term=[0;1;1];dir=[0;-1;-1];
end
function P=electrical_power(g,p,n)
 d=coupled(g,1,p,n,0);P=d.P;
end
function s=switchvalue(g,p,n)
 d=coupled(g,1,p,n,1);s=d.rawI-1.2;
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















function plot_cct_results(results,out)
% This function reads the table only; it never changes or saves simulation results.
assert(height(results)==12,'Expected all 12 voltage sag scenarios');
R=sortrows(results,'Sag_Depth_percent');status=string(R.Status);
verified=status=="NUMERICAL_CCT";lower=status=="NO_TRANSITION_UP_TO_MAX";
unresolved=status=="UNRESOLVED_ENDPOINT_CHECK";
assert(all(verified|lower|unresolved),'Unrecognized status: review plotting semantics before export');
assert(all(isfinite(R.CCT_s(verified))) && all(isfinite(R.CCT_lower_s(lower))));
assert(all(isfinite(R.CCT_lower_s(unresolved)))&&all(isfinite(R.CCT_upper_s(unresolved))));
assert(all(R.CCT_upper_s(unresolved)>=R.CCT_lower_s(unresolved)));
paths=fullfile(out,{'CCT_vs_Voltage_Sag.png','CCT_vs_Voltage_Sag.fig','CCT_vs_Voltage_Sag.pdf'});
assert(~any(cellfun(@isfile,paths)),'CCT figures already exist: no existing figure will be overwritten.');
fig=figure('Color','w','Visible','off','Units','pixels','Position',[80 80 1200 760],'Renderer','painters');
ax=axes(fig,'Position',[.10 .13 .86 .78]);hold(ax,'on');
set(ax,'FontName','Arial','FontSize',13,'LineWidth',1,'TickDir','out','Box','off','GridAlpha',.12);
x=R.Sag_Depth_percent;y=R.CCT_s;y(~verified)=NaN;
% NaNs break the line at unresolved/lower-bound scenarios.
h1=plot(ax,x,y,'-o','Color',[.06 .30 .57],'MarkerFaceColor',[.06 .30 .57],'MarkerSize',7,'LineWidth',1.8,'DisplayName','Numerical CCT (endpoint verified)');
% FCT reference line removed at user request.
h3=plot(ax,x(lower),R.CCT_lower_s(lower),'^','LineStyle','none','Color',[.08 .48 .35],'MarkerFaceColor','w','MarkerSize',9,'LineWidth',1.7,'DisplayName','Lower-bound observation');
ytop=max([R.CCT_s(verified);R.CCT_lower_s(lower);R.CCT_upper_s(unresolved);.25]);
for k=find(lower)'
 quiver(ax,x(k),R.CCT_lower_s(k),0,.055*ytop,0,'Color',[.08 .48 .35],'LineWidth',1.3,'MaxHeadSize',.8,'HandleVisibility','off');
end
for bound=unique(R.CCT_lower_s(lower))'
 ix=lower & R.CCT_lower_s==bound;
 text(ax,mean(x(ix)),bound+.085*ytop,sprintf('No transition observed up to %g s',bound),'HorizontalAlignment','center','FontSize',11,'Color',[.08 .38 .28]);
end
% Marker sits on the saved LOWER bracket endpoint, never on an invented CCT midpoint.
h4=plot(ax,x(unresolved),R.CCT_lower_s(unresolved),'d','LineStyle','none','Color',[.78 .30 .08],'MarkerFaceColor','w','MarkerSize',9,'LineWidth',1.7,'DisplayName','Unresolved (8 s test bracket)');
for k=find(unresolved)'
 lo=R.CCT_lower_s(k);hi=R.CCT_upper_s(k);
 plot(ax,[x(k) x(k)],[lo hi],'-','Color',[.78 .30 .08],'LineWidth',1.7,'HandleVisibility','off');
 plot(ax,x(k)+[-.6 .6],[hi hi],'-','Color',[.78 .30 .08],'LineWidth',1.4,'HandleVisibility','off');
 text(ax,x(k),lo+.035*ytop,'Unresolved','HorizontalAlignment','center','FontSize',11,'Color',[.65 .22 .06]);
end
xlabel(ax,'Voltage Sag Depth (%)');ylabel(ax,'Critical Clearing Time (s)');
title(ax,'Critical Clearing Time vs. Voltage Sag Severity','FontSize',17,'FontWeight','normal');
xlim(ax,[min(x)-5 max(x)+5]);ylim(ax,[0 1.18*ytop]);xticks(ax,10:10:80);grid(ax,'on');
legendHandles=h1;if any(lower),legendHandles(end+1)=h3;end;if any(unresolved),legendHandles(end+1)=h4;end;legend(ax,legendHandles,'Location','northeast','Box','off','FontSize',11);
fig.UserData=struct('verifiedCount',nnz(verified),'lowerBoundCount',nnz(lower),'unresolvedCount',nnz(unresolved),'sourceResults',R,'unresolvedMarker','saved lower bracket endpoint; not a CCT estimate');
% MATLAB PNG/PDF export through an ASCII directory avoids Unicode path issues.
[~,id]=fileparts(tempname);tmp=fullfile(getenv('SystemRoot'),'Temp',['cct_' id]);mkdir(tmp);
cleanup=onCleanup(@()rmdir(tmp,'s'));
exportgraphics(fig,fullfile(tmp,'plot.png'),'Resolution',300,'BackgroundColor','white');
exportgraphics(fig,fullfile(tmp,'plot.pdf'),'ContentType','vector','BackgroundColor','white');
savefig(fig,fullfile(tmp,'plot.fig'));
copyfile(fullfile(tmp,'plot.png'),paths{1});copyfile(fullfile(tmp,'plot.fig'),paths{2});copyfile(fullfile(tmp,'plot.pdf'),paths{3});
close(fig);
fprintf('CCT plot saved: %d numerical, %d lower-bound, %d unresolved cases. No CCT recomputation.\n%s\n',nnz(verified),nnz(lower),nnz(unresolved),out);
end

function refine_endpoint_case(source)
% Targeted extended-horizon audit. All original model/solver functions are reused.
a=load(fullfile(source,'cct_results.mat'));k=find(a.results.Status=="UNRESOLVED_ENDPOINT_CHECK");assert(numel(k)==1);
u=a.results.Uret(k);p=a.p;n=a.n;x0=a.x0;pf=p;pf.U=p.U*u;nf=rebuild(0,pf,n);
out=fullfile(fileparts(source),['results_refined_' char(datetime('now','Format','yyyyMMdd_HHmmss_SSS'))]);mkdir(out);
qcfg=a.cfg;qcfg.tEnd=32;qcfg.RelTol=1e-11;qcfg.AbsTol=[1e-13 1e-14];qcfg.faultMaxStep=.001;qcfg.postMaxStep=.0025;
oldBracket=[a.results.CCT_lower_s(k),a.results.CCT_upper_s(k)];
% Include original bracket and neighboring known baseline trial times.
old=a.trials(a.trials.Uret==u & a.trials.phase=="coarse",:);
left=max(old.tClear_s(old.tClear_s<oldBracket(1)));right=min(old.tClear_s(old.tClear_s>oldBracket(2)));
times=unique([linspace(left,right,9),oldBracket]);audit=table();saved=cell(0,1);
for t=times
 [q,~]=assess(t,qcfg,false,x0,p,n,pf,nf);q.Uret=u;q.phase="extended_coarse";audit=[audit;struct2table(q)];fprintf('32 s: clear %.9f -> %s\n',t,q.classification);
end
r=audit.classification=="RECOVERED";nr=audit.classification=="NOT_RECOVERED";
assert(all(r|nr)&&r(1)&&nr(end)&&nnz(diff(r))==1,'Extended scan is ambiguous; cannot refine a single transition');
c=find(diff(r),1);lo=audit.tClear_s(c);hi=audit.tClear_s(c+1);
while hi-lo>2.5e-4
 mid=(lo+hi)/2;[q,~]=assess(mid,qcfg,false,x0,p,n,pf,nf);q.Uret=u;q.phase="extended_bisection";audit=[audit;struct2table(q)];
 assert(any(q.classification==["RECOVERED","NOT_RECOVERED"]),'Ambiguous midpoint');
 if q.classification=="RECOVERED",lo=mid;else,hi=mid;end
end
ok=true;
for tc=[lo hi]
 want="RECOVERED";if tc==hi,want="NOT_RECOVERED";end
 for horizon=[32 64]
  for precision=1:2
   test=qcfg;test.tEnd=horizon;
   if precision==2,test.RelTol=1e-12;test.AbsTol=[1e-14 1e-15];test.faultMaxStep=.0005;test.postMaxStep=.00125;end
   [q,runData]=assess(tc,test,false,x0,p,n,pf,nf);q.Uret=u;q.phase="endpoint_"+string(horizon)+"s_precision"+string(precision);audit=[audit;struct2table(q)];saved{end+1}=runData;
   ok=ok&&(q.classification==want);fprintf('%s: %.9f -> %s\n',q.phase,tc,q.classification);
  end
 end
end
originalResults=a.results;results=originalResults;
results.CCT_lower_s(k)=lo;results.CCT_upper_s(k)=hi;results.Endpoints_verified(k)=ok;
if ok
 results.CCT_s(k)=(lo+hi)/2;results.Status(k)="NUMERICAL_CCT";
 results.Detail(k)="Extended 32 s local scan: single sampled transition. Bracket <=0.25 ms. Endpoints agree at 32/64 s with RelTol 1e-11 and 1e-12. Finite-horizon numerical threshold, not a mathematical guarantee.";
else
 results.CCT_s(k)=NaN;results.Detail(k)="Extended 32/64 s endpoint verification failed; see refinement audit.";
end
results.Assessment_horizon_s=repmat(a.cfg.tEnd,height(results),1);results.Assessment_horizon_s(k)=32;
results.Verification_horizon_s=repmat(a.search.extendedEnd,height(results),1);results.Verification_horizon_s(k)=64;
a.results=results;a.originalResults=originalResults;a.refinementAudit=audit;a.refinementEndpointRuns=saved;a.refinementConfig=qcfg;a.refinementSource=source;
save(fullfile(out,'cct_results.mat'),'-struct','a','-v7.3');writetable(results,fullfile(out,'cct_results.csv'));writetable(audit,fullfile(out,'refinement_trials.csv'));
plot_cct_results(results,out);disp(results);fprintf('Refined results: %s\n',out);
end


