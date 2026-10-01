%% Fixed-Pstar numerical recovery domain of the unchanged ideal radial benchmark.
% A transient gamma is NOT assigned Pstar=Pe(gamma). Pstar is always 1.
% No contraction certificate, controller redesign, EMT, or metric optimization.
clearvars;clc;close all;format long g
root=fileparts(mfilename('fullpath'));out=fullfile(root,'fixed_Pstar_recovery_SCR5_XR10');if ~isfolder(out),mkdir(out);end
diary(fullfile(out,'run_log.txt'));cleanup=onCleanup(@()diary('off'));
a=load(fullfile(root,'current_limiter_SCR5_XR10','current_limiter_SCR5_XR10_results.mat'),'results');p=a.results.parameters;n=a.results.network;
assert(p.H==3.4&&p.K==20&&p.Kd==5&&p.D==25&&p.mass==6.8&&p.Imax==1.2&&p.Vmax==1.3&&p.AD==.2);
cfg=struct('RelTol',2e-9,'AbsTol',2e-11,'MaxStep',.25,'hold',10,'angleTol',.1,'nuTol',1e-4,'fieldTol',1e-4,'nuFloor',.1,'escapeDeg',360,'eq',NaN);
% Optional refresh reuses saved arrays; it does not rerun the grid sweep.
if strcmp(getenv('GFM_RECOVERY_REFRESH'),'1'),refresh_saved(out,p,n);diary off;return;end
g=(-180:.002:180)';d=net(deg2rad(g),p,n);switches=roots_on_grid(g,abs(d.q)-p.Imax,@(z)abs(net(deg2rad(z),p,n).q)-p.Imax);
% Include stationary points in the root search to detect tangencies as well.
stationary=roots_on_grid(g,d.Pg,@(z)net(deg2rad(z),p,n).Pg);cuts=unique([-180;switches;stationary;180]);eq=[];
for j=1:numel(cuts)-1,gg=linspace(cuts(j),cuts(j+1),1001)';dd=net(deg2rad(gg),p,n);eq=[eq;roots_on_grid(gg,dd.Pe-1,@(z)net(deg2rad(z),p,n).Pe-1)];end %#ok<AGROW>
eq=uniquetol(eq,1e-9,'DataScale',1);ed=net(deg2rad(eq),p,n);lam=complex(zeros(numel(eq),2));status=strings(numel(eq),1);
for j=1:numel(eq),lam(j,:)=eig([0 1;-ed.Pg(j)/p.mass -p.D/p.mass]);if ed.Pg(j)>0,status(j)="LOCALLY_STABLE";else,status(j)="SADDLE";end;end
equilibria=table(eq,ed.active,abs(ed.ic),abs(ed.vc),ed.Pg,lam(:,1),lam(:,2),status,'VariableNames',{'gamma_eq_deg','limiter_active','Iconv_mag','Vc_mag','Pe_gamma','eigenvalue_1','eigenvalue_2','classification'});
cfg.eq=eq(find(~ed.active & ed.Pg>0,1));saddle=eq(find(ed.active & ed.Pg<0 & eq>0,1));assert(~isempty(saddle));
writetable(equilibria,fullfile(out,'fixed_Pstar_equilibria.csv'));disp(equilibria);fprintf('Nominal %.12f deg; switches %s deg\n',cfg.eq,mat2str(switches,13));
% Check the optimized RHS closure against the previously validated electrical map.
prior=load(fullfile(root,'diagnose_radial_limiter_instability_SCR5_XR10','radial_limiter_diagnostics.mat'),'results');check=net(deg2rad(prior.results.scan.gamma_deg),p,n);
closureError=max(abs(check.Pe-prior.results.scan.Pe));assert(closureError<1e-12);
% Primary requested 361 x 51 grid. Each vectorized ODE block uses a shared
% adaptive mesh; every individual switching/recovery event has its own index.
[G,N]=meshgrid(0:.25:90,.95:.002:1.05);initial=[G(:),N(:)];
[scan,episodes]=run_grid(initial,p,n,cfg,100,200);scan.scan_stage=repmat("primary",height(scan),1);
save(fullfile(out,'primary_checkpoint.mat'),'scan','episodes','equilibria','cfg','p','n','-v7.3');
% Adaptive domain extension: test edge rows before adding the next rows.
% The basin separatrix can extend indefinitely in nu; report an open edge,
% rather than falsely claiming a closed/global basin after a finite expansion.
domainHistory=strings(0,1);nuLo=.95;nuHi=1.05;
for step=1:3
 edge=scan.nu0==nuLo | scan.nu0==nuHi;mix=any(scan.classification(edge)=="RECOVERED")&&any(scan.classification(edge)~="RECOVERED");if ~mix,break;end
 span=[.1 .2 .5];newLo=1-span(step);newHi=1+span(step);nv=unique([newLo:.01:nuLo-.005,nuHi+.005:.01:newHi,newLo,newHi]);[gg,nn]=meshgrid(0:.5:90,nv);
 [st,ep]=run_grid([gg(:),nn(:)],p,n,cfg,100,200);st.scan_stage=repmat("adaptive_nu",height(st),1);scan=[scan;st];episodes=[episodes;ep];
 domainHistory(end+1)=sprintf('Boundary intersects nu edges %.3f/%.3f: expanded to %.3f/%.3f.',nuLo,nuHi,newLo,newHi);nuLo=newLo;nuHi=newHi;
end
if any(scan.gamma0_deg==90 & scan.classification=="RECOVERED")
 for hi=[120 150 180]
  lo=max(scan.gamma0_deg);[gg,nn]=meshgrid(lo+.5:.5:hi,unique(scan.nu0));[st,ep]=run_grid([gg(:),nn(:)],p,n,cfg,100,200);st.scan_stage=repmat("adaptive_angle",height(st),1);scan=[scan;st];episodes=[episodes;ep];if ~any(st.gamma0_deg==hi & st.classification=="RECOVERED"),break;end
 end
end
domainHistory(end+1)="The separatrix remains open at the finite nu limits; this is a sampled window, not a globally exhausted basin.";
fprintf('%s\n',domainHistory);
% Fine nu=1 angle slice, much finer than the primary grid.
[slice,sliceEpisodes]=run_grid([(0:.01:90)',ones(9001,1)],p,n,cfg,100,200);
% Refine EVERY recovery/non-recovery bracket on selected nu rows, then add
% midpoint nu rows. Distinguish finite-horizon recovery from the saddle basin.
nuRows=unique([linspace(nuLo,nuHi,101),.95:.002:1.05,1]);
boundary=refine_boundary(nuRows,scan,p,n,cfg,200);writetable(boundary,fullfile(out,'fixed_Pstar_recovery_boundary.csv'));
% A finite settling horizon biases the apparent boundary near a saddle.
% Repeat the same brackets at 2000 analytical time units, and verify a tighter
% recovery ball simultaneously. No finite horizon is called a theorem.
longCfg=cfg;longCfg.MaxStep=.5;
longBoundary=refine_boundary(nuRows,scan,p,n,longCfg,2000);
writetable(longBoundary,fullfile(out,'fixed_Pstar_recovery_boundary_extended.csv'));
one=longBoundary(abs(longBoundary.nu0-1)<1e-12,:);assert(height(one)==1);
% Representatives include categories A/B/C, both sides of the refined bound,
% longest recovered limiter dwell, and farthest recovered penetration.
rec=scan.classification=="RECOVERED";lim=scan.limiter_entry_count>0;ixB=[];candidate=rec & lim & ~scan.initial_limiter_active;if any(candidate),[~,ixB]=max(replace_invalid(scan.total_limiter_active_time,candidate));end
[~,ir]=max(replace_invalid(scan.total_limiter_active_time,rec));[~,ip]=max(replace_invalid(scan.max_gamma_deg,rec));
repIC=[cfg.eq+1,1;20,1;one.gamma_recovered_deg,1;one.gamma_nonrecovered_deg,1;scan.gamma0_deg(ir),scan.nu0(ir);scan.gamma0_deg(ip),scan.nu0(ip)];
repNames=["Inactive recovery";"Initially limited recovery";"Near boundary: recovered";"Near boundary: not recovered";"Longest limited dwell (scan)";"Farthest penetration (scan)"];
if ~isempty(ixB),repIC(end+1,:)=[scan.gamma0_deg(ixB),scan.nu0(ixB)];repNames(end+1)="Inactive start, later limited";end
[~,far]=max(longBoundary.gamma_recovered_deg);repIC(end+1,:)=[longBoundary.gamma_recovered_deg(far),longBoundary.nu0(far)];repNames(end+1)="Farthest extended-boundary recovery";
[~,slow]=max(longBoundary.recovered_limiter_dwell);repIC(end+1,:)=[longBoundary.gamma_recovered_deg(slow),longBoundary.nu0(slow)];repNames(end+1)="Longest extended-boundary dwell";
[repTable,repEpisodes,representatives]=simulate_block(repIC,p,n,longCfg,2000);repTable.name=repNames;
% Solver/tolerance check on the boundary and all representatives.
strict=longCfg;strict.RelTol=2e-11;strict.AbsTol=2e-13;strict.MaxStep=.125;
[strictTable,~,strictRep]=simulate_block(repIC,p,n,strict,2000);
validation=struct('savedClosureMaxError',closureError,'representativeClassificationAgreement',all(repTable.classification==strictTable.classification),'maxRepresentativeFinalAngleDifference_deg',max(abs(repTable.final_gamma_deg-strictTable.final_gamma_deg)),'primaryVsTightBallDifferent',sum(scan.classification~=scan.tight_classification),'extendedBoundaryTightDisagreements',sum(longBoundary.tight_disagreement));
assert(validation.representativeClassificationAgreement,'Refine solver: representative classifications disagree.');
validation.maxRepresentativeDwellDifference=max(abs(repTable.total_limiter_active_time-strictTable.total_limiter_active_time));validation.maxRepresentativeExcursionDifference_deg=max(abs(repTable.max_angle_excursion_deg-strictTable.max_angle_excursion_deg));validation.phaseDepartureThreshold_deg=180;
writetable(scan,fullfile(out,'fixed_Pstar_recovery_scan.csv'));writetable(episodes,fullfile(out,'fixed_Pstar_limiter_episodes.csv'));writetable(slice,fullfile(out,'fixed_Pstar_angle_slice.csv'));writetable(repTable,fullfile(out,'representative_trajectories.csv'));writetable(repEpisodes,fullfile(out,'representative_limiter_episodes.csv'));
figures(scan,boundary,longBoundary,equilibria,switches,representatives,repNames,cfg,p,out);
results=struct('parameters',p,'network',n,'Pstar',1,'settings',cfg,'equilibria',equilibria,'switches_deg',switches,'scan',scan,'episodes',episodes,'slice',slice,'sliceEpisodes',sliceEpisodes,'boundary',boundary,'extendedBoundary',longBoundary,'representativeTable',repTable,'representatives',{representatives},'representativeEpisodes',repEpisodes,'strictRepresentatives',{strictRep},'validation',validation,'domainHistory',domainHistory);
save(fullfile(out,'fixed_Pstar_recovery_results.mat'),'results','-v7.3');report(results,out);disp(validation);fprintf('FIXED-PSTAR RECOVERY STUDY COMPLETE\n');diary off

function d=net(g,p,n)
v=p.E*exp(1i*(g+p.phi));q=(v-p.t*p.U)/p.Z;r=abs(q);active=r>p.Imax;ic=q;ic(active)=p.Imax*q(active)./r(active);vc=p.t*p.U+p.Z*ic;vc(~active)=v(~active);
vg=n.vgE*vc+n.vgU*p.U;ig=n.igE*vc+n.igU*p.U;S=vg.*conj(ig);u=vc+(n.Zconverter+p.AD)*ic;
qg=1i*v/p.Z;icg=qg;icg(active)=p.Imax*(qg(active)./r(active)-q(active).*real(conj(q(active)).*qg(active))./r(active).^3);vcg=p.Z*icg;
Pg=real(n.vgE*vcg.*conj(ig)+vg.*conj(n.igE*vcg));
d=struct('q',q,'ic',ic,'vc',vc,'vg',vg,'ig',ig,'Pe',real(S),'Qe',imag(S),'u',u,'active',active,'Pg',Pg);
end
function dx=rhs(~,x,p,n)
g=x(1:2:end);nu=x(2:2:end);d=net(g,p,n);dx=zeros(size(x));dx(1:2:end)=nu-1;dx(2:2:end)=((1-d.Pe)./nu-p.D*(nu-1))/p.mass;
end
function r=ball(x,p,n,c,tight)
g=x(1:2:end,:);nu=x(2:2:end,:);d=net(g,p,n);acc=((1-d.Pe)./nu-p.D*(nu-1))/p.mass;scale=1;if tight,scale=.1;end
r=max(cat(3,abs(rad2deg(g)-c.eq)/(c.angleTol*scale),abs(nu-1)/(c.nuTol*scale),abs(acc)/(c.fieldTol*scale)),[],3)-1;
end
function [v,term,dir]=events(~,x,p,n,c)
g=x(1:2:end);nu=x(2:2:end);d=net(g,p,n);v=[abs(d.q)-p.Imax;ball(x,p,n,c,false);ball(x,p,n,c,true);c.escapeDeg-abs(rad2deg(g)-c.eq);nu-c.nuFloor];term=zeros(size(v));dir=zeros(size(v));
end
function [tab,ep,full]=simulate_block(ic,p,n,c,T)
m=size(ic,1);x=reshape([deg2rad(ic(:,1))';ic(:,2)'],[],1);opts=odeset('RelTol',c.RelTol,'AbsTol',c.AbsTol,'MaxStep',c.MaxStep,'Events',@(t,x)events(t,x,p,n,c));
sol=ode45(@(t,x)rhs(t,x,p,n),[0 T],x,opts);tab=table();ep=table();full=cell(m,1);
% Evaluate the dense interpolant once for the entire ensemble. Repeating
% deval's interval lookup separately for every trajectory is unnecessary.
commonTimes=unique([0:.1:T,T,sol.x,sol.xe]);commonStates=deval(sol,commonTimes);
for k=1:m
 ev=sol.ie(:);te=sol.xe(:);stopT=T;escape=te(ev==3*m+k);domain=te(ev==4*m+k);if ~isempty(domain),stopT=min(stopT,domain(1));end;if ~isempty(escape),stopT=min(stopT,escape(1));end
 keep=commonTimes<=stopT;ts=commonTimes(keep);xx=commonStates(2*k-1:2*k,keep);d=net(xx(1,:),p,n);dx=[xx(2,:)-1;((1-d.Pe)./xx(2,:)-p.D*(xx(2,:)-1))/p.mass];
 et=te(ev==k & te>1e-8 & te<=stopT);et=uniquetol(et,1e-8,'DataScale',1);bounds=[0;et;stopT];episodes=zeros(0,8);
 for j=1:numel(bounds)-1
  aa=bounds(j);bb=bounds(j+1);mid=deval(sol,(aa+bb)/2,2*k-1);if abs(net(mid,p,n).q)<=p.Imax,continue;end
  xa=deval(sol,aa,2*k-1:2*k);xb=deval(sol,bb,2*k-1:2*k);closed=bb<stopT-1e-7 || (any(abs(et-bb)<1e-7));exit=bb;if ~closed,exit=NaN;end
  episodes(end+1,:)=[aa,rad2deg(xa(1)),xa(2),exit,rad2deg(xb(1)),xb(2),bb-aa,double(closed)]; %#ok<AGROW>
 end
 [rt,confirm]=settled(sol,k,m,stopT,p,n,c,false);[rtt,ctt]=settled(sol,k,m,stopT,p,n,c,true);
 voltage=all(abs(d.u)<=p.Vmax+1e-9);class="NOT_RECOVERED_WITHIN_HORIZON";tightclass=class;
 if isfinite(confirm),class="RECOVERED";end;if isfinite(ctt),tightclass="RECOVERED";end
 % Half-turn departure is a numerical escape criterion, not a theorem.
 % Continue recording episodes until the horizon or full-turn bailout.
 if max(abs(rad2deg(xx(1,:))-c.eq))>=180,class="PHASE_SLIP_OR_ESCAPE";tightclass=class;end
 if ~isempty(escape)&&escape(1)<=stopT+1e-7,class="PHASE_SLIP_OR_ESCAPE";tightclass=class;end
 if ~voltage,class="OUTSIDE_VOLTAGE_VALIDITY";tightclass=class;end
 if ~isempty(domain)&&domain(1)<=stopT+1e-7,class="OUTSIDE_ANALYTICAL_MODEL_DOMAIN";tightclass=class;end
 first=NaN;firstexit=NaN;firstdwell=0;maxdwell=0;total=0;last=NaN;
 if ~isempty(episodes),first=episodes(1,1);firstexit=episodes(1,4);firstdwell=episodes(1,7);maxdwell=max(episodes(:,7));total=sum(episodes(:,7));ex=episodes(isfinite(episodes(:,4)),4);if ~isempty(ex),last=max(ex);end
  etab=array2table(episodes,'VariableNames',{'entry_time','entry_gamma_deg','entry_nu','exit_time','exit_gamma_deg','exit_nu','dwell_duration','episode_closed'});etab.gamma0_deg=repmat(ic(k,1),height(etab),1);etab.nu0=repmat(ic(k,2),height(etab),1);ep=[ep;etab];
 end
 post=NaN;if class=="RECOVERED"&&isfinite(last),post=max(0,confirm-last);end
 row=struct('gamma0_deg',ic(k,1),'nu0',ic(k,2),'initial_limiter_active',d.active(1),'classification',class,'tight_classification',tightclass,'recovery_time',rt,'recovery_confirmed_time',confirm,'tight_recovery_time',rtt,'post_limiter_recovery_time',post,'limiter_entry_count',size(episodes,1),'first_entry_time',first,'first_exit_time',firstexit,'first_active_dwell_time',firstdwell,'max_single_active_dwell_time',maxdwell,'total_limiter_active_time',total,'last_limiter_exit_time',last,'max_gamma_deg',max(rad2deg(xx(1,:))),'min_gamma_deg',min(rad2deg(xx(1,:))),'max_angle_excursion_deg',max(abs(rad2deg(xx(1,:))-c.eq)),'max_nu',max(xx(2,:)),'min_nu',min(xx(2,:)),'max_abs_nu_deviation',max(abs(xx(2,:)-1)),'max_q_mag',max(abs(d.q)),'max_Iconv_mag',max(abs(d.ic)),'min_Vc_mag',min(abs(d.vc)),'max_Vc_mag',max(abs(d.vc)),'max_u_required',max(abs(d.u)),'max_Pe',max(d.Pe),'min_Pe',min(d.Pe),'max_Qe',max(d.Qe),'min_Qe',min(d.Qe),'voltage_valid',voltage,'final_gamma_deg',rad2deg(xx(1,end)),'final_nu',xx(2,end),'integration_horizon',T,'observation_end',stopT);
 tab=[tab;struct2table(row)];if nargout>2,full{k}=struct('t',ts,'x',xx,'dx',dx,'electrical',d,'episodes',episodes,'summary',row);end
end
end
function [start,confirmed]=settled(sol,k,m,T,p,n,c,tight)
offset=m;if tight,offset=2*m;end;ev=sol.xe(sol.ie==offset+k);bounds=unique([0,ev(ev>1e-8 & ev<T),T]);start=NaN;confirmed=NaN;
for j=1:numel(bounds)-1
 if bounds(j+1)-bounds(j)<c.hold-1e-7,continue;end
 xx=deval(sol,(bounds(j)+bounds(j+1))/2,2*k-1:2*k);if ball(xx,p,n,c,tight)<=0,start=bounds(j);confirmed=start+c.hold;return;end
end
end
function [tab,ep]=run_grid(ic,p,n,c,T,extend)
tab=table();ep=table();batch=96;
tic;for b=1:batch:size(ic,1)
 ii=b:min(b+batch-1,size(ic,1));[tt,ee]=simulate_block(ic(ii,:),p,n,c,T);need=tt.classification=="NOT_RECOVERED_WITHIN_HORIZON" | tt.classification=="PHASE_SLIP_OR_ESCAPE" | (tt.classification=="RECOVERED" & tt.tight_classification~="RECOVERED");
 if extend>T&&any(need)
  [xt,xe]=simulate_block(ic(ii(need),:),p,n,c,extend);tt(need,:)=xt;
  if ~isempty(ee),drop=false(height(ee),1);for k=find(need)',drop=drop | (ee.gamma0_deg==ic(ii(k),1)&ee.nu0==ic(ii(k),2));end;ee(drop,:)=[];end;ee=[ee;xe];
 end
 tab=[tab;tt];ep=[ep;ee];if mod(b-1,batch*10)==0,fprintf('Scan %d/%d, recovered %d, elapsed %.1f s\n',min(b+batch-1,size(ic,1)),size(ic,1),sum(tab.classification=="RECOVERED"),toc);end
end
end
function b=refine_boundary(nus,scan,p,n,c,T)
% Initial brackets are verified by actual integrations for each selected row.
lo=zeros(numel(nus),1);hi=repmat(max(scan.gamma0_deg),numel(nus),1);nus=nus(:);
[lt,~]=simulate_block([lo,nus],p,n,c,T);[ht,~]=simulate_block([hi,nus],p,n,c,T);assert(all(lt.classification=="RECOVERED")&&all(ht.classification~="RECOVERED"),'Domain has multiple/missing recovery brackets; inspect before refinement.');
% Verify coarse rows have a single transition; do not assume topology silently.
rows=unique(scan.nu0);for j=1:numel(rows),ss=sortrows(scan(scan.nu0==rows(j),:),'gamma0_deg');rr=ss.classification=="RECOVERED";assert(sum(diff(rr)~=0)==1,'Multiple transitions need separate brackets.');end
while max(hi-lo)>1e-4
 mid=(lo+hi)/2;[mt,~]=simulate_block([mid,nus],p,n,c,T);r=mt.classification=="RECOVERED";lo(r)=mid(r);hi(~r)=mid(~r);
end
[lt,~]=simulate_block([lo,nus],p,n,c,T);[ht,~]=simulate_block([hi,nus],p,n,c,T);
b=table(nus,lo,hi,(lo+hi)/2,hi-lo,repmat(T,numel(nus),1),lt.tight_classification~="RECOVERED",lt.total_limiter_active_time,lt.last_limiter_exit_time,lt.recovery_confirmed_time,ht.classification,'VariableNames',{'nu0','gamma_recovered_deg','gamma_nonrecovered_deg','gamma_critical_mid_deg','bracket_width_deg','horizon','tight_disagreement','recovered_limiter_dwell','recovered_last_exit','recovery_confirmed_time','outer_classification'});
fprintf('Refined %d boundary rows to %.3g deg at horizon %g.\n',height(b),max(b.bracket_width_deg),T);
end
function r=roots_on_grid(g,y,f)
ii=find(y(1:end-1).*y(2:end)<0);r=g(abs(y)<1e-12);for k=ii',r(end+1,1)=fzero(f,g(k+[0 1]),optimset('TolX',1e-12));end;r=uniquetol(r,1e-8,'DataScale',1);
end
function y=replace_invalid(y,valid)
y(~valid)=-Inf;
end
function figures(s,b,lb,eq,sw,reps,names,c,p,out)
colors=[.12 .56 .30;.8 .4 .08;.65 .15 .15;.5 .2 .7;.2 .2 .2];classes=["RECOVERED","NOT_RECOVERED_WITHIN_HORIZON","PHASE_SLIP_OR_ESCAPE","OUTSIDE_VOLTAGE_VALIDITY","OUTSIDE_ANALYTICAL_MODEL_DOMAIN"];
f=figure('Color','w','Position',[60 60 1450 750]);ax=axes(f);hold(ax,'on');
for j=1:5,ix=s.classification==classes(j);h=scatter(ax,s.gamma0_deg(ix),s.nu0(ix),7,colors(j,:),'filled','DisplayName',strrep(classes(j),'_',' '));h.MarkerFaceAlpha=.32;end
plot(ax,b.gamma_critical_mid_deg,b.nu0,'k--','LineWidth',1.8,'DisplayName','Recovery boundary, T=200');plot(ax,lb.gamma_critical_mid_deg,lb.nu0,'k-','LineWidth',2,'DisplayName','Recovery boundary, T=2000');
for j=1:numel(reps),rr=reps{j};plot(ax,rad2deg(rr.x(1,:)),rr.x(2,:),'LineWidth',1.2,'DisplayName',names(j));end
for j=1:height(eq),plot(ax,eq.gamma_eq_deg(j),1,'p','MarkerSize',12,'MarkerFaceColor','y','MarkerEdgeColor','k','DisplayName',sprintf('Equilibrium %.4f deg: %s',eq.gamma_eq_deg(j),eq.classification(j)));end
for j=1:numel(sw),xline(ax,sw(j),':','Switch','HandleVisibility','off');end
% Stable saddle eigenvector orientation only, not a computed global manifold.
ie=find(eq.classification=="SADDLE" & eq.gamma_eq_deg>0,1);[V,L]=eig([0 1;-eq.Pe_gamma(ie)/p.mass -p.D/p.mass]);v=real(V(:,find(real(diag(L))<0,1)));v=v/max(abs(v));rr=linspace(-.04,.04,20);plot(ax,eq.gamma_eq_deg(ie)+rad2deg(v(1)*rr),1+v(2)*rr,'m-','LineWidth',2,'DisplayName','Local stable eigenvector');
xlim(ax,[min(-20,min(s.gamma0_deg)),max(s.gamma0_deg)]);ylim(ax,[min(s.nu0),max(s.nu0)]);xlabel(ax,'Initial / trajectory gamma (deg, unwrapped)');ylabel(ax,'nu (pu)');grid(ax,'on');title(ax,{'Fixed Pstar=1: numerical recovery domain (not contraction certified)','Colored grid: T<=200; solid boundary and representative paths: T=2000'});legend(ax,'Location','eastoutside','FontSize',8);save_plot(f,out,'01_recovery_phase_plane');
for j=1:2
 f=figure('Color','w','Position',[60 60 1150 700]);ax=axes(f);hold(ax,'on');ix=s.classification=="RECOVERED";v=s.total_limiter_active_time;if j==2,v=s.max_angle_excursion_deg;end
 scatter(ax,s.gamma0_deg(ix),s.nu0(ix),12,v(ix),'filled');plot(ax,b.gamma_critical_mid_deg,b.nu0,'k--','LineWidth',1.3);cb=colorbar(ax);if j==1,cb.Label.String='Total limiter-active time (normalized analytical time)';else,cb.Label.String='Maximum unwrapped angle excursion (deg)';end
 xlabel(ax,'Initial gamma (deg)');ylabel(ax,'Initial nu (pu)');grid(ax,'on');title(ax,'RECOVERED points only: fixed Pstar=1, finite simulation horizon');save_plot(f,out,sprintf('0%d_recovered_map',j+1));
end
ids=find(cellfun(@(r)r.summary.classification=="RECOVERED" && r.summary.limiter_entry_count>0,reps));
for jj=1:numel(ids)
 j=ids(jj);r=reps{j};f=figure('Color','w','Position',[60 30 1150 1050]);tl=tiledlayout(f,5,1,'TileSpacing','compact');data={rad2deg(r.x(1,:)),r.x(2,:),abs(r.electrical.q),abs(r.electrical.ic),r.electrical.Pe};labels={'gamma (deg, unwrapped)','nu (pu)','|q| (pu)','|Iconv| (pu)','Pe (pu)'};
 for k=1:5,ax=nexttile(tl);plot(ax,r.t,data{k},'LineWidth',1.4);hold(ax,'on');for z=1:size(r.episodes,1),xline(ax,r.episodes(z,1),'r:','Entry','HandleVisibility','off');if isfinite(r.episodes(z,4)),xline(ax,r.episodes(z,4),'g:','Exit','HandleVisibility','off');end;end
  if ismember(k,[3 4]),yline(ax,1.2,'k--');elseif k==5,yline(ax,1,'k--','Pstar');end;ylabel(ax,labels{k});grid(ax,'on');xlim(ax,[0,min(r.t(end),max(60,r.summary.recovery_confirmed_time+5))]);
 end
 xlabel(ax,'Normalized analytical time');title(tl,names(j)+": limiter entry, exit and recovery");save_plot(f,out,sprintf('04_entry_exit_%02d',j));
end
end
function save_plot(f,out,name)
drawnow;exportgraphics(f,fullfile(out,[name '.png']),'Resolution',160);savefig(f,fullfile(out,[name '.fig']));close(f);
end
function report(r,out)
s=r.scan;rec=s.classification=="RECOVERED";lim=s.limiter_entry_count>0;one=r.extendedBoundary(abs(r.extendedBoundary.nu0-1)<1e-12,:);short=r.boundary(abs(r.boundary.nu0-1)<1e-12,:);[~,ip]=max(replace_invalid(s.max_gamma_deg,rec));[~,it]=max(replace_invalid(s.total_limiter_active_time,rec));
fid=fopen(fullfile(out,'fixed_Pstar_recovery_report.md'),'w','n','UTF-8');guard=onCleanup(@()fclose(fid));
fprintf(fid,'# Fixed-Pstar numerical recovery domain\n\nPstar=1 for EVERY trajectory, including limiter-active transient states. Same saved SCR5/XR10 radial benchmark; H=3.4, K=20, Kd=5, D=25, Imax=1.2. No EMT, metric optimization, or physical model changes. This is not a contraction certificate.\n\n');
fprintf(fid,'## Fixed-power equilibria\n\n| gamma deg | limiter active | Iconv pu | Vc pu | Pe_gamma pu/rad | eigenvalues | classification |\n|---|---|---|---|---|---|---|\n');
for j=1:height(r.equilibria),e=r.equilibria(j,:);fprintf(fid,'| %.12f | %d | %.9f | %.9f | %.9f | %s, %s | %s |\n',e.gamma_eq_deg,e.limiter_active,e.Iconv_mag,e.Vc_mag,e.Pe_gamma,num2str(e.eigenvalue_1,12),num2str(e.eigenvalue_2,12),e.classification);end
fprintf(fid,'\nNominal attractor: %.12f deg. Current switches: %s deg. Roots searched over [-180,180], separately across nonsmooth branches, including stationary-point cuts. Periodic copies are not nominal recovery.\n\n',r.settings.eq,mat2str(r.switches_deg,13));
fprintf(fid,'## Recovery definition and validity\n\nAngle error <=0.1 deg, |nu-1|<=1e-4, |gamma_dot|<=1e-4 and |nu_dot|<=1e-4 continuously for 10 normalized analytical time units. Recovery entry/exit uses adaptive ODE events on the maximum normalized residual; the tighter check divides all tolerances by 10. recovery_time is the start of the qualifying interval; recovery_confirmed_time is its end. Post-limiter recovery time ends at confirmation.\n\nEach trajectory first runs to 100, with unresolved/tighter-ball and branch-departure cases extended to 200. Near-saddle boundary calculations additionally use 2000 because slow saddle departure can dominate finite-horizon classification. PHASE_SLIP_OR_ESCAPE means an observed |gamma-gamma_eq|>=180 deg, a numerical departure from the nominal angle branch. Data/episodes continue through the horizon to show the subsequent branch; a separate full-turn event at 360 deg is an observation bailout. nu<=0.1 leaves the analytical domain. These are stopping/observation criteria, not theorems. Voltage-invalid trajectories are reported separately, never automatically called unstable. No voltage limiter is added.\n\n');
fprintf(fid,'The vectorized ode45 mesh uses RelTol=2e-9, AbsTol=2e-11, MaxStep=0.25 (0.5 for extended studies), independent event functions per trajectory, no state reset/hysteresis. Saved electrical maxima are evaluated on the adaptive mesh plus dense 0.1-time samples and event times; they are numerical extrema, not certified global bounds. Full representative records include phasors, vector fields and episodes. Open episodes are right-censored: exit time NaN, dwell measured through observation end.\n\n');
fprintf(fid,'## Domain and fractions\n\nPrimary scan: gamma=0:0.25:90, nu=0.95:0.002:1.05. Adaptive extensions and the 0.01 deg angle slice are separate.\n\n');
for msg=r.domainHistory',fprintf(fid,'- %s\n',msg);end
for name=["primary","all sampled rows"],ix=true(height(s),1);if name=="primary",ix=s.scan_stage=="primary";end;fprintf(fid,'\n%s: %d states; recovered %.6f%%; limiter-active episode %.6f%%; limited AND recovered %.6f%%. These are point-count fractions, not phase-space area estimates.\n',name,sum(ix),100*mean(rec(ix)),100*mean(lim(ix)),100*mean(rec(ix)&lim(ix)));end
fprintf(fid,'\n| Classification | Count |\n|---|---:|\n');classes=unique(s.classification);for k=1:numel(classes),fprintf(fid,'| %s | %d |\n',classes(k),sum(s.classification==classes(k)));end
fprintf(fid,'\nRecovered categories: inactive/never limited %d; inactive/later limited %d; initially limited %d.\n\n',sum(rec&~lim),sum(rec&lim&~s.initial_limiter_active),sum(rec&s.initial_limiter_active));
fprintf(fid,'## Angle-jump and boundary results\n\nAt nu0=1, T=200 recovery bracket: [%.9f, %.9f] deg. Extended T=2000 bracket: [%.9f, %.9f] deg; Delta gamma bracket [%.9f, %.9f] deg. The lower endpoint recovers; the upper does not satisfy the finite-horizon recovery definition. Brackets are refined below 1e-4 deg; this bracket resolution is not a claim of an exact global stability boundary.\n\n',short.gamma_recovered_deg,short.gamma_nonrecovered_deg,one.gamma_recovered_deg,one.gamma_nonrecovered_deg,one.gamma_recovered_deg-r.settings.eq,one.gamma_nonrecovered_deg-r.settings.eq);
ie=find(r.equilibria.gamma_eq_deg>0 & r.equilibria.classification=="SADDLE",1);sa=r.equilibria.gamma_eq_deg(ie);fprintf(fid,'Positive limiter-active saddle: %.12f deg. Extended nu0=1 bracket midpoint minus saddle = %.9g deg. The 2-D refined boundary tilts with initial nu; the saddle local stable eigenvector is overlaid. This supports its organizing role but no global stable manifold has been rigorously computed.\n\n',sa,one.gamma_critical_mid_deg-sa);
rr=r.representativeTable(3,:);fprintf(fid,'Recovered-side extended near-boundary trajectory: initial gamma %.9f, limiter entry at t=%.9f and gamma=%.9f; max gamma %.9f; total limiter dwell %.9f; last exit %.9f; recovery confirmation %.9f. All times are normalized analytical time.\n\n',rr.gamma0_deg,rr.first_entry_time,rr.gamma0_deg,rr.max_gamma_deg,rr.total_limiter_active_time,rr.last_limiter_exit_time,rr.recovery_confirmed_time);
fprintf(fid,'## Penetration and dwell\n\nAmong the finite-horizon scan recovered cases: farthest positive gamma %.9f deg, from (gamma0,nu0)=(%.9f,%.9f); penetration beyond positive switch %.9f deg. Largest absolute angle excursion %.9f deg. Longest total limiter dwell %.9f, from (%.9f,%.9f); longest individual episode %.9f. These are sampled maxima, not universal allowable dwell times. Near-boundary extended trajectories can dwell much longer.\n\n',s.max_gamma_deg(ip),s.gamma0_deg(ip),s.nu0(ip),s.max_gamma_deg(ip)-max(r.switches_deg),max(s.max_angle_excursion_deg(rec)),s.total_limiter_active_time(it),s.gamma0_deg(it),s.nu0(it),max(s.max_single_active_dwell_time(rec)));
fprintf(fid,'## Recovery mechanism\n\nBetween the positive switch and the positive saddle, Pe>1 even though Pe_gamma<0. At nu=1 this gives nu_dot<0, followed by nu<1 and gamma_dot<0: the state can exit current limitation and return to nominal. For initially nonunit nu, damping also matters: the sign of 1-Pe ALONE is not the sign of nu_dot. Above the saddle Pe<1 creates positive accelerating tendency at nu=1.\n\nA transient at gamma=20 deg is not a new equilibrium with Pstar=Pe(20). Local instability of the Pe(gamma0)-parameterized equilibrium family does not prohibit recovery of transient states at fixed Pstar=1. The empirical recovery boundary separates outcomes at a stated horizon; it is not a contraction-certified boundary.\n\n');
fprintf(fid,'## Numerical checks\n\n');fn=fieldnames(r.validation);for k=1:numel(fn),fprintf(fid,'- %s: %.12g\n',fn{k},r.validation.(fn{k}));end
fprintf(fid,'\nBoundary endpoints also carry tight-ball disagreement flags. Any flagged finite-horizon boundary depends on the settling tolerance and must not be presented as a physical separatrix. Representative classifications were repeated with RelTol=2e-11, AbsTol=2e-13 and MaxStep=0.125.\n');
fprintf(fid,'\n## Additional physical and branch checks\n\nMaximum voltage command across the scan = %.12f pu < 1.3 pu; voltage-invalid count = %d; analytical-domain-invalid count = %d.\n\n',max(s.max_u_required),sum(~s.voltage_valid),sum(s.classification=="OUTSIDE_ANALYTICAL_MODEL_DOMAIN"));
fprintf(fid,'Largest tested recovered-side initial angle among the extended boundary brackets = %.9f deg. Longest tested recovered-side boundary dwell = %.9f normalized analytical time. These extend the T<=200 scan results; neither is a universal bound.\n\n',max(r.extendedBoundary.gamma_recovered_deg),max(r.extendedBoundary.recovered_limiter_dwell));
fprintf(fid,'| Representative | Initial gamma / nu | Pe initially | nu_dot initially | Minimum nu | Final unwrapped gamma | Classification |\n|---|---|---|---|---|---|---|\n');
for j=1:numel(r.representatives),rr=r.representatives{j};ss=rr.summary;fprintf(fid,'| %s | %.9f / %.9f | %.9f | %.9f | %.9f | %.9f | %s |\n',r.representativeTable.name(j),ss.gamma0_deg,ss.nu0,rr.electrical.Pe(1),rr.dx(2,1),ss.min_nu,ss.final_gamma_deg,ss.classification);end
fprintf(fid,'\nThe outer nu0=1 example approaches the next angular copy near 373.810758941 deg, not the original 13.810758941 deg equilibrium. Its branch departure is classified separately even though the reduced model later settles in that other branch.\n');
fprintf('Recovered %d/%d; limited and recovered %d; nu0=1 extended bracket %.9f .. %.9f deg\n',sum(rec),height(s),sum(rec&lim),one.gamma_recovered_deg,one.gamma_nonrecovered_deg);
end
function t=branch_labels(t)
ix=t.max_angle_excursion_deg>=180 & t.voltage_valid & t.classification~="OUTSIDE_ANALYTICAL_MODEL_DOMAIN";
t.classification(ix)="PHASE_SLIP_OR_ESCAPE";t.tight_classification(ix)="PHASE_SLIP_OR_ESCAPE";
end
function refresh_saved(out,p,n)
a=load(fullfile(out,'fixed_Pstar_recovery_results.mat'),'results');r=a.results;
assert(isequaln(p,r.parameters)&&isequaln(n,r.network),'Saved data must match the current verified network.');
r.scan=branch_labels(r.scan);r.slice=branch_labels(r.slice);
c=r.settings;c.MaxStep=.5;
% Reclassify only the extended outer boundary; no new sweep or root refinement.
[outer,~]=simulate_block([r.extendedBoundary.gamma_nonrecovered_deg,r.extendedBoundary.nu0],p,n,c,2000);r.extendedBoundary.outer_classification=outer.classification;
rec=r.scan.classification=="RECOVERED";candidate=rec & r.scan.limiter_entry_count>0 & ~r.scan.initial_limiter_active;[~,ib]=max(replace_invalid(r.scan.total_limiter_active_time,candidate));
ics=[r.representativeTable.gamma0_deg,r.representativeTable.nu0];names=r.representativeTable.name;ics(7,:)=[r.scan.gamma0_deg(ib),r.scan.nu0(ib)];
if size(ics,1)<8,[~,j]=max(r.extendedBoundary.gamma_recovered_deg);ics(8,:)=[r.extendedBoundary.gamma_recovered_deg(j),r.extendedBoundary.nu0(j)];names(8)="Farthest extended-boundary recovery";end
if size(ics,1)<9,[~,j]=max(r.extendedBoundary.recovered_limiter_dwell);ics(9,:)=[r.extendedBoundary.gamma_recovered_deg(j),r.extendedBoundary.nu0(j)];names(9)="Longest extended-boundary dwell";end
[rt,re,rr]=simulate_block(ics,p,n,c,2000);rt.name=names;
strict=c;strict.RelTol=2e-11;strict.AbsTol=2e-13;strict.MaxStep=.125;[st,~,sr]=simulate_block(ics,p,n,strict,2000);
r.representativeTable=rt;r.representativeEpisodes=re;r.representatives=rr;r.strictRepresentatives=sr;
r.validation.representativeClassificationAgreement=all(rt.classification==st.classification);r.validation.maxRepresentativeFinalAngleDifference_deg=max(abs(rt.final_gamma_deg-st.final_gamma_deg));assert(r.validation.representativeClassificationAgreement);
% Compare dwell and extrema too, not only final equilibria.
r.validation.maxRepresentativeDwellDifference=max(abs(rt.total_limiter_active_time-st.total_limiter_active_time));r.validation.maxRepresentativeExcursionDifference_deg=max(abs(rt.max_angle_excursion_deg-st.max_angle_excursion_deg));
r.validation.phaseDepartureThreshold_deg=180;
writetable(r.scan,fullfile(out,'fixed_Pstar_recovery_scan.csv'));writetable(r.slice,fullfile(out,'fixed_Pstar_angle_slice.csv'));writetable(r.extendedBoundary,fullfile(out,'fixed_Pstar_recovery_boundary_extended.csv'));writetable(rt,fullfile(out,'representative_trajectories.csv'));writetable(re,fullfile(out,'representative_limiter_episodes.csv'));
figures(r.scan,r.boundary,r.extendedBoundary,r.equilibria,r.switches_deg,rr,names,r.settings,p,out);
results=r;save(fullfile(out,'fixed_Pstar_recovery_results.mat'),'results','-v7.3');report(r,out);disp(r.validation);fprintf('SAVED RECOVERY OUTPUTS REFRESHED AND VERIFIED\n');
end
