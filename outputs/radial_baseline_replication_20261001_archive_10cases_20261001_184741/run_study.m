%% Baseline figure replication: five initial states, two source sags, fixed Pstar.
% All prior baseline and certified metric files are read-only inputs.
out=fileparts(mfilename('fullpath'));root=fileparts(out);addpath(out);
assert(~isfile(fullfile(out,'results.mat')),'Use an isolated fresh directory for a new run.');
diary(fullfile(out,'console.txt'));
a=load(fullfile(root,'radial_active_transient_20261001_152918','results.mat'),'p','n');p=a.p;n=a.n;
assert(abs(p.mass-6.8)<1e-12&&p.D==25&&p.Imax==1.2&&p.Vmax==1.3&&abs(imag(n.Zeq)-.35)<1e-12);
p.Pstar=1;p.omega=100*pi;
mi=jsondecode(fileread(fullfile(root,'radial_inactive_metric_20261001_170815','metric.json')));
ma=jsondecode(fileread(fullfile(root,'radial_quadratic_robust_20261001_161903','metric.json')));
% jsondecode yields a 6-by-2-by-2 array: retain exact coefficient orientation.
mi.C=permute(mi.coefficients,[2 3 1]);ma.C=permute(ma.coefficients,[2 3 1]);
assert(max(abs(mi.C(:,:,1)-[.0019393908480797441 .00481967414426573;.00481967414426573 .9980606091519203]),[],'all')<1e-15);
terminal=jsondecode(fileread(fullfile(root,'radial_inactive_metric_20261001_170815','invariant_subset.json')));
cfg=struct('tOn',0,'tClear',.1,'tEnd',8,'RelTol',1e-9,'AbsTol',[1e-11 1e-12], ...
 'faultMaxStep',.002,'postMaxStep',.005,'switchTol',1e-9,'gamma0',[14 20.633904875622 24.373430419712 27 36], ...
 'retained',[.7 .5],'caseNames',["A","B","C","D","E"], ...
 'gammaEq',fzero(@(g)electrical_power(g,p,n)-1,[.3 .4]),'gammaSwitch',fzero(@(g)switchvalue(g,p,n),[.42 .43]));
cfg.R1=[16.426345885756*pi/180,cfg.gammaSwitch;.99652713910517932,1.0025];
cfg.R2=[cfg.gammaSwitch,33.858410380369*pi/180;.996174514620,1.005178107392];
cfg.mi=mi;cfg.ma=ma;cfg.terminal=terminal;
cfg.energyGamma=unique(sort([linspace(-pi,pi,144001) cfg.gammaEq]));
ed=coupled(cfg.energyGamma,1,p,n,0);cfg.energyW=cumtrapz(cfg.energyGamma,ed.P-1);cfg.energyW=cfg.energyW-interp1(cfg.energyGamma,cfg.energyW,cfg.gammaEq);
fprintf('Baseline Figure 3/4/5A/5B source and exports inspected before figure creation.\n');
fprintf('Old Usag=.3 meant 70%% reduction. New tests use retained .7 and .5 (30%% and 50%% reductions).\n');
fprintf('Old normalized duration=1 has no consistent direct seconds conversion. Fault=[0,.1] s, end=8 s.\n');
fprintf('Initial displacement imposed at 0-; source sag at 0+; no state reset at clearing. Pstar=1 throughout.\n');
initial=table();runs=cell(5,2);summary=table();verification=table();
for k=1:5
 x0=[cfg.gamma0(k)*pi/180;1];d=coupled(x0(1),1,p,n,0);
 state="inactive";if d.rawI>1.2+cfg.switchTol,state="active";elseif abs(d.rawI-1.2)<=cfg.switchTol,state="boundary";end
 initial=[initial;table(cfg.caseNames(k),cfg.gamma0(k),d.rawI,abs(d.ic),d.P,abs(d.u),state,'VariableNames',{'Case','gamma0_deg','rawCurrent','limitedCurrent','Pe','requiredVoltage','limiterState'})]; %#ok<AGROW>
 for j=1:2
  pf=p;pf.U=p.U*cfg.retained(j);nf=rebuild(0,pf,n); % n already contains the physical LV reactor; do not add it twice.
  assert(abs(nf.Zeq-n.Zeq)<1e-12&&abs(nf.Vth/n.Vth-cfg.retained(j))<1e-12);
  [tr,row]=experiment(x0,k,j,p,n,pf,nf,cfg,false);runs{k,j}=tr;summary=[summary;row]; %#ok<AGROW>
  % Independent tighter integration, same ten experiments and event surfaces.
  [tight,~]=experiment(x0,k,j,p,n,pf,nf,cfg,true);
  common=min(tr.t(end),tight.t(end));mainState=sample_state(tr,common);tightState=sample_state(tight,common);
  err=max(abs(mainState-tightState));terr=abs(tr.t(end)-tight.t(end));
  assert(err<2e-6&&terr<2e-5,'Tighter integration disagrees.');
  assert(tr.valid==tight.valid&&height(tr.currentEvents)==height(tight.currentEvents));
  eventError=0;if ~isempty(tr.currentEvents),eventError=max(abs(tr.currentEvents.time_s-tight.currentEvents.time_s));end
  assert(eventError<2e-5);
  verification=[verification;table(cfg.caseNames(k),100*(1-cfg.retained(j)),err,terr,eventError,'VariableNames',{'Case','sagDepth_percent','stateDifference','stopTimeDifference_s','currentEventTimeDifference_s'})]; %#ok<AGROW>
  tag=sprintf('case_%s_sag%02d',cfg.caseNames(k),round(100*(1-cfg.retained(j))));
  writetable(tr.data,fullfile(out,[tag '_trajectory.csv']));writetable(tr.currentEvents,fullfile(out,[tag '_current_events.csv']));writetable(tr.domainEvents,fullfile(out,[tag '_domain_events_sampled.csv']));writetable(tr.metricTransitions,fullfile(out,[tag '_metric_transitions.csv']));
  fprintf('%s sag=%d%% clear=(%.6f,%.9f), active %.6f s, valid=%d recovery=%s\n',cfg.caseNames(k),round(100*(1-cfg.retained(j))),row.clearGamma_deg,row.clearNu,row.activeDuration_s,tr.valid,row.recovery);
 end
end
di=coupled(cfg.gammaSwitch,1,p,n,1);da=coupled(cfg.gammaSwitch,1,p,n,2);
assert(abs(di.P-da.P)<1e-10&&abs(di.E-da.E)<1e-10&&abs(di.rawI-1.2)<1e-10);
writetable(initial,fullfile(out,'initial_conditions.csv'));writetable(summary,fullfile(out,'comparison_summary.csv'));writetable(verification,fullfile(out,'integration_verification.csv'));
save(fullfile(out,'results.mat'),'p','n','cfg','runs','summary','initial','verification','-v7.3');
disp(initial);disp(summary);diary off;
plot_baseline_figures(out);
run(fullfile(out,'audit_study.m'));

function [tr,row]=experiment(x0,k,j,p,n,pf,nf,cfg,tighter)
 rel=cfg.RelTol;ab=cfg.AbsTol;mf=cfg.faultMaxStep;mp=cfg.postMaxStep;
 if tighter,rel=1e-11;ab=[1e-13 1e-14];mf=mf/2;mp=mp/2;end
 s0=struct('time_s',{},'type',{},'stage',{},'gamma_deg',{},'nu',{});ev=s0;
 pre=coupled(x0(1),1,p,n,0);df0=coupled(x0(1),1,pf,nf,0);
 sf=segment([0 cfg.tClear],x0,pf,nf,rel,ab,mf);sp=[];clearState=[NaN NaN];reason="";
 if abs(sf.t(end)-cfg.tClear)<1e-9
  clearState=sf.x(end,:);sp=segment([cfg.tClear cfg.tEnd],clearState',p,n,rel,ab,mp);
 end
 valid=~isempty(sp)&&abs(sp.t(end)-cfg.tEnd)<1e-9;
 if ~valid,reason="voltage-command/model validity stop";end
 t=[0;sf.t];x=[x0';sf.x];stage=["prefault";repmat("fault",numel(sf.t),1)];
 if ~isempty(sp),t=[t;sp.t];x=[x;sp.x];stage=[stage;repmat("postfault",numel(sp.t),1)];end
 m=numel(t);P=zeros(m,1);rawI=P;limI=P;u=P;E=P;hE=P;Pg=P;for ss=["prefault","fault","postfault"]
  ix=stage==ss;pp=p;nn=n;if ss=="fault",pp=pf;nn=nf;end
  d=coupled(x(ix,1),1,pp,nn,0);P(ix)=d.P;rawI(ix)=d.rawI;limI(ix)=abs(d.ic);u(ix)=abs(d.u);E(ix)=d.E;hE(ix)=d.hE;Pg(ix)=d.totalPg;
 end
 assert(max(limI)<=1.2+1e-10&&min(hE)>0&&all(diff(t)>=0));
 active=rawI>1.2;boundary=abs(rawI-1.2)<=cfg.switchTol;
 if abs(pre.u)>1.3,error('Prefault initial point infeasible; do not run a fault from it.');end
 if pre.rawI<=1.2+cfg.switchTol&&df0.rawI>1.2+cfg.switchTol,ev=add_event(ev,0,"activation","source jump",x0');end
 if pre.rawI>1.2+cfg.switchTol&&df0.rawI<1.2-cfg.switchTol,ev=add_event(ev,0,"release","source jump",x0');end
 ev=segment_events(ev,sf,pf,nf,"fault");duration=sf.activeDuration;
 if ~isempty(sp)
  before=coupled(clearState(1),1,pf,nf,0);after=coupled(clearState(1),1,p,n,0);
  if before.rawI>1.2+cfg.switchTol&&after.rawI<1.2-cfg.switchTol,ev=add_event(ev,cfg.tClear,"release","clearing jump",clearState);end
  if before.rawI<1.2-cfg.switchTol&&after.rawI>1.2+cfg.switchTol,ev=add_event(ev,cfg.tClear,"activation","clearing jump",clearState);end
  ev=segment_events(ev,sp,p,n,"postfault");duration=duration+sp.activeDuration;
 end
 currentEvents=event_table(ev);
 domain=zeros(m,1);minM=nan(m,1);maxQ=minM;L1=minM;L2=minM;symMin=minM;symMax=minM;
 post=stage=="postfault";
 for ii=find(post)'
  domain(ii)=region(x(ii,:),rawI(ii),cfg);
  if domain(ii)>0
   pp=P(ii);pg=Pg(ii);if boundary(ii),side=coupled(x(ii,1),1,p,n,domain(ii));pp=side.P;pg=side.totalPg;end
   [M,MD,rate]=metric_at(x(ii,:),domain(ii),cfg,p,pp);J=jacobian(x(ii,:),p,pp,pg);Q0=J'*M+M*J+MD;
   minM(ii)=min(eig(M));maxQ(ii)=max(eig(Q0+2*rate*M));R=chol(M);ss=eig((R'\Q0/R)/2);symMin(ii)=min(ss);symMax(ii)=max(ss);
   e=x(ii,:)'-[cfg.gammaEq;1];L=sqrt(max(0,e'*M*e));if domain(ii)==1,L1(ii)=L;else,L2(ii)=L;end
  end
 end
 assert(all(minM(post&domain>0)>1e-5)&&all(maxQ(post&domain>0)<-1e-5));
 W=interp1(cfg.energyGamma,cfg.energyW,x(:,1),'pchip',NaN);w=x(:,2)-1;energy=W+6.8*p.omega*(w.^2/2+w.^3/3);
 omegaMember=post & energy<=.0063 & x(:,1)>=cfg.R1(1,1)&x(:,1)<=cfg.R2(1,2)&x(:,2)>=.9975&x(:,2)<=1.0025;
 e=x-[cfg.gammaEq 1];terminalV=sum((e*cfg.terminal.P).*e,2);terminalMember=post&terminalV<=str2double(cfg.terminal.rho);
 domainEvents=table();postix=find(post);if ~isempty(postix)
  transitions=postix([true;diff(domain(postix))~=0]);domainEvents=table(t(transitions),domain(transitions),rad2deg(x(transitions,1)),x(transitions,2),'VariableNames',{'sampleTime_s','region','gamma_deg','nu'});
 end
 if isempty(domainEvents),domainEvents=table(zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1),'VariableNames',{'sampleTime_s','region','gamma_deg','nu'});end
 % Both metric quadratic forms at each continuous postfault limiter switch.
 metricTransitions=table();if ~isempty(sp)
  for ii=1:numel(sp.sol.xe)
   if sp.sol.ie(ii)~=1,continue;end
   xx=sp.sol.ye(:,ii)';if xx(2)<cfg.R1(2,1)||xx(2)>cfg.R1(2,2),continue;end
   dd=coupled(xx(1),1,p,n,0);[m1,~,~]=metric_at(xx,1,cfg,p,dd.P);[m2,~,~]=metric_at(xx,2,cfg,p,dd.P);ee=xx'-[cfg.gammaEq;1];la=sqrt(ee'*m2*ee);li=sqrt(ee'*m1*ee);
   gain=256;if xx(2)>=.9992&&xx(2)<=1.0008,gain=4;end
   metricTransitions=[metricTransitions;table(sp.sol.xe(ii),la,li,li/max(la,eps),sqrt(gain),'VariableNames',{'time_s','L_active','L_inactive','chordRatio_activeToInactive','validatedLengthGainBound'})]; %#ok<AGROW>
  end
 end
 if isempty(metricTransitions),metricTransitions=table(zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1),'VariableNames',{'time_s','L_active','L_inactive','chordRatio_activeToInactive','validatedLengthGainBound'});end
 tail=post&t>=.9*cfg.tEnd;recovery="N";
 if valid&&any(tail)&&all(abs(rad2deg(x(tail,1)-cfg.gammaEq))<.1&abs(w(tail))<1e-4)&&abs(P(end)-1)<1e-4,recovery="Y";end
 if ~valid,recovery="Not assessed: validity stop";end
 firstOmega=first_time(t,omegaMember);firstTerminal=first_time(t,terminalMember);contained=~isempty(postix)&&all(domain(postix)>0);
 clearRegion="not reached";if ~isempty(sp),clearRegion=region_name(domain(postix(1)));end
 exitTime=NaN;reentry=NaN;if ~isempty(postix),inside=domain(postix)>0;tt=t(postix);ex=find(diff(inside)==-1,1);if ~isempty(ex),exitTime=tt(ex+1);re=find(diff(inside(ex+1:end))==1,1);if ~isempty(re),reentry=tt(ex+re+1);end;end;end
 data=table(t,stage,rad2deg(x(:,1)),x(:,2),50*w,P,rawI,limI,active,boundary,u,1.3-u,E,domain,minM,maxQ,symMin,symMax,L1,L2,energy,omegaMember,terminalMember, ...
  'VariableNames',{'time_s','stage','gamma_deg_unwrapped','nu','frequencyDeviation_Hz','Pe_pu','rawCurrent_pu','limitedCurrent_pu','limiterActive','atSwitchTolerance','requiredVoltage_pu','voltageMargin_pu','Eref_pu','certifiedRegion','lambdaMinM','lambdaMaxQc','lambdaMinSymF','lambdaMaxSymF','L_inactive_chord','L_active_chord','energyNumerical','inOmegaNumerical','inTerminalNumerical'});
 activated=any(active(stage~="prefault"));tr=struct('t',t,'x',x,'data',data,'fault',sf,'post',sp,'valid',valid,'reason',reason,'currentEvents',currentEvents,'domainEvents',domainEvents,'metricTransitions',metricTransitions,'caseIndex',k,'sagIndex',j,'retained',cfg.retained(j),'clearState',clearState);
 row=table(cfg.caseNames(k),cfg.gamma0(k),100*(1-cfg.retained(j)),rad2deg(clearState(1)),clearState(2),activated,duration,valid,clearRegion,contained,isfinite(firstOmega),isfinite(firstTerminal),recovery,rad2deg(x(end,1)),x(end,2), ...
  max(rad2deg(x(:,1))),min(rad2deg(x(:,1))),max(50*x(:,2)),min(50*x(:,2)),max(u),t(end),rad2deg(x(end,1)-cfg.gammaEq),x(end,2)-1,P(end)-1,firstOmega,firstTerminal,exitTime,reentry,reason, ...
  'VariableNames',{'Case','gamma0_deg','sagDepth_percent','clearGamma_deg','clearNu','limiterActivated','activeDuration_s','voltageFeasible','clearRegion','postfaultContained_sampled','entersOmega_sampled','entersTerminal_sampled','recovery','finalGamma_deg','finalNu','maxGamma_deg','minGamma_deg','maxFrequency_Hz','minFrequency_Hz','maxRequiredVoltage','lastValidTime_s','finalGammaError_deg','finalNuError','finalPowerResidual','firstOmegaEntry_sampled_s','firstTerminalEntry_sampled_s','firstDomainExit_sampled_s','firstDomainReentry_sampled_s','modelValidityNote'});
end

function seg=segment(span,x0,p,n,rel,ab,maxstep)
 d=coupled(x0(1),1,p,n,0);
 if abs(d.u)>p.Vmax||x0(2)<=.5
  seg=struct('t',span(1),'x',x0','sol',struct('xe',[],'ye',zeros(2,0),'ie',[]),'activeDuration',0);return;
 end
 opt=odeset('RelTol',rel,'AbsTol',ab,'MaxStep',maxstep,'Events',@(t,x)events(t,x,p,n));
 sol=ode45(@(t,x)rhs(x,p,n),span,x0,opt);
 count=6001;if span(1)==0,count=501;end
 t=unique(sort([linspace(span(1),sol.x(end),count) sol.xe]));x=deval(sol,t)';
 cut=unique(sort([span(1),sol.xe(sol.ie==1),sol.x(end)]));duration=0;
 for j=1:numel(cut)-1,mid=mean(cut(j:j+1));xx=deval(sol,mid);dd=coupled(xx(1),1,p,n,0);if dd.rawI>1.2,duration=duration+diff(cut(j:j+1));end;end
 seg=struct('t',t','x',x,'sol',sol,'activeDuration',duration);
end
function e=segment_events(e,s,p,n,stage)
 for j=1:numel(s.sol.xe)
  if s.sol.ie(j)~=1,continue;end
  xx=s.sol.ye(:,j)';h=1e-6;gp=coupled(xx(1)+h,1,p,n,0);gm=coupled(xx(1)-h,1,p,n,0);sd=(gp.rawI-gm.rawI)/(2*h)*p.omega*(xx(2)-1);type="activation";if sd<0,type="release";end
  if abs(sd)<1e-8,type="tangent";end;e=add_event(e,s.sol.xe(j),type,stage,xx);
 end
end
function e=add_event(e,t,type,stage,x)
 e(end+1)=struct('time_s',t,'type',type,'stage',stage,'gamma_deg',rad2deg(x(1)),'nu',x(2));
end
function T=event_table(e)
 if isempty(e),T=table(zeros(0,1),strings(0,1),strings(0,1),zeros(0,1),zeros(0,1),'VariableNames',{'time_s','type','stage','gamma_deg','nu'});else,T=sortrows(struct2table(e),'time_s');end
end
function x=sample_state(tr,t)
 if ~isempty(tr.post)&&t>=tr.post.t(1)&&isfield(tr.post.sol,'x'),x=deval(tr.post.sol,t)';elseif isfield(tr.fault.sol,'x'),x=deval(tr.fault.sol,t)';else,x=tr.x(end,:);end
end
function rr=region(x,raw,cfg)
 rr=0;tol=cfg.switchTol;
 b=cfg.R1;if raw<=1.2+tol&&x(1)>=b(1,1)-1e-12&&x(1)<=b(1,2)+1e-12&&x(2)>=b(2,1)&&x(2)<=b(2,2),rr=1;end
 b=cfg.R2;if raw>1.2+tol&&x(1)>=b(1,1)-1e-12&&x(1)<=b(1,2)+1e-12&&x(2)>=b(2,1)&&x(2)<=b(2,2),rr=2;end
end
function s=region_name(r)
 s="outside";if r==1,s="R1";elseif r==2,s="R2";end
end
function [M,MD,rate]=metric_at(x,branch,cfg,p,Pe)
 m=cfg.mi;if branch==2,m=cfg.ma;end
 z=(x-m.center')./m.scale';f=rhs_power(x,p,Pe);dz=f'./m.scale';b=[1 z z(1)^2 z(1)*z(2) z(2)^2];bd=[0 dz 2*z(1)*dz(1) z(2)*dz(1)+z(1)*dz(2) 2*z(2)*dz(2)];M=zeros(2);MD=M;
 for j=1:6,M=M+b(j)*m.C(:,:,j);MD=MD+bd(j)*m.C(:,:,j);end;rate=m.lambda;
end
function J=jacobian(x,p,P,Pg)
 J=[0 p.omega;-Pg/(p.mass*x(2)),(-(1-P)/x(2)^2-p.D)/p.mass];
end
function t=first_time(tt,mask)
 t=NaN;j=find(mask,1);if ~isempty(j),t=tt(j);end
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






