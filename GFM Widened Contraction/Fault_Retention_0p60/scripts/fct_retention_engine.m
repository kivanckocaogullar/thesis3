function varargout=fct_retention_engine(action,varargin)
% Shared extension of the original five-FCT script. Physics/recovery functions
% below are copied verbatim; no tolerance reruns, searches or optimization.
switch action
 case 'run',run_cases(varargin{:});
 case 'postprocess',postprocess(varargin{:});
 case 'electrical',varargout{1}=coupled(varargin{1},1,varargin{2},varargin{3},0);
 case 'combined',combined(varargin{:});
 case 'currentplots',currentplots(varargin{:});
 otherwise,error('Unknown action');
end
end
function run_cases(root,Uret,FCT)
input=fullfile(fileparts(mfilename('fullpath')),'inputs');out=fullfile(root,'results');
assert(~isfile(fullfile(out,'fct_analysis.mat')),'Saved results exist: use plots/postprocess; no automatic rerun.');
base=load(fullfile(input,'baseline_results.mat'));a=load(fullfile(input,'model_inputs.mat'));w=load(fullfile(input,'region_model_inputs.mat'));ref=load(fullfile(input,'reference_040.mat'),'p','n','cfg','x0');
p=base.p;n=base.n;cfg=base.cfg;a.p.Pstar=1;a.p.omega=100*pi;w.p.Pstar=1;w.p.omega=100*pi;
assert(isequaln(p,a.p)&&isequaln(n,a.n)&&isequaln(p,w.p)&&isequaln(n,w.n)&&isequaln(p,ref.p)&&isequaln(n,ref.n));
for f={'tEnd','RelTol','AbsTol','faultMaxStep','postMaxStep'},assert(isequaln(cfg.(f{1}),ref.cfg.(f{1})));end
cfg.gammaEq=fzero(@(g)electrical_power(g,p,n)-p.Pstar,[.3 .4],optimset('TolX',1e-12));x0=[cfg.gammaEq;1];
assert(norm(x0-ref.x0,inf)<1e-12&&norm(rhs(x0,p,n),inf)<1e-10);
assert(cfg.tEnd==8&&p.mass==6.8&&p.D==25&&p.Imax==1.2&&p.Vmax==1.3);
if isfield(cfg,'tClear'),cfg=rmfield(cfg,'tClear');end % assess receives each FCT explicitly.
cfg.selectedFCT=FCT;cfg.retention=Uret;cfg.faultApplicationTime_s=0;cfg.observationConvention='8 seconds from fault inception (reference convention)';
pf=p;pf.U=p.U*Uret;nf=rebuild(0,pf,n);
metrics=cell(2,1);certificates=metrics;names={'inactive','active'};
for j=1:2
 m=load(fullfile(input,names{j},'metric.mat'));orig=load(fullfile(input,['original_' names{j} '_metric.mat']));
 assert(isequaln(m.C,orig.C)&&isequaln(m.bounds,orig.bounds)&&m.lambda==orig.lambda);
 cert=jsondecode(fileread(fullfile(input,names{j},'interval_verification.json')));
 jm=jsondecode(fileread(fullfile(input,names{j},'metric.json')));
 assert(strcmp(cert.status,'VERIFIED')&&cert.unresolved_boxes==0&&cert.counterexample_boxes==0);
 assert(max(abs(m.C-permute(jm.coefficients,[2 3 1])),[],'all')<1e-12);
 assert(jm.U==p.U&&jm.mass==p.mass&&jm.D==p.D&&jm.Imax==p.Imax&&jm.umax==p.Vmax&&jm.AD==p.AD&&abs(jm.omega-p.omega)<1e-12);
 f=fieldnames(jm.network);for q=1:numel(f),z=jm.network.(f{q});assert(abs(n.(f{q})-complex(z(1),z(2)))<1e-12);end
 expected=[10 24.3234;.98 1.02]';if j==2,expected=[24.4234 43;.98 1.02]';end
 assert(max(abs(rad2deg(m.bounds(:,1))-expected(:,1)))<1e-8&&max(abs(m.bounds(:,2)-expected(:,2)))<1e-12);
 rates=[1 .1];assert(m.lambda==rates(j)&&m.upper==10);metrics{j}=m;certificates{j}=cert;
end
cct=readtable(fullfile(input,'cct_results.csv'));cct=cct(abs(cct.Uret-Uret)<1e-12,:);assert(height(cct)==1);
writetable(cct,fullfile(out,'cct_reference.csv'));
switchGamma=fzero(@(g)switchvalue(g,p,n),deg2rad([24 25]),optimset('TolX',1e-12));
results=table();runs=cell(numel(FCT),1);
for k=1:numel(FCT)
 fprintf('RUN Uret=%.2f FCT=%.2f (%d/5)\n',Uret,FCT(k),k);
 [q,runs{k}]=assess(FCT(k),cfg,false,x0,p,n,pf,nf);results=[results;struct2table(q)];
 fprintf(' -> %s gamma=%.12f nu=%.12f\n',q.classification,q.gamma_clear_deg,q.nu_clear);
end
save(fullfile(out,'fct_analysis.mat'),'results','runs','p','n','pf','nf','cfg','x0','Uret','FCT','metrics','certificates','cct','switchGamma','-v7.3');
writetable(results,fullfile(out,'assessment_details.csv'));
postprocess(root);
end
function postprocess(root)
out=fullfile(root,'results');a=load(fullfile(out,'fct_analysis.mat'));details=table();limiterIntervals=table();
for k=1:numel(a.FCT)
 r=a.runs{k};q=a.results(k,:);tc=a.FCT(k);
 row=struct('Uret',a.Uret,'VoltageDrop_percent',100*(1-a.Uret),'FCT_s',tc,'GammaClear_deg',q.gamma_clear_deg,'NuClear_pu',q.nu_clear, ...
 'MaxPostfaultGammaUnwrapped_deg',NaN,'MinNu_pu',NaN,'MaxNu_pu',NaN,'MaxRequestedCurrent_pu',NaN,'MaxActualCurrent_pu',NaN,'MaxVoltageCommand_pu',NaN, ...
 'LimiterEverActive',false,'LimiterActiveAtClearing',false,'FaultLimiterDuration_s',NaN,'PostfaultLimiterDuration_s',NaN, ...
 'Recovery',q.classification,'FinalGamma_deg',q.final_gamma_deg,'FinalNu_pu',q.final_nu,'InsideInactive',false,'InsideActive',false,'ContractionRegion',"OUTSIDE", ...
 'AllPostfaultSamplesInside',false,'LeavesClearingModeRectangle',false,'FirstClearingModeExit_s_after_clear',NaN,'LeavesCertifiedUnion',false,'FirstUnionExit_s_after_clear',NaN,'FirstOutsideUnion_s_after_clear',NaN,'TimeInside_percent',NaN,'MinMetricEigAtClearing',NaN,'MaxMetricEigAtClearing',NaN,'MetricConditionAtClearing',NaN);
 if isfield(r,'postfault')&&any(q.classification==["RECOVERED","NOT_RECOVERED"])
 sf=r.fault;sp=r.postfault;df=coupled(sf.x(:,1),1,a.pf,a.nf,0);dp=coupled(sp.x(:,1),1,a.p,a.n,0);
 [ii,ia]=membership(sp.x,dp,a.metrics);inside=ii|ia;
 row.MaxPostfaultGammaUnwrapped_deg=max(rad2deg(sp.x(:,1)));row.MinNu_pu=min([sf.x(:,2);sp.x(:,2)]);row.MaxNu_pu=max([sf.x(:,2);sp.x(:,2)]);
 row.MaxRequestedCurrent_pu=max([df.rawI;dp.rawI]);row.MaxActualCurrent_pu=max(abs([df.ic;dp.ic]));row.MaxVoltageCommand_pu=max(abs([df.u;dp.u]));
 row.LimiterEverActive=any([df.active;dp.active]);row.LimiterActiveAtClearing=dp.active(1);row.FaultLimiterDuration_s=sf.activeDuration;row.PostfaultLimiterDuration_s=sp.activeDuration;
 row.InsideInactive=ii(1);row.InsideActive=ia(1);row.AllPostfaultSamplesInside=all(inside);
 if ii(1),row.ContractionRegion="INACTIVE";elseif ia(1),row.ContractionRegion="ACTIVE";end
 ix=find(~inside,1);if ~isempty(ix),row.FirstOutsideUnion_s_after_clear=sp.t(ix)-tc;end
 if inside(1)
  row.LeavesCertifiedUnion=any(~inside);row.FirstUnionExit_s_after_clear=row.FirstOutsideUnion_s_after_clear;
  own=ii;if ia(1),own=ia;end;ix=find(~own,1);row.LeavesClearingModeRectangle=~isempty(ix);
  if ~isempty(ix),row.FirstClearingModeExit_s_after_clear=sp.t(ix)-tc;end
 end
 row.TimeInside_percent=100*trapz(sp.t,double(inside))/(sp.t(end)-tc);
 m=a.metrics{1+double(dp.active(1))};ev=metric_eigs(sp.x(1,:),m);row.MinMetricEigAtClearing=min(ev);row.MaxMetricEigAtClearing=max(ev);row.MetricConditionAtClearing=max(ev)/min(ev);
 for stage=1:2
  s=sf;d=df;p=a.pf;n=a.nf;label="FAULT";if stage==2,s=sp;d=dp;p=a.p;n=a.n;label="POSTFAULT";end
  tt=table(s.t,rad2deg(s.x(:,1)),s.x(:,2),d.rawI,abs(d.ic),abs(d.u),d.active,'VariableNames',{'Time_s','GammaUnwrapped_deg','Nu_pu','RequestedCurrent_pu','ActualCurrent_pu','VoltageCommand_pu','LimiterActive'});
  if stage==2,tt.InsideInactive=ii;tt.InsideActive=ia;else,tt.InsideInactive=false(height(tt),1);tt.InsideActive=false(height(tt),1);end
  writetable(tt,fullfile(out,sprintf('trajectory_FCT_%.2f_%s.csv',tc,label)));
  cuts=unique(sort([s.t(1),s.sol.xe(s.sol.ie==1),s.t(end)]));
  for z=1:numel(cuts)-1
   xx=deval(s.sol,mean(cuts(z:z+1)));dd=coupled(xx(1),1,p,n,0);
   if dd.active
    tr=table(a.Uret,tc,label,cuts(z),cuts(z+1),cuts(z+1)-cuts(z),'VariableNames',{'Uret','FCT_s','Stage','Start_s','End_s','Duration_s'});limiterIntervals=[limiterIntervals;tr];
   end
  end
 end
 end
 details=[details;struct2table(row)];
end
summary=details(:,{'FCT_s','GammaClear_deg','NuClear_pu','Recovery','LimiterActiveAtClearing','ContractionRegion'});
writetable(details,fullfile(out,'case_details.csv'));writetable(summary,fullfile(out,'five_FCT_summary.csv'));
if isempty(limiterIntervals),limiterIntervals=table('Size',[0 6],'VariableTypes',{'double','double','string','double','double','double'},'VariableNames',{'Uret','FCT_s','Stage','Start_s','End_s','Duration_s'});end
writetable(limiterIntervals,fullfile(out,'limiter_active_intervals.csv'));save(fullfile(out,'fct_analysis.mat'),'details','summary','limiterIntervals','-append');
export_configuration(a,out);
disp(summary);
end
function export_configuration(a,out)
selected=struct('mass',a.p.mass,'D',a.p.D,'Imax',a.p.Imax,'Vmax',a.p.Vmax,'Pstar',a.p.Pstar,'omega',a.p.omega);
config=struct('cfg',a.cfg,'real_model_parameters',selected,'initial_state',a.x0,'nominal_restored_voltage',a.p.U,'fault_voltage',a.pf.U,'switch_angle_deg',rad2deg(a.switchGamma));
fid=fopen(fullfile(out,'numerical_configuration.json'),'w','n','UTF-8');cleanup=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(config,'PrettyPrint',true));
end
function [ii,ia]=membership(x,d,metrics)
ii=~d.active&all(x>=metrics{1}.bounds(1,:)&x<=metrics{1}.bounds(2,:),2);
ia=d.active&all(x>=metrics{2}.bounds(1,:)&x<=metrics{2}.bounds(2,:),2);
end
function ev=metric_eigs(x,m)
z=(x-m.center)./m.scale;b=[1 z z(1)^2 z(1)*z(2) z(2)^2];M=zeros(2);for j=1:6,M=M+b(j)*m.C(:,:,j);end;ev=eig((M+M')/2);
end
function combined(base)
shared=fileparts(mfilename('fullpath'));ref=load(fullfile(shared,'inputs','reference_040.mat'));new=load(fullfile(base,'Fault_Retention_0p60','results','fct_analysis.mat'));
combinedTable=table();
for k=1:numel(ref.FCT)
 q=ref.results(k,:);x=ref.runs{k}.postfault.x(1,:);d=coupled(x(1),1,ref.p,ref.n,0);[ii,ia]=membership(x,d,new.metrics);reg="OUTSIDE";if ii,reg="INACTIVE";elseif ia,reg="ACTIVE";end
 combinedTable=[combinedTable;table(ref.Uret,100*(1-ref.Uret),ref.FCT(k),q.gamma_clear_deg,q.nu_clear,q.classification,d.active,reg,'VariableNames',{'Uret','VoltageDrop_percent','FCT_s','GammaClear_deg','NuClear_pu','Recovery','LimiterActiveAtClearing','ContractionRegion'})];
end
for code={'60','80'}
 t=readtable(fullfile(base,['Fault_Retention_0p' code{1}],'results','case_details.csv'),'TextType','string');combinedTable=[combinedTable;t(:,combinedTable.Properties.VariableNames)];
end
for code={'60','80'}
 out=fullfile(base,['Fault_Retention_0p' code{1}],'results');writetable(combinedTable,fullfile(out,'combined_040_060_080.csv'));save(fullfile(out,'combined_040_060_080.mat'),'combinedTable');
end
disp(combinedTable);
end
function currentplots(root)
out=fullfile(root,'results');a=load(fullfile(out,'fct_analysis.mat'));colors=[0 .45 .70;0 .62 .45;.80 .47 .05;.65 .25 .65;.85 .25 .10];
f=figure('Visible','off','Color','w','Position',[60 60 1300 1050]);tl=tiledlayout(f,4,1,'TileSpacing','compact','Padding','compact');
for kind=1:4
 ax=nexttile(tl);hold(ax,'on');hs=gobjects(5,1);labs=strings(5,1);
 for k=1:5
  r=a.runs{k};df=coupled(r.fault.x(:,1),1,a.pf,a.nf,0);dp=coupled(r.postfault.x(:,1),1,a.p,a.n,0);
  switch kind,case 1,yf=df.rawI;yp=dp.rawI;case 2,yf=abs(df.ic);yp=abs(dp.ic);case 3,yf=abs(df.u);yp=abs(dp.u);case 4,yf=double(df.active);yp=double(dp.active);end
  % NaN break at voltage restoration: no artificial diagonal through a jump.
  t=[r.fault.t;NaN;r.postfault.t];y=[yf;NaN;yp];ss='-';if a.results.classification(k)~="RECOVERED",ss='--';end
  hs(k)=plot(ax,t,y,ss,'Color',colors(k,:),'LineWidth',1.2);labs(k)=sprintf('FCT %.2f s',a.FCT(k));
 end
 if kind<=2,yline(ax,a.p.Imax,'k:','I_{max}','HandleVisibility','off');elseif kind==3,yline(ax,a.p.Vmax,'k:','V_{max}','HandleVisibility','off');else,ylim(ax,[-.1 1.1]);yticks(ax,[0 1]);end
 yl={'Requested current (pu)','Actual current (pu)','Voltage command (pu)','Limiter active'};ylabel(ax,yl{kind});grid(ax,'on');xlim(ax,[0 a.cfg.tEnd]);set(ax,'FontName','Arial','FontSize',11);
 if kind==1,legend(ax,hs,labs,'Location','northoutside','NumColumns',5);end
end
xlabel(ax,'Time from fault inception (s)');title(tl,sprintf('Current-limiter behavior | U_{ret}=%.2f pu (%.0f%% voltage drop)',a.Uret,100*(1-a.Uret)));
export_three(f,fullfile(root,'figures'),'05_current_limiter');close(f);
end
function export_three(fig,out,name)
tmp=tempname(fullfile(getenv('SystemRoot'),'Temp'));mkdir(tmp);cleanup=onCleanup(@()rmdir(tmp,'s'));
drawnow;exportgraphics(fig,fullfile(tmp,'plot.png'),'Resolution',300);exportgraphics(fig,fullfile(tmp,'plot.pdf'),'ContentType','vector');savefig(fig,fullfile(tmp,'plot.fig'));
for ext={'.png','.pdf','.fig'},copyfile(fullfile(tmp,['plot' ext{1}]),fullfile(out,[name ext{1}]));end
end

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



















