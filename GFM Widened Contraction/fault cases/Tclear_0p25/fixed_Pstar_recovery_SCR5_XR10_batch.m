%% Fixed-Pstar nominal-equilibrium retained-voltage sweep (physical seconds).
clearvars;clc;
out=fileparts(mfilename('fullpath'));
a=load(fullfile(out,'model_inputs.mat'),'p','n');p=a.p;n=a.n;
save(fullfile(out,'model_inputs.mat'),'p','n');
p.Pstar=1;p.omega=100*pi;
cfg=struct('tClear',.25,'tEnd',8,'RelTol',1e-9,'AbsTol',[1e-11 1e-12],'faultMaxStep',.002,'postMaxStep',.005);
cfg.gammaEq=fzero(@(g)electrical_power(g,p,n)-p.Pstar,[.3 .4]);
assert(abs(rad2deg(cfg.gammaEq)-20.633904875622)<1e-10);
cfg.gammaEq=deg2rad(20.633904875622);x0=[cfg.gammaEq;1];assert(norm(rhs(x0,p,n),inf)<1e-10);
Uret=[.9 .8 .7 .65 .6 .55 .5 .45 .4 .35 .3 .25];
[Cpost,Cinactive]=load_widened_regions();
summary=table();verification=table();runs=cell(numel(Uret),1);
for j=1:numel(Uret)
 pf=p;pf.U=p.U*Uret(j);nf=rebuild(0,pf,n);
 [r,row]=one_case(x0,p,n,pf,nf,cfg,Cpost,false);runs{j}=r;
 [rr,rt]=one_case(x0,p,n,pf,nf,cfg,Cpost,true);
 assert(row.classification==rt.classification&&row.inside_Cpost==rt.inside_Cpost&&row.limiter_encountered==rt.limiter_encountered);
 err=max(abs(r.clearState-rr.clearState));fprintf('Clearing verification error %.12g\n',err);
 row.Uret=Uret(j);row.sag_depth_percent=100*(1-Uret(j));row.Pstar=p.Pstar;row.gamma0_deg=rad2deg(x0(1));row.nu0=x0(2);row.Tclear_s=cfg.tClear;
 row.inside_inactive_Cpost=~row.limiter_active_postclear&&row.gamma_clear_deg>=Cinactive(1,1)&&row.gamma_clear_deg<=Cinactive(1,2)&&row.nu_clear>=Cinactive(2,1)&&row.nu_clear<=Cinactive(2,2);row.inside_any_Cpost=row.inside_Cpost||row.inside_inactive_Cpost;
 summary=[summary;struct2table(row)];
 verification=[verification;table(Uret(j),err,abs(row.max_gamma_deg-rt.max_gamma_deg),abs(row.max_raw_current_pu-rt.max_raw_current_pu),'VariableNames',{'Uret','clearing_state_max_error','max_gamma_difference_deg','max_raw_current_difference_pu'})];
 writetable(r.data,fullfile(out,sprintf('trajectory_Uret_%.2f.csv',Uret(j))));
 writetable(r.currentEvents,fullfile(out,sprintf('limiter_events_Uret_%.2f.csv',Uret(j))));
 fprintf('Uret %.2f: clear %.9f deg / %.9f; %s; Cpost %d\n',Uret(j),row.gamma_clear_deg,row.nu_clear,row.classification,row.inside_Cpost);
end
assert(all(summary.gamma0_deg==rad2deg(x0(1)))&&all(summary.Pstar==1)&&all(summary.nu0==1));
writetable(summary,fullfile(out,'fixed_Pstar_sag_sweep.csv'));writetable(verification,fullfile(out,'integration_verification.csv'));
save(fullfile(out,'fixed_Pstar_sag_sweep_results.mat'),'p','n','cfg','x0','Uret','Cpost','Cinactive','summary','verification','runs','-v7.3');
fid=fopen(fullfile(out,'fixed_Pstar_sag_sweep_summary.txt'),'w');guard=onCleanup(@()fclose(fid));
fprintf(fid,'Fixed-Pstar nominal-equilibrium sag-depth assessment\nPstar = 1 throughout; gamma0 = %.12f deg; nu0 = 1; Tclear = %.6g s.\n',rad2deg(x0(1)),cfg.tClear);
fprintf(fid,'Model: saved radial Xeq=0.35, SCR5/XR10, algebraic E+0.1Q=1 coupling, radial Imax=1.2. Source U is reduced during fault; Q-V reference stays 1. omega=100*pi rad/s, mass=6.8, D=25.\n');
fprintf(fid,'Widened active Cpost: gamma=[%.8g,%.8g] deg, nu=[%.8g,%.8g]. Correct-mode membership is not a recovery certificate.\n',Cpost(1,:),Cpost(2,:));
fprintf(fid,'Postfault observation ends at 8 s. Recovered: last 10%% of postfault samples within 0.1 deg and 1e-4 nu of the SAME nominal equilibrium, final |Pe-1|<1e-4, and no validity stop. Other complete runs are non-recovered within horizon. Validity stops are not assessed.\n');
fprintf(fid,'Voltage-command >1.3 or nu<=0.5 ends model-valid observation. Maxima use adaptive nodes plus 501 fault / 6001 postfault dense samples and events; numerical extrema only.\n');
fprintf(fid,'Clearing activity is recorded immediately before and after cfg.tClear. Release time is first observed release; last release and all event files also exported. NaN means unavailable/not encountered. Times are absolute seconds from fault onset.\n');
fprintf(fid,'Every case repeated at RelTol=1e-11, AbsTol=[1e-13,1e-14], half maximum steps; classifications, membership and limiter encounter agree. Maximum clearing-state difference %.3g.\n\n',max(verification.clearing_state_max_error));
fprintf(fid,'Uret sag%% gamma_clear_deg nu_clear classification limiter clear- clear+ first_release_s inside_Cpost max_gamma_deg min_nu max_nu max_raw_I max_limited_I\n');
for j=1:height(summary),s=summary(j,:);fprintf(fid,'%.2f %.0f %.9f %.9f %s %d %d %d %.9g %d %.9f %.9f %.9f %.9f %.9f\n',s.Uret,s.sag_depth_percent,s.gamma_clear_deg,s.nu_clear,s.classification,s.limiter_encountered,s.limiter_active_at_clearing,s.limiter_active_postclear,s.limiter_release_time_s,s.inside_Cpost,s.max_gamma_deg,s.min_nu,s.max_nu,s.max_raw_current_pu,s.max_limited_current_pu);end
rec=summary.recovered;assessed=summary.classification~="NOT_ASSESSED_VALIDITY_STOP";
if all(rec),fprintf(fid,'\nAll requested cases recovered. Lowest tested retained voltage %.2f (maximum tested sag %.0f%%) recovered; no critical threshold bracket established below this value.\n',min(Uret),100*(1-min(Uret)));
elseif any(rec)&&any(assessed&~rec)
 lo=min(summary.Uret(rec));outer=max(summary.Uret(assessed&~rec&summary.Uret<lo));
 if ~isempty(outer),fprintf(fid,'\nSampled transition: retained voltage in [%.2f, %.2f], lower non-recovered and upper recovered; sag-depth in [%.0f%%, %.0f%%]. This is a grid bracket at the configured clearing time, not an exact critical voltage.\n',outer,lo,100*(1-lo),100*(1-outer));end
else,fprintf(fid,'\nNo assessed recovery/non-recovery bracket was established.\n');end
fprintf(fid,'This sweep assesses critical retained-voltage / maximum sag-depth at the configured clearing time. Initial angle and Pstar never vary.\n');
disp(summary);
comparisonScript=fullfile(fileparts(out),'compare_clearing_times.m');if isfile(comparisonScript),run(comparisonScript);end

function [r,row]=one_case(x0,p,n,pf,nf,cfg,Cpost,tighter)
rel=cfg.RelTol;ab=cfg.AbsTol;mf=cfg.faultMaxStep;mp=cfg.postMaxStep;
if tighter,rel=1e-11;ab=[1e-13 1e-14];mf=mf/2;mp=mp/2;end
sf=segment([0 cfg.tClear],x0,pf,nf,rel,ab,mf);sp=[];clear=[NaN NaN];
ev=struct('time_s',{},'type',{},'stage',{},'gamma_deg',{},'nu',{});
d0=coupled(x0(1),1,p,n,0);df=coupled(x0(1),1,pf,nf,0);
if ~d0.active&&df.active,ev=add_event(ev,0,"activation","source jump",x0');end
ev=segment_events(ev,sf,pf,nf,"fault");
if abs(sf.t(end)-cfg.tClear)<1e-9
 clear=sf.x(end,:);before=coupled(clear(1),1,pf,nf,0);after=coupled(clear(1),1,p,n,0);
 sp=segment([cfg.tClear cfg.tEnd],clear',p,n,rel,ab,mp);
 if before.active&&~after.active,ev=add_event(ev,cfg.tClear,"release","clearing jump",clear);end
 if ~before.active&&after.active,ev=add_event(ev,cfg.tClear,"activation","clearing jump",clear);end
 ev=segment_events(ev,sp,p,n,"postfault");
else,before=struct('active',false);after=before;end
t=sf.t;x=sf.x;d=coupled(x(:,1),1,pf,nf,0);rawI=d.rawI;limited=abs(d.ic);u=abs(d.u);P=d.P;stage=repmat("fault",numel(t),1);
if ~isempty(sp)
 dp=coupled(sp.x(:,1),1,p,n,0);t=[t;sp.t];x=[x;sp.x];rawI=[rawI;dp.rawI];limited=[limited;abs(dp.ic)];u=[u;abs(dp.u)];P=[P;dp.P];stage=[stage;repmat("postfault",numel(sp.t),1)];
end
valid=~isempty(sp)&&abs(sp.t(end)-cfg.tEnd)<1e-9;
tail=stage=="postfault"&t>=.9*cfg.tEnd;
rec=valid&&any(tail)&&all(abs(rad2deg(x(tail,1)-cfg.gammaEq))<.1&abs(x(tail,2)-1)<1e-4)&&abs(P(end)-1)<1e-4;
class="NON_RECOVERED_WITHIN_HORIZON";if rec,class="RECOVERED";end;if ~valid,class="NOT_ASSESSED_VALIDITY_STOP";end
events=event_table(ev);releases=events.time_s(events.type=="release");first=NaN;last=NaN;if ~isempty(releases),first=releases(1);last=releases(end);end
inside=after.active && rad2deg(clear(1))>=Cpost(1,1)&&rad2deg(clear(1))<=Cpost(1,2)&&clear(2)>=Cpost(2,1)&&clear(2)<=Cpost(2,2);
row=struct('gamma_clear_deg',rad2deg(clear(1)),'nu_clear',clear(2),'recovered',rec,'classification',class,'limiter_encountered',any(rawI>p.Imax),'limiter_active_at_clearing',before.active,'limiter_active_postclear',after.active,'limiter_release_time_s',first,'limiter_last_release_time_s',last,'inside_Cpost',inside,'max_gamma_deg',max(rad2deg(x(:,1))),'min_nu',min(x(:,2)),'max_nu',max(x(:,2)),'max_raw_current_pu',max(rawI),'max_limited_current_pu',max(limited),'max_required_voltage_pu',max(u),'model_valid',valid,'observation_end_s',t(end),'final_gamma_deg',rad2deg(x(end,1)),'final_nu',x(end,2));
data=table(t,stage,rad2deg(x(:,1)),x(:,2),P,rawI,limited,u,'VariableNames',{'time_s','stage','gamma_deg','nu','Pe','raw_current_pu','limited_current_pu','required_voltage_pu'});
r=struct('clearState',clear,'fault',sf,'post',sp,'data',data,'currentEvents',events);
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












