clearvars;clc;
out=fileparts(mfilename('fullpath'));a=load(fullfile(out,'fixed_Pstar_sag_sweep_results.mat'));cfg=a.cfg;
cfg.RelTol=1e-12;cfg.AbsTol=[1e-14 1e-15];cfg.faultMaxStep=.0005;cfg.postMaxStep=.00125;
check=table();
for j=find(a.verification.clearing_state_max_error>2e-7)'
 pf=a.p;pf.U=a.p.U*a.Uret(j);nf=rebuild(0,pf,a.n);
 [r,s]=one_case(a.x0,a.p,a.n,pf,nf,cfg,a.Cpost,false);
 err=max(abs(r.clearState-a.runs{j}.clearState));
 assert(s.classification==a.summary.classification(j)&&s.inside_Cpost==a.summary.inside_Cpost(j));
 assert(err<2e-6);
 check=[check;table(a.Uret(j),err,s.classification,'VariableNames',{'Uret','main_vs_extra_tight_clear_state_error','classification'})];
end
writetable(check,fullfile(out,'extra_tight_verification.csv'));disp(check);
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













