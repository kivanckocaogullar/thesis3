function plot_baseline_figures(out)
% Layouts replicated from baseline_gfm_contraction_SCR5_XR10.m:
% base_figures (3 stacked spectra), trajectory_study (2 x N),
% figure5_views (region bars and 3 x 2 temporal overview + four examples).
% State-dependent M requires Mdot; no equilibrium-family Pstar substitution.
load(fullfile(out,'results.mat'),'p','n','cfg','runs','summary');colors=lines(5);
set(groot,'defaultAxesFontName','Arial','defaultAxesFontSize',10,'defaultLineLineWidth',1.5);
labels=cfg.caseNames+" : "+compose('%.4f',cfg.gamma0)+" deg";
%% Figure 3: same J / M / Sym(F) three-panel matrix analysis.
gd=unique(sort([linspace(0,180,3601) cfg.gamma0 rad2deg(cfg.R1(1,:)) rad2deg(cfg.R2(1,:))]));gg=gd*pi/180;dd=coupled(gg,1,p,n,0);
ej=nan(numel(gd),2);em=ej;es=ej;
for k=1:numel(gd)
 x=[gg(k) 1];J=[0 p.omega;-dd.totalPg(k)/p.mass,(-(1-dd.P(k))-p.D)/p.mass];ej(k,:)=sort(real(eig(J)))';
 r=0;if gg(k)>=cfg.R1(1,1)&&gg(k)<=cfg.R1(1,2),r=1;elseif gg(k)>cfg.R2(1,1)&&gg(k)<=cfg.R2(1,2),r=2;end
 if r>0
  side=coupled(gg(k),1,p,n,r);J=[0 p.omega;-side.totalPg/p.mass,(-(1-side.P)-p.D)/p.mass];[M,MD]=metric_value(x,r,cfg,p,side.P);R=chol(M);em(k,:)=sort(eig(M))';es(k,:)=sort(eig((R'\(J'*M+M*J+MD)/R)/2))';
 end
end
fig=figure('Visible','off','Color','w','Position',[50 30 1200 1100]);tl=tiledlayout(fig,3,1,'TileSpacing','compact');vals={ej,em,es};labs={'Real eigenvalues of J (1/s)','Eigenvalues of M','Eigenvalues of Sym(F) (1/s)'};
for j=1:3
 ax=nexttile(tl);hold(ax,'on');left=gd<=rad2deg(cfg.gammaSwitch);right=~left;
 h=plot(ax,gd(left),vals{j}(left,:),'LineWidth',1.4);for q=1:2,plot(ax,gd(right),vals{j}(right,q),'Color',h(q).Color,'LineWidth',1.4,'HandleVisibility','off');end
 yline(ax,0,'k:','HandleVisibility','off');bounds_vertical(ax,cfg);grid(ax,'on');xlabel(ax,'gamma (deg), nu = 1, fixed Pstar = 1');ylabel(ax,labs{j});xlim(ax,[0 180]);
 if j==3,ylim(ax,[-12 8]);end
 for k=1:5
  xx=cfg.gamma0(k);[~,ix]=min(abs(gd-xx));if all(isfinite(vals{j}(ix,:))),plot(ax,[xx xx],vals{j}(ix,:),'o','Color',colors(k,:),'MarkerFaceColor',colors(k,:),'HandleVisibility','off');end
 end
 if j==1
  hh=gobjects(1,5);for k=1:5,hh(k)=plot(ax,NaN,NaN,'o','Color',colors(k,:),'MarkerFaceColor',colors(k,:));end
  legend(ax,[h;hh'],['eig 1','eig 2',labels],'Location','northeast','NumColumns',2,'FontSize',8);
 end
end
title(tl,{'Figure 3 - ConstantRatioLimiter | SCR5/XR10, Xeq=0.350, Pstar=1','Fixed-Pstar state slice; only B is equilibrium. M and Sym(F) shown only on validated branches.'},'FontSize',13);
saveall(fig,out,'Figure3_ConstantRatioLimiter');
writetable(table(gd',ej(:,1),ej(:,2),em(:,1),em(:,2),es(:,1),es(:,2),'VariableNames',{'gamma_deg','eigJ1_real','eigJ2_real','eigM1','eigM2','eigSymF1','eigSymF2'}),fullfile(out,'Figure3_matrix_slice.csv'));
%% Figure 4: phase plane above metric displacement, one column per case.
for s=1:2
 fig=figure('Visible','off','Color','w','Position',[40 40 1800 850]);tl=tiledlayout(fig,2,5,'TileSpacing','compact');
 for k=1:5
  tr=runs{k,s};T=tr.data;ax=nexttile(tl,k);hold(ax,'on');
  plot(ax,T.gamma_deg_unwrapped,T.nu-1,'Color',colors(k,:));
  plot(ax,rad2deg(cfg.gammaEq),0,'ro','MarkerSize',10,'LineWidth',1.7,'DisplayName','Nominal equilibrium');plot(ax,cfg.gamma0(k),0,'ks','MarkerFaceColor','w','DisplayName','Initial point');
  if all(isfinite(tr.clearState)),plot(ax,rad2deg(tr.clearState(1)),tr.clearState(2)-1,'^','Color',[.8 .5 0],'MarkerFaceColor',[1 .75 0]);end
  if tr.valid,plot(ax,T.gamma_deg_unwrapped(end),T.nu(end)-1,'gd','MarkerFaceColor','g');else,plot(ax,T.gamma_deg_unwrapped(end),T.nu(end)-1,'mx','MarkerSize',9,'LineWidth',2);end
  b1=cfg.R1;b2=cfg.R2;rectangle(ax,'Position',[rad2deg(b1(1,1)),b1(2,1)-1,rad2deg(diff(b1(1,:))),diff(b1(2,:))],'EdgeColor',[.3 .5 .9],'LineStyle',':');rectangle(ax,'Position',[rad2deg(b2(1,1)),b2(2,1)-1,rad2deg(diff(b2(1,:))),diff(b2(2,:))],'EdgeColor',[.9 .6 .3],'LineStyle',':');
  xline(ax,rad2deg(cfg.gammaSwitch),'b--','HandleVisibility','off');grid(ax,'on');xlabel(ax,'gamma unwrapped (deg)');ylabel(ax,'nu-1');ax.YAxis.Exponent=0;ytickformat(ax,'%.3f');title(ax,{char(labels(k)),status_title(tr)},'FontSize',9);
  ax=nexttile(tl,5+k);hold(ax,'on');plot(ax,T.time_s,T.L_inactive_chord,'-','Color',colors(k,:));plot(ax,T.time_s,T.L_active_chord,'--','Color',colors(k,:));
  mt=tr.metricTransitions;if ~isempty(mt),plot(ax,mt.time_s,mt.L_active,'v','Color',colors(k,:),'MarkerFaceColor','w');plot(ax,mt.time_s,mt.L_inactive,'o','Color',colors(k,:),'MarkerFaceColor','w');end
  if ~any(isfinite(T.L_inactive_chord)|isfinite(T.L_active_chord)),text(ax,.5,.5,'No validated post-fault metric samples','Units','normalized','HorizontalAlignment','center','FontSize',8);end
  grid(ax,'on');xlabel(ax,'Physical time (s)');ylabel(ax,'L_i: local weighted displacement');xlim(ax,[0 8]);fault_mark(ax,cfg);
 end
 title(tl,{sprintf('Figure 4 - ConstantRatioLimiter | %d%% source sag, 100 ms',round(100*(1-cfg.retained(s)))),'Top: red equilibrium, black initial, gold clearing, green final / magenta validity stop','Bottom: solid M inactive; dashed M active; breaks at switches/outside domains; not one global distance'},'FontSize',12);
 saveall(fig,out,sprintf('Figure4_ConstantRatioLimiter_Sag%d',round(100*(1-cfg.retained(s)))));
end
%% Figure 5A: original bar-map methodology, now initial states (not equilibria).
fig=figure('Visible','off','Color','w','Position',[40 50 1300 650]);ax=axes(fig);hold(ax,'on');
plot(ax,rad2deg(cfg.R1(1,:)),[1 1],'LineWidth',12,'DisplayName','R1: inactive');plot(ax,rad2deg(cfg.R2(1,:)),[2 2],'LineWidth',12,'DisplayName','R2: active');
for k=1:5,plot(ax,cfg.gamma0(k),0,'o','Color',colors(k,:),'MarkerFaceColor',colors(k,:),'HandleVisibility','off');text(ax,cfg.gamma0(k),-.1,char(labels(k)),'Rotation',35,'FontSize',10,'Color',colors(k,:));end
xline(ax,rad2deg(cfg.gammaSwitch),'b--','Switch','HandleVisibility','off');xlim(ax,[10 42]);ylim(ax,[-1 2.5]);grid(ax,'on');xlabel(ax,'Initial gamma0 (deg), nu0=1, Pstar=1');ylabel(ax,'Certificate region');legend(ax,'Location','northeast');
title(ax,{'Figure 5A - ConstantRatioLimiter | initial-state locations','Angular projection only: full nu bounds are required; only Case B is equilibrium'},'FontSize',13);saveall(fig,out,'Figure5A_ConstantRatioLimiter');
%% Figure 5B: original wide top + four lower panels and early-time inset.
for s=1:2
 fig=figure('Visible','off','Color','w','Position',[30 30 1450 1100]);tl=tiledlayout(fig,3,2,'TileSpacing','compact');ax=nexttile(tl,[1 2]);hold(ax,'on');
 for k=1:5,T=runs{k,s}.data;plot(ax,T.time_s,T.gamma_deg_unwrapped,'Color',colors(k,:),'DisplayName',char(labels(k)));end
 gamma_guides(ax,cfg);grid(ax,'on');xlim(ax,[0 8]);xlabel(ax,'Physical time (s)');ylabel(ax,'gamma unwrapped (deg)');legend(ax,'Location','northwest','NumColumns',3,'FontSize',8);title(ax,'Five initial angles; same fixed Pstar=1');fault_mark(ax,cfg);
 selected=[1 3 4 5];
 for j=1:4
  k=selected(j);tr=runs{k,s};T=tr.data;a=nexttile(tl,2+j);hold(a,'on');plot(a,T.time_s,T.gamma_deg_unwrapped,'Color',colors(k,:));yline(a,rad2deg(cfg.gammaEq),'r--','Nominal','LabelHorizontalAlignment','right');
  gamma_guides(a,cfg);fault_mark(a,cfg);grid(a,'on');xlim(a,[0 8]);xlabel(a,'Physical time (s)');ylabel(a,'gamma unwrapped (deg)');title(a,{char(labels(k)),status_title(tr)},'FontSize',10);
 end
 title(tl,{sprintf('Figure 5B - ConstantRatioLimiter | %d%% source sag (retained %.2f)',round(100*(1-cfg.retained(s))),cfg.retained(s)),'Blue dashed: nominal-voltage domain gamma bounds; valid only with branch and nu conditions'},'FontSize',13);
 drawnow;pos=ax.Position;ins=axes(fig,'Position',[pos(1)+.63*pos(3),pos(2)+.15*pos(4),.32*pos(3),.45*pos(4)]);hold(ins,'on');
 for k=1:5,T=runs{k,s}.data;plot(ins,T.time_s,T.gamma_deg_unwrapped,'Color',colors(k,:));end
 gamma_guides(ins,cfg);xlim(ins,[0 .6]);grid(ins,'on');title(ins,'Early time: 0-0.6 s','FontSize',8);ins.FontSize=7;
 saveall(fig,out,sprintf('Figure5B_ConstantRatioLimiter_Sag%d',round(100*(1-cfg.retained(s)))));
end
%% Separate diagnostics: do not replace any baseline figure with these.
fig=figure('Visible','off','Color','w','Position',[60 40 1250 800]);tl=tiledlayout(fig,2,1,'TileSpacing','compact');
for s=1:2,ax=nexttile(tl);hold(ax,'on');for k=1:5,T=runs{k,s}.data;plot(ax,T.time_s,T.frequencyDeviation_Hz,'Color',colors(k,:),'DisplayName',char(labels(k)));end;grid(ax,'on');fault_mark(ax,cfg);xlabel(ax,'Physical time (s)');ylabel(ax,'50(nu-1) (Hz)');xlim(ax,[0 8]);title(ax,sprintf('%d%% source sag',round(100*(1-cfg.retained(s)))));legend(ax,'Location','northeast','NumColumns',3);end
title(tl,'Supplement - Frequency response, same five case colours');saveall(fig,out,'Supplement_Frequency_ConstantRatioLimiter');
for s=1:2
 fig=figure('Visible','off','Color','w','Position',[40 30 1800 1000]);tl=tiledlayout(fig,3,5,'TileSpacing','compact');
 for k=1:5
  T=runs{k,s}.data;
  ax=nexttile(tl,k);plot(ax,T.time_s,T.rawCurrent_pu,'--','Color',colors(k,:));hold(ax,'on');plot(ax,T.time_s,T.limitedCurrent_pu,'-','Color',colors(k,:));yline(ax,1.2,'r:');grid(ax,'on');ylabel(ax,'Current (pu)');title(ax,labels(k));fault_mark(ax,cfg);xlim(ax,[0 8]);if k==1,legend(ax,'Raw reference','Limited','Imax','FontSize',8);end
  ax=nexttile(tl,5+k);plot(ax,T.time_s,T.requiredVoltage_pu,'Color',colors(k,:));hold(ax,'on');yline(ax,1.3,'r--');grid(ax,'on');ylabel(ax,'Required command (pu)');fault_mark(ax,cfg);xlim(ax,[0 8]);
  ax=nexttile(tl,10+k);stairs(ax,T.time_s,double(T.limiterActive),'Color',colors(k,:));grid(ax,'on');ylabel(ax,'Limiter active');ylim(ax,[-.1 1.15]);yticks(ax,[0 1]);xlabel(ax,'Physical time (s)');fault_mark(ax,cfg);xlim(ax,[0 8]);
 end
 title(tl,sprintf('Supplement - Current and voltage validity | %d%% source sag',round(100*(1-cfg.retained(s)))));saveall(fig,out,sprintf('Supplement_Current_ConstantRatioLimiter_Sag%d',round(100*(1-cfg.retained(s)))));
end
end

function [M,MD]=metric_value(x,branch,cfg,p,Pe)
m=cfg.mi;if branch==2,m=cfg.ma;end;z=(x-m.center')./m.scale';f=[p.omega*(x(2)-1),((1-Pe)/x(2)-p.D*(x(2)-1))/p.mass];dz=f./m.scale';b=[1 z z(1)^2 z(1)*z(2) z(2)^2];bd=[0 dz 2*z(1)*dz(1) z(2)*dz(1)+z(1)*dz(2) 2*z(2)*dz(2)];M=zeros(2);MD=M;for k=1:6,M=M+b(k)*m.C(:,:,k);MD=MD+bd(k)*m.C(:,:,k);end
end
function gamma_guides(ax,cfg)
for g=unique([cfg.R1(1,:) cfg.R2(1,:)]),yline(ax,rad2deg(g),'b--','HandleVisibility','off','Alpha',.35);end
end
function bounds_vertical(ax,cfg)
for g=unique([cfg.R1(1,:) cfg.R2(1,:)]),xline(ax,rad2deg(g),'b--','HandleVisibility','off','Alpha',.35);end
end
function fault_mark(ax,cfg)
xline(ax,cfg.tClear,'k:','HandleVisibility','off');yl=ylim(ax);h=patch(ax,[0 cfg.tClear cfg.tClear 0],[yl(1) yl(1) yl(2) yl(2)],[1 .88 .55],'FaceAlpha',.18,'EdgeColor','none','HandleVisibility','off');uistack(h,'bottom');ylim(ax,yl);
end
function s=status_title(tr)
if ~tr.valid,s=sprintf('Validity stops at %.5f s',tr.t(end));else
 ix=tr.t>=7.2;s='No recovery within horizon';if all(abs(tr.data.gamma_deg_unwrapped(ix)-20.633904875622)<.1 & abs(tr.x(ix,2)-1)<1e-4),s='Recovered numerically';end
end
end
function saveall(fig,out,name)
drawnow;exportgraphics(fig,fullfile(out,[name '.png']),'Resolution',300);exportgraphics(fig,fullfile(out,[name '.pdf']),'ContentType','vector');savefig(fig,fullfile(out,[name '.fig']));close(fig);
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






