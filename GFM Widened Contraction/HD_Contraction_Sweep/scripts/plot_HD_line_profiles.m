function plot_HD_line_profiles(root,tag)
%PLOT_HD_LINE_PROFILES Bounds versus D, one line per H (suited to dense D grids).
% (a) active gamma_max (Study A); (b) active symmetric half-width delta* (Study B);
% (c) inactive symmetric half-width delta*; (d) inactive feasibility versus D = 4*lambda*H.
% Open markers: edge reached the search-domain limit. Grey crosses: an inconclusive
% probe lay inside the final bisection bracket (bound possibly conservative).
% Raw computed values only; no smoothing or interpolation.
s=hd_style();r=load(fullfile(root,'results','diagnostics',['sweep_results_' tag '.mat']),'res','cfg','dom');
res=r.res;cfg=r.cfg;H=cfg.H;D=cfg.D;nH=numel(H);out=fullfile(root,'figures','region_comparisons');
cm=[0.267 0.005 0.329;0.229 0.322 0.546;0.128 0.567 0.551;0.369 0.789 0.383];if nH>4,cm=parula(nH+1);cm=cm(1:nH,:);end
fig=figure('Visible','off','Color','w','Units','centimeters','Position',[2 2 26 18]);tl=tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
pan={'(a) Limiter active: \gamma_{max}  (\nu \in [0.98, 1.02])','(b) Limiter active: symmetric band 1 \pm \delta*  (\gamma \in [24.42°, 43°])', ...
     '(c) Limiter inactive: symmetric band 1 \pm \delta*  (\gamma \in [10°, 24.32°])','(d) Limiter inactive: feasibility at \lambda = 1 s^{-1}'};
ylab={'\gamma_{max} (deg)','\delta* (pu)','\delta* (pu)','Inertia {\itH} (s)'};
hl=gobjects(nH,1);
for q=1:3
 ax=nexttile(tl);hold(ax,'on');
 for i=1:nH
  v=nan(1,numel(D));lim=false(1,numel(D));inc=lim;
  for k=1:numel(D)
   x=res(i,k);
   switch q
    case 1,S=x.A.active;if S.status=="feasible",v(k)=rad2deg(S.bound);end
    case 2,S=x.B.active.sym;if S.status=="feasible",v(k)=S.bound;end
    case 3,S=x.B.inactive.sym;if S.status=="feasible",v(k)=S.bound;end
   end
   if S.status=="feasible",lim(k)=S.limited_by=="search_domain_limit";inc(k)=S.inconclusive_in_bracket;end
  end
  hl(i)=plot(ax,D,v,'-','Color',cm(i,:),'LineWidth',1.4);
  plot(ax,D(~lim),v(~lim),'o','Color',cm(i,:),'MarkerFaceColor',cm(i,:),'MarkerSize',3.5,'HandleVisibility','off');
  plot(ax,D(lim),v(lim),'o','Color',cm(i,:),'MarkerFaceColor','w','MarkerSize',4.5,'HandleVisibility','off');
  plot(ax,D(inc),v(inc),'x','Color',[.45 .45 .45],'MarkerSize',6,'LineWidth',1,'HandleVisibility','off');
 end
 [~,ib]=min(abs(H-3.4));[~,kb]=min(abs(D-25));x=res(ib,kb);
 switch q
  case 1,yb=rad2deg(x.A.active.bound);case 2,yb=x.B.active.sym.bound;case 3,yb=x.B.inactive.sym.bound;
 end
 plot(ax,D(kb),yb,'kp','MarkerFaceColor','k','MarkerSize',10,'HandleVisibility','off');
 if q>1,yline(ax,min(1-cfg.nuDomain(1),cfg.nuDomain(2)-1),':','Color',s.fixed,'HandleVisibility','off');ylim(ax,[0 .055]);end
 hold(ax,'off');sty(ax,s);xlabel(ax,'Damping {\itD} = {\itK} + {\itK}_d');ylabel(ax,ylab{q});title(ax,pan{q},'FontSize',s.fs,'FontWeight','normal');xlim(ax,[min(D) max(D)]);
end
% (d) inactive feasibility map along D for each H, with the analytic line D = 4*lambda*H.
ax=nexttile(tl);hold(ax,'on');
for i=1:nH
 st=arrayfun(@(k)res(i,k).A.inactive.status,1:numel(D));f=st=="feasible";
 plot(ax,D(f),H(i)+0*D(f),'s','MarkerSize',4.5,'Color',s.feasible,'MarkerFaceColor',s.feasible,'HandleVisibility','off');
 plot(ax,D(~f),H(i)+0*D(~f),'x','MarkerSize',6,'Color',[.35 .35 .35],'LineWidth',1.1,'HandleVisibility','off');
end
dd=linspace(min(D),max(D),50);hn=plot(ax,dd,dd/(4*cfg.lambda(1)),'k--','LineWidth',1.2);
hf=plot(ax,NaN,NaN,'s','Color',s.feasible,'MarkerFaceColor',s.feasible);hx=plot(ax,NaN,NaN,'x','Color',[.35 .35 .35],'LineWidth',1.1);
hold(ax,'off');sty(ax,s);xlabel(ax,'Damping {\itD} = {\itK} + {\itK}_d');ylabel(ax,ylab{4});title(ax,pan{4},'FontSize',s.fs,'FontWeight','normal');
xlim(ax,[min(D)-.5 max(D)+.5]);ylim(ax,[min(H)-.6 max(H)+.6]);yticks(ax,H);
legend(ax,[hf hx hn],{'Feasible rectangle found','No feasible rectangle','D = 4\lambda_{inactive}H (necessary)'},'Location','southeast','FontSize',s.fsSmall,'Box','on','Color','w');
lg=legend(nexttile(tl,1),hl,compose('H = %g s',H(:)),'NumColumns',nH,'FontSize',s.fsSmall,'Box','off');lg.Layout.Tile='south';
title(tl,{'Contraction bounds versus damping for fixed inertia (K = 20, D = 20 + K_d)', ...
 '\fontsize{9}Filled: limited by contraction feasibility; open: search-domain limit; grey \times: inconclusive probe in final bracket; black star: baseline'},'FontSize',s.fsTitle,'FontName',s.font);
annotation(fig,'textbox',[0 0 1 .025],'String',s.note,'EdgeColor','none','HorizontalAlignment','right','FontSize',8,'Color',[.35 .35 .35],'FontName',s.font);
hd_export_figure(fig,fullfile(out,['profiles_vs_D_' tag]));close(fig);
end
function sty(ax,s)
set(ax,'FontName',s.font,'FontSize',s.fs-1,'TickDir','out','Box','on','LineWidth',.8,'GridAlpha',.12);grid(ax,'on');
end
