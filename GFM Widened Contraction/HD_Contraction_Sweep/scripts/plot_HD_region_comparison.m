function plot_HD_region_comparison(root,tag)
%PLOT_HD_REGION_COMPARISON Figure D: gamma-nu rectangles for selected (H, D) versus baseline.
% Filled: Study A rectangle (gamma searched, nu band fixed). Solid outline: Study B
% rectangle (gamma fixed, nu searched; combined [nu_min, nu_max] when jointly feasible,
% otherwise the symmetric band). Dashed black: baseline interval-verified rectangles.
s=hd_style();r=load(fullfile(root,'results','diagnostics',['sweep_results_' tag '.mat']),'res','cfg','dom');
res=r.res;cfg=r.cfg;dom=r.dom;H=cfg.H;D=cfg.D;out=fullfile(root,'figures','region_comparisons');
pickD=unique(arrayfun(@(d)nearest(D,d),[20 25 35 50]));pickH=unique(arrayfun(@(h)nearest(H,h),[1 3.4 5.5 8]));
sets={'D_effect',nearest(H,3.4)+0*pickD,pickD,sprintf('Effect of damping at H = %g s',nearest(H,3.4)); ...
      'H_effect',pickH,nearest(D,50)+0*pickH,sprintf('Effect of inertia at D = %g',nearest(D,50))};
for q=1:2
 hs=sets{q,2};ds=sets{q,3};m=numel(hs);
 fig=figure('Visible','off','Color','w','Units','centimeters','Position',[2 2 6.2*m+2 9.5]);tl=tiledlayout(fig,1,m,'TileSpacing','compact','Padding','compact');
 for j=1:m
  i=find(H==hs(j));k=find(D==ds(j));x=res(i,k);ax=nexttile(tl);hold(ax,'on');h=gobjects(0);L={};
  for b=1:2
   c={s.inactive,s.active};cl={s.inactiveLight,s.activeLight};nm={'inactive','active'};
   A=x.A.(nm{b});
   if A.status=="feasible",h(end+1)=rect(ax,A.best.bounds,cl{b},c{b},'-',.55);L{end+1}=['Study A, ' nm{b}];end %#ok<AGROW>
   B=x.B.(nm{b});
   if B.status=="feasible",R=B.combined;if R.status~="feasible",R=B.sym.best;end
    bb=R.bounds;h(end+1)=plot(ax,rad2deg(bb([1 2 2 1 1],1)),bb([1 1 2 2 1],2),'-','Color',c{b}*.7,'LineWidth',1.6);L{end+1}=['Study B, ' nm{b}];end %#ok<AGROW>
   bb=cfg.baseBounds{b};hb=plot(ax,rad2deg(bb([1 2 2 1 1],1)),bb([1 1 2 2 1],2),'k--','LineWidth',1);
  end
  h(end+1)=hb;L{end+1}='Baseline (interval-verified)';
  xline(ax,rad2deg(dom.gamma_switch),':','Color',s.fixed,'LineWidth',1,'HandleVisibility','off');
  h(end+1)=plot(ax,rad2deg(dom.gamma_eq),1,'kp','MarkerFaceColor','k','MarkerSize',9);L{end+1}='Equilibrium';
  if x.A.inactive.status~="feasible"&&x.B.inactive.status~="feasible"
   text(ax,0,1.045,'inactive: no feasible rectangle','FontSize',8,'Color',s.inactive,'HorizontalAlignment','center','FontName',s.font);
  end
  hold(ax,'off');axis(ax,[-30 70 cfg.nuDomain(1)-.005 cfg.nuDomain(2)+.005]);
  set(ax,'FontName',s.font,'FontSize',s.fs-1,'TickDir','out','Box','on','Layer','top','LineWidth',.8);grid(ax,'on');set(ax,'GridAlpha',.12);
  xlabel(ax,'Angle \gamma (deg)');if j==1,ylabel(ax,'Normalized frequency \nu (pu)');end
  title(ax,sprintf('H = %g s, D = %g (K_d = %g)',x.H,x.D,x.Kd),'FontSize',s.fs,'FontWeight','normal');
  if j==m,hl=h;Ll=L;end
 end
 lg=legend(nexttile(tl,m),hl,Ll,'NumColumns',3,'FontSize',s.fsSmall,'Box','off');lg.Layout.Tile='south';
 title(tl,{sets{q,4},s.note},'FontSize',s.fsTitle,'FontName',s.font);
 hd_export_figure(fig,fullfile(out,sprintf('regions_%s_%s',sets{q,1},tag)));close(fig);
end
end
function v=nearest(g,x),[~,i]=min(abs(g-x));v=g(i);end
function h=rect(ax,bb,fc,ec,ls,fa)
g=rad2deg(bb(:,1));h=patch(ax,g([1 2 2 1]),bb([1 1 2 2],2),fc,'EdgeColor',ec,'LineStyle',ls,'FaceAlpha',fa,'LineWidth',1);
end
