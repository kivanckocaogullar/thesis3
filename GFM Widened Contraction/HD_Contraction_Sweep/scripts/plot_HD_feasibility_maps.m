function plot_HD_feasibility_maps(root,tag)
%PLOT_HD_FEASIBILITY_MAPS Figure C: H-D feasibility maps per branch and search family.
% Cells span the midpoints between grid values. The dashed line D = 4*lambda*H is the
% analytic necessary condition for any rectangle containing the equilibrium (at the
% equilibrium Mdot = 0, so trace Sym(F) = trace J = -D/(2H) must be <= -2*lambda);
% it applies to the inactive branch, whose rectangles always contain the equilibrium.
s=hd_style();a=load(fullfile(root,'results','angle_bounds',['angle_bounds_' tag '.mat']));f=load(fullfile(root,'results','frequency_bounds',['frequency_bounds_' tag '.mat']));
out=fullfile(root,'figures','feasibility_maps');H=a.Hgrid;D=a.Dgrid;he=edges(H);de=edges(D);
fig=figure('Visible','off','Color','w','Units','centimeters','Position',[2 2 22 18]);tl=tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
studies={'A','B'};branch={'Limiter inactive','Limiter active'};
for st=1:2
 for b=1:2
  ax=nexttile(tl);hold(ax,'on');
  if st==1,S=a.SA(:,:,b);else,S=f.SB(:,:,b);end
  for i=1:numel(H)
   for k=1:numel(D)
    switch S(i,k)
     case "feasible",c=s.feasible;case "no_feasible_rectangle",c=s.nofeasible;case "inconclusive",c=s.inconclusive;otherwise,c=s.noteval;
    end
    if numel(D)<=8,patch(ax,he([i i+1 i+1 i]),de([k k k+1 k+1]),c,'EdgeColor','w','LineWidth',1.2);
    else,patch(ax,he([i i+1 i+1 i]),de([k k k+1 k+1]),c,'EdgeColor',c,'LineWidth',.1);end
    txt=cell_text(a,f,st,b,i,k);if numel(D)<=8&&~isempty(txt),text(ax,H(i),D(k),txt,'HorizontalAlignment','center','FontSize',7.5,'FontName',s.font);end
   end
  end
  if numel(D)>8,for i=2:numel(H),xline(ax,he(i),'w','LineWidth',1.5,'HandleVisibility','off');end,end
  [~,ib]=min(abs(H-3.4));[~,kb]=min(abs(D-25));
  patch(ax,he([ib ib+1 ib+1 ib]),de([kb kb kb+1 kb+1]),'k','FaceColor','none','EdgeColor','k','LineWidth',2);
  if b==1
   hh=linspace(he(1),he(end),100);plot(ax,hh,4*a.cfg.lambda(1)*hh,'k--','LineWidth',1.2);
  end
  hold(ax,'off');axis(ax,[he(1) he(end) de(1) de(end)]);set(ax,'FontName',s.font,'FontSize',s.fs-1,'TickDir','out','Box','on','Layer','top','LineWidth',.8);
  xticks(ax,H);if numel(D)<=8,yticks(ax,D);end;xlabel(ax,'Inertia {\itH} (s)');ylabel(ax,'Damping {\itD} = {\itK} + {\itK}_d');
  if st==1,ttl=sprintf('Study A (\\gamma searched, \\nu \\in [%.2f, %.2f])',a.cfg.nuBandA);else,ttl='Study B (\nu searched, \gamma fixed)';end
  title(ax,{branch{b},ttl},'FontSize',s.fs,'FontWeight','normal');
 end
end
% Legend proxies.
ax=nexttile(tl,1);hold(ax,'on');hp=gobjects(5,1);
hp(1)=patch(ax,NaN,NaN,s.feasible);hp(2)=patch(ax,NaN,NaN,s.nofeasible);hp(3)=patch(ax,NaN,NaN,s.inconclusive);
hp(4)=patch(ax,NaN,NaN,'k','FaceColor','none','LineWidth',2);hp(5)=plot(ax,NaN,NaN,'k--','LineWidth',1.2);hold(ax,'off');
lg=legend(ax,hp,{'Numerically feasible rectangle','No feasible rectangle found','Inconclusive solver outcome','Baseline (H = 3.4 s, D = 25)', ...
 sprintf('D = 4\\lambda_{inactive}H: necessary condition (inactive)')},'NumColumns',3,'FontSize',s.fsSmall,'Box','off');lg.Layout.Tile='south';
ttl={'Contraction-rectangle feasibility in the {\itH}–{\itD} plane',sprintf('\\lambda_{inactive} = %g s^{-1},  \\lambda_{active} = %g s^{-1},  M < 10I,  K = %g',a.cfg.lambda,a.cfg.K)};
if numel(D)<=8,ttl{end+1}='\fontsize{9}Study B cells: [\nu_{min}, \nu_{max}];  * one-sided extents not jointly feasible (combined test inconclusive), jointly feasible band 1 \pm \delta* given below';end
title(tl,ttl,'FontSize',s.fsTitle,'FontName',s.font);
annotation(fig,'textbox',[0 0 1 .025],'String',s.note,'EdgeColor','none','HorizontalAlignment','right','FontSize',8,'Color',[.35 .35 .35],'FontName',s.font);
hd_export_figure(fig,fullfile(out,['feasibility_maps_' tag]));close(fig);
end
function e=edges(v)
v=v(:)';if numel(v)==1,e=v+[-.5 .5];return,end
m=(v(1:end-1)+v(2:end))/2;e=[v(1)-(m(1)-v(1)) m v(end)+(v(end)-m(end))];
end
function t=cell_text(a,f,st,b,i,k)
t='';
if st==1
 if b==1,v=a.G.min_inactive(i,k);if ~isnan(v),t=sprintf('%.1f°',v);end
 else,v=a.G.max_active(i,k);if ~isnan(v),t=sprintf('%.1f°',v);end,end
else
 lab={'inactive','active'};lo=f.U.(['min_' lab{b}])(i,k);hi=f.U.(['max_' lab{b}])(i,k);
 if ~isnan(lo)
  t=sprintf('[%.3f,\n %.3f]',lo,hi);
  if f.U.(['combined_' lab{b}])(i,k)~="feasible",t=sprintf('[%.3f,\n %.3f]*\n1 \\pm %.4f',lo,hi,f.U.(['delta_' lab{b}])(i,k));end
 end
end
end
