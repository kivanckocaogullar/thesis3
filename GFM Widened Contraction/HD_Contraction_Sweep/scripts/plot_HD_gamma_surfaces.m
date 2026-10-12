function plot_HD_gamma_surfaces(root,tag)
%PLOT_HD_GAMMA_SURFACES Figure A: Study A angle bounds versus H and D (nu band fixed).
% Surfaces are drawn only through computed points; infeasible/inconclusive points are
% left as holes (NaN) and marked on the floor. No interpolation or smoothing.
s=hd_style();a=load(fullfile(root,'results','angle_bounds',['angle_bounds_' tag '.mat']));
[X,Y]=ndgrid(a.Hgrid,a.Dgrid);G=a.G;SA=a.SA;sw=rad2deg(a.dom.gamma_switch);out=fullfile(root,'figures','angle_surfaces');
gi=rad2deg(a.cfg.gammaFixedInactiveHi);ga=rad2deg(a.cfg.gammaFixedActiveLo);
dense=numel(a.Dgrid)>8;ms=22;if dense,ms=9;end
cases={'combined',[1 2];'inactive',1;'active',2};
for c=1:size(cases,1)
 br=cases{c,2};fig=figure('Visible','off','Color','w','Units','centimeters','Position',[2 2 24 13]);ax=axes(fig);hold(ax,'on');h=gobjects(0);lab={};
 zl=[];
 if any(br==1)
  Z=G.min_inactive;h(end+1)=surf(ax,X,Y,Z,'FaceColor',s.inactive,'FaceAlpha',.80,'EdgeColor',edgec(s.inactive*.6,dense),'LineWidth',.7);lab{end+1}='\gamma_{min} (limiter inactive), searched';
  hlines(ax,X,Y,Z,s.inactive*.6,dense);scatter3(ax,X(:),Y(:),Z(:),ms,s.inactive*.6,'filled','HandleVisibility','off');zl=[zl;Z(:)];
  if c==2,h(end+1)=fixed_plane(ax,X,Y,gi,s);lab{end+1}=sprintf('\\gamma_{max} (inactive) = %.2f° — fixed by limiter switch',gi);zl=[zl;gi];end
 end
 if any(br==2)
  Z=G.max_active;h(end+1)=surf(ax,X,Y,Z,'FaceColor',s.active,'FaceAlpha',.80,'EdgeColor',edgec(s.active*.6,dense),'LineWidth',.7);lab{end+1}='\gamma_{max} (limiter active), searched';
  hlines(ax,X,Y,Z,s.active*.6,dense);scatter3(ax,X(:),Y(:),Z(:),ms,s.active*.6,'filled','HandleVisibility','off');zl=[zl;Z(:)];
  if c==3,h(end+1)=fixed_plane(ax,X,Y,ga,s);lab{end+1}=sprintf('\\gamma_{min} (active) = %.2f° — fixed by limiter switch',ga);zl=[zl;ga];end
 end
 if c==1,h(end+1)=fixed_plane(ax,X,Y,sw,s);lab{end+1}=sprintf('Limiter switch \\gamma_{sw} = %.2f° (fixed; edges %.2f° / %.2f°)',sw,gi,ga);zl=[zl;sw];end
 zlo=floor(min(zl,[],'omitnan')/5)*5-5;zhi=ceil(max(zl,[],'omitnan')/5)*5+5;
 % Baseline (H = 3.4 s, D = 25) on each searched surface.
 [~,ib]=min(abs(a.Hgrid-3.4));[~,kb]=min(abs(a.Dgrid-25));
 for b=br
  if b==1,zb=G.min_inactive(ib,kb);else,zb=G.max_active(ib,kb);end
  if ~isnan(zb),hb=plot3(ax,X(ib,kb),Y(ib,kb),zb,'kp','MarkerSize',12,'MarkerFaceColor','k');end
 end
 if exist('hb','var'),h(end+1)=hb;lab{end+1}='Baseline (H = 3.4 s, D = 25)';clear hb,end
 % Points without a feasible rectangle, marked on the floor.
 [hx,lx]=floor_markers(ax,X,Y,SA,br,zlo,s);h=[h hx];lab=[lab lx];
 hold(ax,'off');style_axes(ax,s);view(ax,[-38 24]);set(ax,'Projection','orthographic');zlim(ax,[zlo zhi]);
 xlabel(ax,'Inertia {\itH} (s)');ylabel(ax,'Damping {\itD} = {\itK} + {\itK}_d');zlabel(ax,'Angle \gamma (deg)');
 xticks(ax,a.Hgrid);if numel(a.Dgrid)<=8,yticks(ax,a.Dgrid);end;xlim(ax,[min(a.Hgrid) max(a.Hgrid)]);ylim(ax,[min(a.Dgrid) max(a.Dgrid)]);
 ttl={'Contraction angle bounds versus inertia and damping',sprintf('\\nu \\in [%.2f, %.2f] fixed;  \\lambda_{inactive} = %g s^{-1},  \\lambda_{active} = %g s^{-1};  K = %g',a.cfg.nuBandA,a.cfg.lambda,a.cfg.K)};
 title(ax,ttl,'FontSize',s.fsTitle,'FontWeight','normal');
 legend(ax,h,lab,'Location','eastoutside','FontSize',s.fsSmall,'Box','off');
 annotation(fig,'textbox',[0 0 1 .03],'String',s.note,'EdgeColor','none','HorizontalAlignment','right','FontSize',8,'Color',[.35 .35 .35],'FontName',s.font);
 hd_export_figure(fig,fullfile(out,sprintf('gamma_bounds_%s_%s',cases{c,1},tag)));close(fig);
end
T=table(X(:),Y(:),Y(:)-a.cfg.K,G.min_inactive(:),G.max_inactive(:),G.min_active(:),G.max_active(:),reshape(SA(:,:,1),[],1),reshape(SA(:,:,2),[],1), ...
 'VariableNames',{'H_s','D','Kd','gamma_min_inactive_deg','gamma_max_inactive_deg_fixed','gamma_min_active_deg_fixed','gamma_max_active_deg','status_inactive','status_active'});
writetable(T,fullfile(out,['gamma_surface_data_' tag '.csv']));
end

function c=edgec(c,dense),if dense,c='none';end,end
function hlines(ax,X,Y,Z,c,dense)
% Dense D grids: draw one line per H instead of the full mesh.
if dense,for i=1:size(X,1),plot3(ax,X(i,:),Y(i,:),Z(i,:),'-','Color',c,'LineWidth',.9,'HandleVisibility','off');end,end
end
function h=fixed_plane(ax,X,Y,z,s)
Z=z+0*X;h=surf(ax,X,Y,Z,'FaceColor',s.fixed,'FaceAlpha',.12,'EdgeColor','none');
xe=[min(X(:)) max(X(:)) max(X(:)) min(X(:)) min(X(:))];ye=[min(Y(:)) min(Y(:)) max(Y(:)) max(Y(:)) min(Y(:))];
plot3(ax,xe,ye,z+0*xe,'--','Color',s.fixed,'LineWidth',1,'HandleVisibility','off');
end
function [h,lab]=floor_markers(ax,X,Y,S,br,zlo,s)
h=gobjects(0);lab={};
nf=false(size(X));ic=nf;for b=br,nf=nf|S(:,:,b)=="no_feasible_rectangle";ic=ic|S(:,:,b)=="inconclusive";end
if any(nf(:)),h(end+1)=plot3(ax,X(nf),Y(nf),zlo+0*X(nf),'kx','MarkerSize',9,'LineWidth',1.5);lab{end+1}='No feasible rectangle found (marked on floor)';end
if any(ic(:)),h(end+1)=plot3(ax,X(ic),Y(ic),zlo+0*X(ic),'^','MarkerSize',7,'Color',s.inconclusive,'MarkerFaceColor',s.inconclusive);lab{end+1}='Inconclusive solver outcome';end
end
function style_axes(ax,s)
set(ax,'FontName',s.font,'FontSize',s.fs,'LineWidth',.8,'TickDir','out','Box','on','BoxStyle','back','GridAlpha',.15);grid(ax,'on');
end
