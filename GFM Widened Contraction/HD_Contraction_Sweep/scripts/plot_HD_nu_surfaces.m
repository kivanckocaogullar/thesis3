function plot_HD_nu_surfaces(root,tag)
%PLOT_HD_NU_SURFACES Figure B: Study B frequency bounds versus H and D (gamma fixed).
% Two variants are exported:
%  'onesided'  : nu_min and nu_max, each a one-sided extension from the largest
%                symmetric band (the other edge held at 1 -/+ delta*);
%  'symmetric' : the largest jointly feasible symmetric band 1 +/- delta*.
% The active branch shows a strong trade-off between the lower and upper edges, so the
% one-sided values are conditional on the held edge and can jump with D; the symmetric
% band is the jointly feasible quantity.
% Filled markers: limited by contraction feasibility. Open markers: search-domain limit
% reached. Holes: no feasible rectangle. No interpolation or smoothing.
s=hd_style();a=load(fullfile(root,'results','frequency_bounds',['frequency_bounds_' tag '.mat']));
[X,Y]=ndgrid(a.Hgrid,a.Dgrid);U=a.U;SB=a.SB;out=fullfile(root,'figures','frequency_surfaces');lab={'inactive','active'};
names={'Limiter inactive','Limiter active'};col={s.inactive,s.active};dense=numel(a.Dgrid)>8;
gB=cellfun(@rad2deg,a.cfg.gammaB,'UniformOutput',false);dlim=min(1-a.cfg.nuDomain(1),a.cfg.nuDomain(2)-1);
for variant=["onesided","symmetric"]
 fig=figure('Visible','off','Color','w','Units','centimeters','Position',[2 2 26 13.5]);tl=tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
 for b=1:2
  ax=nexttile(tl);hold(ax,'on');c=col{b};h=gobjects(0);L={};dl=U.(['delta_' lab{b}]);
  if variant=="onesided"
   Zl=U.(['min_' lab{b}]);Zu=U.(['max_' lab{b}]);limL=U.(['limited_min_' lab{b}])=="search_domain_limit";limU=U.(['limited_max_' lab{b}])=="search_domain_limit";
   nm={'\nu_{max} (one-sided)','\nu_{min} (one-sided)'};
  else
   Zl=1-dl;Zu=1+dl;limL=dl>=dlim-1e-12;limU=limL;nm={'1 + \delta*','1 - \delta*'};
  end
  ec=c*.6;if dense,ec='none';end
  h(end+1)=surf(ax,X,Y,Zu,'FaceColor',c,'FaceAlpha',.80,'EdgeColor',ec,'LineWidth',.7);L{end+1}=nm{1};
  h(end+1)=surf(ax,X,Y,Zl,'FaceColor',c*.55+.45,'FaceAlpha',.80,'EdgeColor',ec,'LineWidth',.7);L{end+1}=nm{2};
  if dense,for i=1:numel(a.Hgrid),plot3(ax,X(i,:),Y(i,:),Zu(i,:),'-','Color',c*.6,'LineWidth',.8,'HandleVisibility','off');plot3(ax,X(i,:),Y(i,:),Zl(i,:),'-','Color',c*.6,'LineWidth',.8,'HandleVisibility','off');end,end
  h(end+1)=surf(ax,X,Y,1+0*X,'FaceColor',s.fixed,'FaceAlpha',.10,'EdgeColor','none');L{end+1}='\nu = 1 (nominal)';
  ms=24;if dense,ms=10;end
  for e=1:2
   if e==1,Z=Zl;lim=limL;else,Z=Zu;lim=limU;end;ok=~isnan(Z);
   scatter3(ax,X(ok&~lim),Y(ok&~lim),Z(ok&~lim),ms,c*.6,'filled','HandleVisibility','off');
   if any(ok(:)&lim(:)),hs=scatter3(ax,X(ok&lim),Y(ok&lim),Z(ok&lim),ms*1.25,c*.6,'LineWidth',1);end
  end
  if exist('hs','var'),h(end+1)=hs;L{end+1}=sprintf('Search-domain limit reached (\\nu \\in [%.2f, %.2f])',a.cfg.nuDomain);clear hs,end
  [~,ib]=min(abs(a.Hgrid-3.4));[~,kb]=min(abs(a.Dgrid-25));
  if ~isnan(Zu(ib,kb)),h(end+1)=plot3(ax,X(ib,kb)*[1 1],Y(ib,kb)*[1 1],[Zl(ib,kb) Zu(ib,kb)],'kp','MarkerSize',11,'MarkerFaceColor','k');L{end+1}='Baseline (H = 3.4 s, D = 25)';end
  nf=SB(:,:,b)=="no_feasible_rectangle";ic=SB(:,:,b)=="inconclusive";zf=a.cfg.nuDomain(1)-.005;
  if any(nf(:)),h(end+1)=plot3(ax,X(nf),Y(nf),zf+0*X(nf),'kx','MarkerSize',7,'LineWidth',1.2);L{end+1}='No feasible rectangle found';end
  if any(ic(:)),h(end+1)=plot3(ax,X(ic),Y(ic),zf+0*X(ic),'^','MarkerSize',7,'Color',s.inconclusive,'MarkerFaceColor',s.inconclusive);L{end+1}='Inconclusive';end
  hold(ax,'off');set(ax,'FontName',s.font,'FontSize',s.fs,'LineWidth',.8,'TickDir','out','Box','on','BoxStyle','back','GridAlpha',.15,'Projection','orthographic');grid(ax,'on');view(ax,[-38 22]);
  zlim(ax,[a.cfg.nuDomain(1)-.005 a.cfg.nuDomain(2)+.005]);xticks(ax,a.Hgrid);if ~dense,yticks(ax,a.Dgrid);end;xlim(ax,[min(a.Hgrid) max(a.Hgrid)]);ylim(ax,[min(a.Dgrid) max(a.Dgrid)]);
  xlabel(ax,'Inertia {\itH} (s)');ylabel(ax,'Damping {\itD} = {\itK} + {\itK}_d');zlabel(ax,'Normalized frequency \nu (pu)');
  title(ax,{names{b},sprintf('\\gamma \\in [%.2f°, %.2f°] fixed,  \\lambda = %g s^{-1}',gB{b},a.cfg.lambda(b))},'FontSize',s.fs,'FontWeight','normal');
  legend(ax,h,L,'Location','southoutside','FontSize',s.fsSmall,'Box','off');
 end
 if variant=="onesided"
  ttl='Contraction frequency bounds versus inertia and damping (one-sided extensions)';nt=';  \nu_{min}, \nu_{max}: one-sided extensions from the largest symmetric band';fn=['nu_bounds_' tag];
 else
  ttl='Contraction frequency bounds versus inertia and damping (largest symmetric band)';nt=';  jointly feasible band 1 \pm \delta*';fn=['nu_bounds_symmetric_' tag];
 end
 title(tl,ttl,'FontSize',s.fsTitle,'FontName',s.font);
 annotation(fig,'textbox',[0 0 1 .04],'String',[s.note nt],'VerticalAlignment','bottom','EdgeColor','none','HorizontalAlignment','right','FontSize',8,'Color',[.35 .35 .35],'FontName',s.font);
 hd_export_figure(fig,fullfile(out,fn));close(fig);
end
T=table(X(:),Y(:),Y(:)-a.cfg.K,U.min_inactive(:),U.max_inactive(:),U.delta_inactive(:),U.combined_inactive(:),U.min_active(:),U.max_active(:),U.delta_active(:),U.combined_active(:), ...
 'VariableNames',{'H_s','D','Kd','nu_min_inactive','nu_max_inactive','sym_halfwidth_inactive','combined_status_inactive','nu_min_active','nu_max_active','sym_halfwidth_active','combined_status_active'});
writetable(T,fullfile(out,['nu_surface_data_' tag '.csv']));
end
