function plot_saved_widened_FCT(out)
% Reuses the previous FCT plotting procedure, with the widened metrics.
load(fullfile(out,'widened_fct_analysis.mat'),'results','runs','cfg','FCT','x0','metrics');
colors=[0 .45 .70;0 .62 .45;.80 .47 .05;.65 .25 .65;.85 .25 .10];
markers={'o','s','d','^','v'};
for kind=1:4
 fig=figure('Visible','off','Color','w','Position',[60 60 1250 850],'Renderer','painters');
 ax=axes(fig,'Position',[.10 .25 .85 .65]);hold(ax,'on');
 set(ax,'FontName','Arial','FontSize',13,'LineWidth',1,'Box','off','TickDir','out','GridAlpha',.12);
 if kind==4
  shades=[.85 .93 .99;.99 .91 .80];
  for j=1:2
   b=metrics{j}.bounds;gg=rad2deg(b(:,1));
   patch(ax,[gg(1) gg(2) gg(2) gg(1)],[b(1,2) b(1,2) b(2,2) b(2,2)],shades(j,:), ...
   'EdgeColor',[.5 .5 .5],'HandleVisibility','off');
  end
 end
 hs=gobjects(5,1);labels=strings(5,1);
 for k=1:5
  r=runs{k};style='--';tag="not recovered (8 s)";
  if results.classification(k)=="RECOVERED",style='-';tag="recovered";
  elseif results.classification(k)~="NOT_RECOVERED",style=':';tag=results.classification(k);end
  labels(k)=sprintf('FCT %.2f s | %s',FCT(k),tag);
  if kind<=2
   t=r.fault.t;y=r.fault.x(:,kind);
   if isfield(r,'postfault'),t=[t;r.postfault.t];y=[y;r.postfault.x(:,kind)];end
   if kind==1,y=rad2deg(y);yc=results.gamma_clear_deg(k);else,yc=results.nu_clear(k);end
   hs(k)=plot(ax,t,y,style,'Color',colors(k,:),'LineWidth',1.6);
   plot(ax,FCT(k),yc,markers{k},'Color',colors(k,:),'MarkerFaceColor','w','LineWidth',1.5,'MarkerSize',8,'HandleVisibility','off');
  else
   plot(ax,rad2deg(r.fault.x(:,1)),r.fault.x(:,2),':','Color',colors(k,:),'LineWidth',1.1,'HandleVisibility','off');
   if isfield(r,'postfault')
    hs(k)=plot(ax,rad2deg(r.postfault.x(:,1)),r.postfault.x(:,2),style,'Color',colors(k,:),'LineWidth',1.6);
   else
    hs(k)=plot(ax,NaN,NaN,style,'Color',colors(k,:),'LineWidth',1.6);
   end
   plot(ax,results.gamma_clear_deg(k),results.nu_clear(k),markers{k},'Color',colors(k,:),'MarkerFaceColor','w','MarkerSize',9,'LineWidth',1.8,'HandleVisibility','off');
  end
 end
 if kind<=2
  xlabel(ax,'Time from fault inception (s)');xlim(ax,[0 cfg.tEnd]);
  if kind==1,ylabel(ax,'Unwrapped angle, \gamma (deg)');ttl='Rotor angle';fn='01_gamma_vs_time';ref=rad2deg(x0(1));
  else,ylabel(ax,'Normalized frequency, \nu (pu)');ttl='Normalized frequency';fn='02_nu_vs_time';ref=1;end
  yline(ax,ref,':','Color',[.3 .3 .3],'HandleVisibility','off');
 else
  xlabel(ax,'Unwrapped angle, \gamma (deg)');ylabel(ax,'Normalized frequency, \nu (pu)');
  plot(ax,rad2deg(x0(1)),1,'kp','MarkerFaceColor','k','MarkerSize',12,'HandleVisibility','off');
  ttl='Nonlinear phase plane';fn='03_phase_plane';
  if kind==4
   ttl='Clearing states and widened contraction rectangles';fn='04_contraction_zoom';
   xlim(ax,[8 max(results.gamma_clear_deg)+5]);ylim(ax,[.978 1.022]);
   text(ax,11,.981,'Inactive: widened rectangle','FontSize',10);
   text(ax,25,.981,'Active: widened rectangle','FontSize',10);
  end
 end
 title(ax,{[ttl ' | 60% voltage sag'], 'U_{ret} = 0.40; common nominal equilibrium; unchanged recovery test'},'FontWeight','normal');
 grid(ax,'on');legend(ax,hs,labels,'Location','southoutside','NumColumns',2,'Box','off','FontSize',11);
 note='Markers: fault clearing. Solid: recovered; dashed: finite-horizon non-recovery, not proven instability.';
 if kind>=3,note=[note ' Dotted: fault-on; black star: equilibrium.'];end
 if kind==4,note=[note ' Rectangle membership alone is not stability certification.'];end
 annotation(fig,'textbox',[.08 .025 .86 .075],'String',note,'EdgeColor','none','FontSize',11,'HorizontalAlignment','center');
 if kind==1 || kind==3
  pos=[.16 .61 .32 .22];if kind==3,pos=[.48 .46 .40 .20];end
  az=axes(fig,'Position',pos,'Color','w','FontName','Arial','FontSize',10,'Box','on');hold(az,'on');
  for kk=1:5
   rr=runs{kk};ss='--';if results.classification(kk)=="RECOVERED",ss='-';end
   if kind==1
    tt=[rr.fault.t;rr.postfault.t];yy=rad2deg([rr.fault.x(:,1);rr.postfault.x(:,1)]);
    plot(az,tt,yy,ss,'Color',colors(kk,:),'LineWidth',1.2);
    plot(az,FCT(kk),results.gamma_clear_deg(kk),markers{kk},'Color',colors(kk,:),'MarkerFaceColor','w','MarkerSize',5);
   else
    plot(az,rad2deg(rr.fault.x(:,1)),rr.fault.x(:,2),':','Color',colors(kk,:));
    plot(az,rad2deg(rr.postfault.x(:,1)),rr.postfault.x(:,2),ss,'Color',colors(kk,:),'LineWidth',1.2);
    plot(az,results.gamma_clear_deg(kk),results.nu_clear(kk),markers{kk},'Color',colors(kk,:),'MarkerFaceColor','w','MarkerSize',5);
   end
  end
  if kind==1,xlim(az,[0 3]);ylim(az,[10 100]);xlabel(az,'Time (s)');ylabel(az,'\gamma (deg)');
  else,xlim(az,[14 75]);ylim(az,[.995 1.015]);xlabel(az,'\gamma (deg)');ylabel(az,'\nu (pu)');plot(az,rad2deg(x0(1)),1,'kp','MarkerFaceColor','k');end
  title(az,'Zoom: clearing and recovery','FontWeight','normal','FontSize',10);grid(az,'on');
 end
 export_figure(fig,out,fn);close(fig);
end

end
function export_figure(fig,out,name)
[~,id]=fileparts(tempname);tmp=fullfile(getenv('SystemRoot'),'Temp',['fct_' id]);mkdir(tmp);
cleanup=onCleanup(@()rmdir(tmp,'s'));
exportgraphics(fig,fullfile(tmp,'plot.png'),'Resolution',300,'BackgroundColor','white');
exportgraphics(fig,fullfile(tmp,'plot.pdf'),'ContentType','vector','BackgroundColor','white');
savefig(fig,fullfile(tmp,'plot.fig'));
for ext={'.png','.pdf','.fig'},copyfile(fullfile(tmp,['plot' ext{1}]),fullfile(out,[name ext{1}]));end
end




