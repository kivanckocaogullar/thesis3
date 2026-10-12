function plot_widened_containment(out)
% Reuses saved samples and the established mode-aware containment procedure.
a=load(fullfile(out,'widened_fct_analysis.mat'));mi=a.metrics{1};ma=a.metrics{2};
colors=[0 .45 .70;0 .62 .45;.80 .47 .05;.65 .25 .65;.85 .25 .10];
markers={'o','s','d','^','v'};
summary=table();contained=cell(5,1);firstOut=nan(5,1);reentry=nan(5,1);
for k=1:5
 s=a.runs{k}.postfault;t=s.t(:);x=s.x;tc=a.FCT(k);
 assert(all(diff(t)>0) && abs(t(1)-tc)<1e-12);
 d=coupled(x(:,1),1,a.p,a.n,0); % Existing restored-nominal electrical closure.
 inI=~d.active & all(x>=mi.bounds(1,:) & x<=mi.bounds(2,:),2);
 inA=d.active & all(x>=ma.bounds(1,:) & x<=ma.bounds(2,:),2);
 inside=inI|inA;contained{k}=inside;
 idx=find(~inside,1);re="No";
 if isempty(idx),first=NaN;re="Not applicable";
 else
  first=t(idx)-tc;firstOut(k)=idx;
  next=find(inside & (1:numel(t))'>idx,1);
  if ~isempty(next)
   reentry(k)=next;
   if inside(1),re=sprintf('Yes (%.6f s)',t(next)-tc);
   else,re=sprintf('Initial entry (%.6f s)',t(next)-tc);end
  end
 end
 % Time-weighted estimate for nonuniform saved samples; no new trajectory samples.
 pct=100*trapz(t,double(inside))/(t(end)-t(1));
 summary=[summary;table(tc,a.results.classification(k),inside(1),all(inside),first,string(re),pct, ...
 'VariableNames',{'FCT_s','Recovery','Clearing_Inside','All_Samples_Inside','First_Exit_s_after_clearing','Re_entry','Time_Inside_percent'})];
 fprintf('FCT %.2f: %d/%d samples inside, first outside %.9f s, time inside %.9f%%\n',tc,nnz(inside),numel(t),first,pct);
end
format long g;disp(summary);
writetable(summary,fullfile(out,'containment_summary.csv'));
fig=figure('Visible','off','Color','w','Position',[50 50 1450 1000],'Renderer','painters');
ax=axes(fig,'Position',[.085 .35 .87 .56]);hold(ax,'on');
set(ax,'FontName','Arial','FontSize',13,'LineWidth',1,'TickDir','out','Box','off','GridAlpha',.12);
col=[.08 .40 .68;.83 .35 .08];regions={mi,ma};
h=gobjects(0);labs={};
for j=1:2
 b=regions{j}.bounds;g=rad2deg(b(:,1));v=b(:,2);
 h(end+1)=plot(ax,g([1 2 2 1 1]),v([1 1 2 2 1]),'-','Color',col(j,:),'LineWidth',2.2);
 modes={'inactive','active'};labs{end+1}=['Certified ' modes{j} ' rectangle'];
end
for k=1:5
 s=a.runs{k}.postfault;
 plot_segments(ax,rad2deg(s.x(:,1)),s.x(:,2),contained{k},colors(k,:));
 h(end+1)=plot(ax,rad2deg(s.x(1,1)),s.x(1,2),markers{k},'Color',colors(k,:),'MarkerFaceColor','w','LineWidth',1.7,'MarkerSize',9);
 labs{end+1}=sprintf('FCT %.2f s | clearing',a.FCT(k));
end
h(end+1)=plot(ax,rad2deg(a.x0(1)),1,'kp','MarkerFaceColor','k','MarkerSize',12);labs{end+1}='Nominal equilibrium';
h(end+1)=plot(ax,NaN,NaN,'k-','LineWidth',1.6);labs{end+1}='Inside correct-mode rectangle';
h(end+1)=plot(ax,NaN,NaN,'k--','LineWidth',1.6);labs{end+1}='Outside / uncertified gap';
xlim(ax,[8 75]);ylim(ax,[.975 1.024]);grid(ax,'on');
xlabel(ax,'Unwrapped angle, \gamma (deg)');ylabel(ax,'Normalized frequency, \nu (pu)');
title(ax,{'Saved post-fault trajectories | 60% voltage sag', ...
 'Widened interval-certified rectangles (M < 10I); sampled containment only'},'FontWeight','normal');
lg=legend(ax,h,labs,'NumColumns',3,'FontSize',10,'Box','off');
lg.Units='normalized';lg.Position=[.20 .135 .65 .10];
ax.Position=[.085 .32 .87 .58];
% Full unwrapped trajectories in one inset within the same figure.
az=axes(fig,'Position',[.64 .355 .29 .135]);hold(az,'on');
set(az,'FontName','Arial','FontSize',9,'Box','on','Color','w');
for k=1:5
 s=a.runs{k}.postfault;
 plot_segments(az,rad2deg(s.x(:,1)),s.x(:,2),contained{k},colors(k,:));
 plot(az,rad2deg(s.x(1,1)),s.x(1,2),markers{k},'Color',colors(k,:),'MarkerSize',4);
end
plot(az,rad2deg(a.x0(1)),1,'kp','MarkerFaceColor','k','MarkerSize',5);
xlim(az,[0 5300]);ylim(az,[.99 1.056]);grid(az,'on');
title(az,'Full saved trajectories (unwrapped)','FontWeight','normal','FontSize',10);
xlabel(az,'\gamma (deg)');ylabel(az,'\nu (pu)');
annotation(fig,'textbox',[.085 .02 .87 .10],'String', ...
 {'Solid links join two inside samples; links touching an outside sample are dashed. The uncertified switching gap is excluded.', ...
 'Time inside: trapezoidal time-weighted estimate from saved sample membership. Main view is zoomed; inset retains complete excursions.', ...
 'Sampled numerical containment is not a continuous-time containment proof or full stability certification.'}, ...
 'HorizontalAlignment','center','EdgeColor','none','FontSize',10.5);
[~,id]=fileparts(tempname);temp=fullfile(getenv('SystemRoot'),'Temp',['containment_' id '.png']);
exportgraphics(fig,temp,'Resolution',300,'BackgroundColor','white');
copyfile(temp,fullfile(out,'trajectory_containment.png'));delete(temp);
savefig(fig,fullfile(out,'trajectory_containment.fig'));
exportgraphics(fig,fullfile(out,'trajectory_containment.pdf'),'ContentType','vector');
save(fullfile(out,'containment_data.mat'),'summary','contained','mi','ma','firstOut','reentry');
close(fig);
fprintf('First exit = 0 if already outside at clearing; NaN means no sampled exit.\n');
fprintf('Saved only summary CSV and PNG in: %s\n',out);

end
function plot_segments(ax,g,v,inside,color)
% Preserve every saved sample; use conservative dashed transition links.
good=inside(1:end-1)&inside(2:end);
for type=1:2
 keep=good;if type==2,keep=~good;end
 edges=diff([false;keep(:);false]);starts=find(edges==1);stops=find(edges==-1);
 X=[];Y=[];
 for r=1:numel(starts)
  ix=starts(r):stops(r);X=[X;g(ix);NaN];Y=[Y;v(ix);NaN];
 end
 style='-';if type==2,style='--';end
 plot(ax,X(:),Y(:),style,'Color',color,'LineWidth',1.5,'HandleVisibility','off');
end
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




