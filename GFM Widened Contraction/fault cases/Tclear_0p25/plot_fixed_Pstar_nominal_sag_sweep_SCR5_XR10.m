function plot_fixed_Pstar_nominal_sag_sweep_SCR5_XR10(caseStart)
if nargin<1,caseStart=0;end
%% Plot saved nominal-equilibrium sag cases; no simulation is rerun.
clc;root=fileparts(mfilename('fullpath'));out=fullfile(root,'grafikler');caseOut=fullfile(out,'case_by_case');if ~isfolder(out),mkdir(out);end;if ~isfolder(caseOut),mkdir(caseOut);end
S=readtable(fullfile(root,'fixed_Pstar_sag_sweep.csv'),'TextType','string');N=height(S);assert(N==12&&all(S.Pstar==1)&&all(S.nu0==1)&&all(S.Tclear_s==S.Tclear_s(1)));
T=cell(N,1);events=cell(N,1);colors=turbo(N);tc=S.Tclear_s(1);te=8;eq=S.gamma0_deg(1);Imax=1.2;Vmax=1.3;
for k=1:N
 T{k}=readtable(fullfile(root,sprintf('trajectory_Uret_%.2f.csv',S.Uret(k))),'TextType','string');events{k}=readtable(fullfile(root,sprintf('limiter_events_Uret_%.2f.csv',S.Uret(k))),'TextType','string');d=T{k};
 assert(all(diff(d.time_s)>=0)&&abs(d.time_s(1))<1e-12&&abs(d.gamma_deg(1)-eq)<1e-9&&abs(d.nu(1)-1)<1e-12);
 assert(max(d.limited_current_pu)<=Imax+1e-9&&max(abs(d.limited_current_pu-min(d.raw_current_pu,Imax)))<1e-8);
 ix=find(d.stage=="fault",1,'last');assert(abs(d.time_s(ix)-tc)<1e-9&&abs(d.gamma_deg(ix)-S.gamma_clear_deg(k))<1e-7);
 assert(abs(d.gamma_deg(end)-S.final_gamma_deg(k))<1e-7);
end
rec=S.classification=="RECOVERED";non=S.classification=="NON_RECOVERED_WITHIN_HORIZON";assert(all(rec|non|S.classification=="NOT_ASSESSED_VALIDITY_STOP"));
set(groot,'defaultAxesFontName','Arial','defaultAxesFontSize',11,'defaultLineLineWidth',1.35);
if caseStart==0
% 1: fault onset and clearing, with all cases on the same time axis.
f=figure('Visible','off','Color','w','Position',[30 30 1450 1000]);tiledlayout(2,2,'TileSpacing','compact','Padding','compact');fields={'gamma_deg','nu','raw_current_pu','limited_current_pu'};yl={'gamma (deg)','nu (pu)','Raw current demand (pu)','Limited current (pu)'};
for j=1:4,ax=nexttile;hold on;h=gobjects(N,1);for k=1:N,d=T{k};q=d.time_s<=.8;h(k)=plot(d.time_s(q),d.(fields{j})(q),'Color',colors(k,:));end;ylabel(yl{j});if j>2,yline(Imax,'k:','1.2 pu','HandleVisibility','off');end;decorate(ax,tc,.8);if j==1,legend(h,compose('Uret=%.2f',S.Uret),'Location','eastoutside');end;end
sgtitle(sprintf('Same nominal initial state | fault 0-%.3g s | all 12 retained-voltage cases',tc));saveplot(f,out,'01_fault_and_clearing_zoom');
% 2 and 3: separate the scales of recovered and non-recovered cases.
groupplot(find(rec),T,S,colors,tc,te,eq,out,'02_recovered_cases','Recovered within 8 s');
groupplot(find(non),T,S,colors,tc,te,eq,out,'03_nonrecovered_cases','Not recovered within 8 s (unwrapped angle)');
% 4: limiter activity, exact event times from saved solver event records.
f=figure('Visible','off','Color','w','Position',[30 30 1450 700]);hold on;
for k=1:N
 ev=events{k};active=T{k}.raw_current_pu(1)>Imax+1e-9;start=0;
 for j=1:height(ev)
  t=ev.time_s(j);if ev.type(j)=="activation",if ~active,start=t;active=true;end;elseif ev.type(j)=="release",if active,plot([start t],[k k],'-','Color',colors(k,:),'LineWidth',7,'HandleVisibility','off');active=false;end;end
 end
 if active,plot([start T{k}.time_s(end)],[k k],'-','Color',colors(k,:),'LineWidth',7,'HandleVisibility','off');end
end
xline(tc,'k--','Fault cleared','LabelVerticalAlignment','bottom');xlim([0 te]);ylim([.5 N+.5]);yticks(1:N);yticklabels(compose('Uret %.2f | %s',S.Uret,repmat("not recovered",N,1)));labs=compose('Uret %.2f | not recovered',S.Uret);labs(rec)=compose('Uret %.2f | recovered',S.Uret(rec));yticklabels(labs);set(gca,'YDir','reverse');xlabel('Time (s)');title({'Limiter-active intervals from saved activation/release events','Colored segment: limiter active; blank: inactive'});grid on;saveplot(f,out,'04_limiter_activity');
% 5: fault paths and clearing points; no recovery guarantee inferred.
f=figure('Visible','off','Color','w','Position',[30 30 1250 800]);hold on;draw_regions();h=gobjects(N,1);
for k=1:N,d=T{k};q=d.stage=="fault";h(k)=plot(d.gamma_deg(q),d.nu(q),'Color',colors(k,:));marker='x';if rec(k),marker='o';end;plot(S.gamma_clear_deg(k),S.nu_clear(k),marker,'Color',colors(k,:),'MarkerSize',9,'LineWidth',1.8,'HandleVisibility','off');text(S.gamma_clear_deg(k)+.3,S.nu_clear(k),sprintf('%.2f',S.Uret(k)),'FontSize',9,'Color',colors(k,:));end
plot(eq,1,'kp','MarkerFaceColor','k','MarkerSize',12,'HandleVisibility','off');xlabel('gamma (deg)');ylabel('nu (pu)');grid on;title({sprintf('Fault-on paths and widened postfault rectangles | Tclear=%.3g s',tc),'Clearing markers: circle = recovered, cross = not recovered within 8 s'});legend(h,compose('Uret=%.2f',S.Uret),'Location','eastoutside');saveplot(f,out,'05_fault_paths_and_clearing_states');
% 6: transition between the last recovered and first non-recovered test.
idx=find(abs(S.Uret-.60)<1e-9|abs(S.Uret-.55)<1e-9);f=figure('Visible','off','Color','w','Position',[30 30 1450 1000]);tiledlayout(2,2,'TileSpacing','compact');fields={'gamma_deg','nu','Pe','limited_current_pu'};yl={'Unwrapped gamma (deg)','nu (pu)','Electrical power Pe (pu)','Limited current (pu)'};
for j=1:4,ax=nexttile;hold on;h=gobjects(2,1);for q=1:2,k=idx(q);d=T{k};mask=d.time_s<=2;h(q)=plot(d.time_s(mask),d.(fields{j})(mask),'Color',colors(k,:));end;ylabel(yl{j});if j==3,yline(1,'k:','Pstar=1','HandleVisibility','off');elseif j==4,yline(Imax,'k:','Imax','HandleVisibility','off');end;decorate(ax,tc,2);if j==1,legend(h,compose('Uret=%.2f: %s',S.Uret(idx),S.classification(idx)),'Location','northwest');end;end
sgtitle(sprintf('Saved Uret=0.60 and 0.55 cases | Tclear=%.3g s',tc));saveplot(f,out,'06_recovery_transition_060_vs_055');
end
% Individual six-panel reports for each saved case.
for k=max(1,caseStart):N
 d=T{k};f=figure('Visible','off','Color','w','Position',[30 30 1450 1200]);tiledlayout(3,2,'TileSpacing','compact','Padding','compact');
 ax=nexttile;plot(d.time_s,d.gamma_deg,'Color',colors(k,:));ylabel('Unwrapped gamma (deg)');yline(eq,'k:','Nominal','HandleVisibility','off');decorate(ax,tc,te);
 ax=nexttile;plot(d.time_s,d.nu,'Color',colors(k,:));ylabel('nu (pu)');yline(1,'k:','Nominal','HandleVisibility','off');decorate(ax,tc,te);
 ax=nexttile;hold on;h1=plot(d.time_s,d.raw_current_pu,'Color',[.7 .2 .2]);h2=plot(d.time_s,d.limited_current_pu,'Color',[.1 .35 .8]);yline(Imax,'k:','Imax=1.2','HandleVisibility','off');ylabel('Current (pu)');legend([h1 h2],{'Raw demand','Limited current'},'Location','best');decorate(ax,tc,te);
 ax=nexttile;plot(d.time_s,d.Pe,'Color',colors(k,:));ylabel('Electrical power Pe (pu)');yline(1,'k:','Pstar=1','HandleVisibility','off');decorate(ax,tc,te);
 ax=nexttile;plot(d.time_s,d.required_voltage_pu,'Color',colors(k,:));ylabel('Required voltage command (pu)');yline(Vmax,'r:','Vmax=1.3','HandleVisibility','off');decorate(ax,tc,te);
 nexttile;hold on;draw_regions();q=d.stage=="fault";h1=plot(d.gamma_deg(q),d.nu(q),'Color',[.85 .3 .1]);h2=plot(d.gamma_deg(~q),d.nu(~q),'Color',[.1 .35 .8]);h3=plot(S.gamma_clear_deg(k),S.nu_clear(k),'kx','MarkerSize',10,'LineWidth',2);h4=plot(eq,1,'kp','MarkerFaceColor','k','MarkerSize',11);xlabel('Unwrapped gamma (deg)');ylabel('nu (pu)');grid on;legend([h1 h2 h3 h4],{'Fault-on','Postfault','Clearing','Nominal'},'Location','best');
 status='NOT RECOVERED WITHIN 8 s';if rec(k),status='RECOVERED';end
 sgtitle(sprintf('Uret=%.2f | sag %.0f%% | Tclear=%.3g s | %s',S.Uret(k),S.sag_depth_percent(k),tc,status));saveplot(f,caseOut,sprintf('case_Uret_%03d',round(100*S.Uret(k))));
end
writetable(S(:,{'Uret','sag_depth_percent','classification','gamma_clear_deg','nu_clear','limiter_encountered','limiter_release_time_s','max_gamma_deg','max_raw_current_pu','max_limited_current_pu','Tclear_s','inside_Cpost','inside_inactive_Cpost','inside_any_Cpost'}),fullfile(out,'plot_case_index.csv'));
fprintf('Completed: 6 overview figures and 12 individual case figures; PNG and FIG.\n');
end
function decorate(ax,tc,te)
axes(ax);hold on;xlim([0 te]);yl=ylim;hh=patch([0 tc tc 0],[yl(1) yl(1) yl(2) yl(2)],[.9 .7 .3],'FaceAlpha',.12,'EdgeColor','none','HandleVisibility','off');uistack(hh,'bottom');ylim(yl);xline(tc,'k--','Clear','HandleVisibility','off','LabelVerticalAlignment','bottom');xlabel('Time (s)');grid on;
end
function saveplot(f,out,name)
[~,id]=fileparts(tempname);tmp=fullfile(getenv('SystemRoot'),'Temp',[id '.png']);drawnow;print(f,tmp,'-dpng','-r160');copyfile(tmp,fullfile(out,[name '.png']));delete(tmp);savefig(f,fullfile(out,[name '.fig']));close(f);
end
function groupplot(idx,T,S,colors,tc,te,eq,out,name,label)
if isempty(idx),f=figure('Visible','off','Color','w');axis off;text(.1,.5,sprintf('No cases in this category at Tclear=%.3g s',tc));title(label);saveplot(f,out,name);return;end
f=figure('Visible','off','Color','w','Position',[30 30 1450 1050]);tiledlayout(2,2,'TileSpacing','compact','Padding','compact');fields={'gamma_deg','nu','limited_current_pu'};yl={'Unwrapped gamma (deg)','nu (pu)','Limited current (pu)'};
for j=1:3,ax=nexttile;hold on;h=gobjects(numel(idx),1);for q=1:numel(idx),k=idx(q);d=T{k};h(q)=plot(d.time_s,d.(fields{j}),'Color',colors(k,:));end;ylabel(yl{j});if j==1,yline(eq,'k:','Nominal','HandleVisibility','off');elseif j==2,yline(1,'k:','Nominal','HandleVisibility','off');else,yline(1.2,'k:','Imax','HandleVisibility','off');end;decorate(ax,tc,te);if j==1,legend(h,compose('Uret=%.2f',S.Uret(idx)),'Location','eastoutside');end;end
nexttile;hold on;draw_regions();for k=idx',d=T{k};plot(d.gamma_deg,d.nu,'Color',colors(k,:));plot(S.gamma_clear_deg(k),S.nu_clear(k),'x','Color',colors(k,:),'MarkerSize',8,'HandleVisibility','off');end;plot(eq,1,'kp','MarkerFaceColor','k','MarkerSize',11);xlabel('Unwrapped gamma (deg)');ylabel('nu (pu)');grid on;title('Phase plane; crosses mark clearing');sgtitle(sprintf('%s | same nominal state | Tclear=%.3g s',label,tc));saveplot(f,out,name);
end

function draw_regions()
[A,I]=load_widened_regions();regions={I,A};colors=[.12 .4 .85;.9 .42 .1];for j=1:2,b=regions{j};patch(b(1,[1 2 2 1]),b(2,[1 1 2 2]),colors(j,:),'FaceAlpha',.1,'EdgeColor',colors(j,:),'LineWidth',1.2,'HandleVisibility','off');end
end
