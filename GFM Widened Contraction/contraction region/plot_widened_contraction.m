%% Current-limited nominal postfault contraction: WIDENED interval-verified rectangles.
% Run this file to check both widened metrics and produce region and J/M/Sym(F) plots.
% Adapted from GFM Contraction Analsys/limiter region/plot_current_limited_contraction.m.
% No new fault simulation, no energy sublevel, no assertion of common metric.
clearvars;clc;root=fileparts(mfilename('fullpath'));out=fullfile(root,'plots');if ~isfolder(out),mkdir(out);end
src={fullfile('metrics','inactive'),fullfile('metrics','active')};old={fullfile('reference_old_metrics','inactive'),fullfile('reference_old_metrics','active')};
a=load(fullfile(root,'inputs','model_inputs.mat'),'p','n');p=a.p;n=a.n;p.Pstar=1;p.omega=100*pi;
assert(abs(p.mass-6.8)<1e-12&&p.D==25&&p.Imax==1.2&&p.Vmax==1.3);
metrics=cell(2,1);certificates=cell(2,1);labels={'Limiter inactive','Limiter active'};checks=table();rng(42);
for j=1:2
 metrics{j}=load(fullfile(root,src{j},'metric.mat'));m=metrics{j};certificate=jsondecode(fileread(fullfile(root,src{j},'interval_verification.json')));assert(strcmp(certificate.status,'VERIFIED')&&certificate.unresolved_boxes==0&&certificate.counterexample_boxes==0);certificates{j}=certificate;
 jsonMetric=jsondecode(fileread(fullfile(root,src{j},'metric.json')));assert(max(abs(m.C-permute(jsonMetric.coefficients,[2 3 1])),[],'all')<1e-12,'MAT/JSON metric mismatch');
 [g,v]=ndgrid(linspace(-1,1,301));z=verify(m.C,[g(:),v(:)],m.lambda,m.center,m.scale,p,n);
 assert(min(z.minM)>m.epsilon&&max(z.maxM)<m.upper&&max(z.maxQ)<-m.epsilon);
 if j==1,assert(max(z.rawI)<p.Imax);else,assert(min(z.rawI)>p.Imax);end
 assert(max(z.u)<p.Vmax);
 maxMdError=0;maxJError=0;maxSymF=-Inf;
 for k=1:30
  zz=-.9+1.8*rand(1,2);x=m.center+m.scale.*zz;d=coupled(x(1),1,p,n,0);nu=x(2);f=[p.omega*(nu-1),((1-d.P)/nu-p.D*(nu-1))/p.mass];J=jac(x,p,d);[M,MD]=metric(m.C,zz,f./m.scale);h=1e-7;MP=metric(m.C,zz+h*f./m.scale,[0 0]);MM=metric(m.C,zz-h*f./m.scale,[0 0]);maxMdError=max(maxMdError,norm((MP-MM)/(2*h)-MD,'fro')/(1+norm(MD,'fro')));
  Jn=zeros(2);for c=1:2,dx=zeros(1,2);dx(c)=1e-6;Jn(:,c)=(field(x+dx,p,n)-field(x-dx,p,n))'/2e-6;end;maxJError=max(maxJError,norm(J-Jn,'fro'));R=chol(M);S=(R'\(J'*M+M*J+MD)/R)/2;maxSymF=max(maxSymF,max(eig(S)));
 end
 assert(maxMdError<1e-6&&maxJError<1e-6&&maxSymF<-m.lambda);
 checks=[checks;table(string(labels{j}),m.lambda,rad2deg(m.bounds(1,1)),rad2deg(m.bounds(2,1)),m.bounds(1,2),m.bounds(2,2),min(z.minM),max(z.maxM),max(z.maxQ),max(z.u),maxMdError,maxJError,maxSymF,'VariableNames',{'branch','lambda','gamma_lo_deg','gamma_hi_deg','nu_lo','nu_hi','dense_min_M','dense_max_M','dense_max_Q','dense_max_voltage','Mdot_relative_error','Jacobian_error','sample_max_SymF'})];
end
writetable(checks,fullfile(out,'metric_checks.csv'));disp(checks);
ge=fzero(@(g)coupled(g,1,p,n,0).P-1,[.3 .4]);gs=fzero(@(g)coupled(g,1,p,n,1).rawI-p.Imax,[.42 .43]);
colors=[.12 .4 .85;.9 .42 .1];fig=figure('Visible','off','Color','w','Position',[50 50 1350 800]);hold on;handles=gobjects(4,1);
for j=1:2,b=metrics{j}.bounds;gd=rad2deg(b(:,1));handles(j)=patch([gd(1) gd(2) gd(2) gd(1)],[b(1,2) b(1,2) b(2,2) b(2,2)],colors(j,:),'FaceAlpha',.2,'EdgeColor',colors(j,:),'LineWidth',2);end
handles(3)=xline(rad2deg(gs),'k--','LineWidth',1.4);handles(4)=plot(rad2deg(ge),1,'kp','MarkerFaceColor','k','MarkerSize',12);
xlabel('gamma (deg)');ylabel('nu (pu)');grid on;box on;xlim([8 45]);ylim([.975 1.025]);set(gca,'Layer','top');
for j=1:2,b=metrics{j}.bounds;text(mean(rad2deg(b(:,1))),b(2,2)-.0015,{sprintf('gamma [%.2f, %.2f] deg',rad2deg(b(:,1))),sprintf('nu [%.2f, %.2f]',b(:,2)),sprintf('lambda = %g s^{-1}',metrics{j}.lambda)},'HorizontalAlignment','center','VerticalAlignment','top','FontSize',10,'Color',colors(j,:)*.8);end
title({'Current-limited nominal postfault model | Pstar=1 | widened rectangles','Interval-verified on the continuous rectangles: 0 < M < 10I, J^TM+MJ+Mdot+2\lambdaM < 0, limiter mode and Vmax'});
legend(handles,{'Interval-verified inactive contraction rectangle','Interval-verified active contraction rectangle','Limiter switching boundary','Nominal equilibrium'},'Location','southoutside','NumColumns',2);
text(.015,.03,{'Rectangles are not claimed invariant; no common switching or region-of-attraction certificate is inferred.','Active branch: strict inequalities certified; the 1e-5 optimization margin is not certified in every box.'},'Units','normalized','VerticalAlignment','bottom','FontSize',9,'BackgroundColor','w','Margin',3);
exportpng(fig,fullfile(out,'current_limited_contraction_regions.png'),'Resolution',180);savefig(fig,fullfile(out,'current_limited_contraction_regions.fig'));close(fig);
% Old (16.4-33.9 deg) versus widened rectangles.
fig=figure('Visible','off','Color','w','Position',[50 50 1350 800]);hold on;handles=gobjects(6,1);
for j=1:2
 b=metrics{j}.bounds;gd=rad2deg(b(:,1));handles(j)=patch([gd(1) gd(2) gd(2) gd(1)],[b(1,2) b(1,2) b(2,2) b(2,2)],colors(j,:),'FaceAlpha',.15,'EdgeColor',colors(j,:),'LineWidth',2);
 o=load(fullfile(root,old{j},'metric.mat'),'bounds');b=o.bounds;gd=rad2deg(b(:,1));handles(j+2)=patch([gd(1) gd(2) gd(2) gd(1)],[b(1,2) b(1,2) b(2,2) b(2,2)],colors(j,:),'FaceAlpha',.35,'EdgeColor','k','LineStyle','--','LineWidth',1.4);
end
handles(5)=xline(rad2deg(gs),'k--','LineWidth',1.4);handles(6)=plot(rad2deg(ge),1,'kp','MarkerFaceColor','k','MarkerSize',12);
xlabel('gamma (deg)');ylabel('nu (pu)');grid on;box on;xlim([8 45]);ylim([.975 1.025]);set(gca,'Layer','top');
title({'Old versus widened interval-verified contraction rectangles','Same model, same lambda values and M < 10I; separate metrics per limiter mode'});
legend(handles,{'Widened inactive (interval-verified)','Widened active (interval-verified)','Old inactive (interval-verified)','Old active (interval-verified)','Limiter switching boundary','Nominal equilibrium'},'Location','southoutside','NumColumns',3);
exportpng(fig,fullfile(out,'old_vs_widened_rectangles.png'),'Resolution',180);savefig(fig,fullfile(out,'old_vs_widened_rectangles.fig'));close(fig);
% Inside/outside samples at nu=1; use the same saved metric without refitting.
fig=figure('Visible','off','Color','w','Position',[50 50 1500 980]);tiledlayout(2,3,'TileSpacing','compact','Padding','compact');spectra=table();pointColors=lines(3);
for j=1:2
 m=metrics{j};if j==1,angles=[0 rad2deg(ge) 5];pointLabels={'far outside below','inside','outside below'};else,angles=[50 mean(rad2deg(m.bounds(:,1))) 60];pointLabels={'outside above','inside','far outside above'};end;
 EV=cell(1,3);for c=1:3,EV{c}=zeros(2,3);end
 for k=1:3
  x=[deg2rad(angles(k)),1];zz=(x-m.center)./m.scale;d=coupled(x(1),1,p,n,0);J=jac(x,p,d);f=[0,(1-d.P)/p.mass];[M,MD]=metric(m.C,zz,f./m.scale);assert(all(d.active==(j==2)),'Sample belongs to the wrong limiter mode');EV{1}(:,k)=eig(J);EV{2}(:,k)=sort(eig(M));[R,cholFlag]=chol(M);EV{3}(:,k)=NaN;if cholFlag==0,S=(R'\(J'*M+M*J+MD)/R)/2;EV{3}(:,k)=sort(eig(S));end;if k==2,assert(min(EV{2}(:,k))>0&&max(EV{3}(:,k))<-m.lambda);end
  spectra=[spectra;table(string(labels{j}),angles(k),1,real(EV{1}(1,k)),imag(EV{1}(1,k)),real(EV{1}(2,k)),imag(EV{1}(2,k)),EV{2}(1,k),EV{2}(2,k),EV{3}(1,k),EV{3}(2,k),'VariableNames',{'branch','gamma_deg','nu','J1_real','J1_imag','J2_real','J2_imag','M1','M2','SymF1','SymF2'})];
 end
 names={'Eigenvalues of J','Eigenvalues of M','Eigenvalues of Sym(F), with Mdot'};
 for c=1:3
  nexttile;hold on;for k=1:3,plot(real(EV{c}(:,k)),imag(EV{c}(:,k)),'o','Color',pointColors(k,:),'LineWidth',1.5,'MarkerSize',7);end
  if c==2,if all(EV{2}>0,'all'),set(gca,'XScale','log');xlabel('Real part (log scale)');else,xline(0,'k:','HandleVisibility','off');xlabel('Real part');end;ylim([-.1 .1]);else,xline(0,'k:','HandleVisibility','off');xlabel('Real part (1/s)');end
  if c==3,xline(-m.lambda,'k--','HandleVisibility','off');ylim([-.1 .1]);end
  ylabel('Imaginary part');grid on;title({labels{j},names{c}});
  if c==1,legend(compose('gamma=%.3f deg | %s',angles(:),string(pointLabels(:))),'Location','southoutside');end
 end
end
sgtitle({'Widened metrics | Pstar=1, nu=1 | inside/outside interval-verified rectangles','Outside does not imply non-contraction; Sym(F) omitted if M is not positive definite'});spectra.inside_rectangle=repmat([false;true;false],2,1);spectra.metric_positive=spectra.M1>0;spectra.pointwise_contracting=spectra.metric_positive & spectra.SymF2<0;spectra.required_rate_met=spectra.metric_positive & spectra.SymF2<repelem(-checks.lambda,3);disp(spectra);
exportpng(fig,fullfile(out,'current_limited_J_M_SymF.png'),'Resolution',180);savefig(fig,fullfile(out,'current_limited_J_M_SymF.fig'));close(fig);writetable(spectra,fullfile(out,'matrix_spectra.csv'));
% Continuous eigenvalue slices, with a deliberate gap at unvalidated switching band.
fig=figure('Visible','off','Color','w','Position',[50 50 1300 920]);tiledlayout(3,1,'TileSpacing','compact');allSlices=cell(2,1);
for j=1:2
 m=metrics{j};g=linspace(m.bounds(1,1),m.bounds(2,1),401);data=zeros(numel(g),6);
 for k=1:numel(g),x=[g(k),1];d=coupled(g(k),1,p,n,0);J=jac(x,p,d);[M,MD]=metric(m.C,(x-m.center)./m.scale,[0,(1-d.P)/p.mass]./m.scale);R=chol(M);S=(R'\(J'*M+M*J+MD)/R)/2;data(k,:)=[sort(real(eig(J)))',sort(eig(M))',sort(eig(S))'];end
 allSlices{j}=struct('gamma',rad2deg(g),'data',data);
end
for c=1:3
 nexttile;hold on;h=gobjects(2,1);for j=1:2,s=allSlices{j};h(j)=plot(s.gamma,s.data(:,2*c-1),'Color',colors(j,:),'LineWidth',1.5);plot(s.gamma,s.data(:,2*c),'--','Color',colors(j,:),'LineWidth',1.5,'HandleVisibility','off');end;xline(rad2deg(gs),'k:','HandleVisibility','off');grid on;xlabel('gamma (deg)');
 if c==1,ylabel('Re eig(J) (1/s)');elseif c==2,set(gca,'YScale','log');ylabel('eig(M), log scale');else,ylabel('eig(Sym(F)) (1/s)');yline(0,'k:','HandleVisibility','off');end
 legend(h,labels,'Location','eastoutside');
end
sgtitle('Current-limited model: eigenvalue slices at nu=1 (widened interval-verified branches)');exportpng(fig,fullfile(out,'current_limited_matrix_slices.png'),'Resolution',180);savefig(fig,fullfile(out,'current_limited_matrix_slices.fig'));close(fig);
save(fullfile(out,'plot_data.mat'),'metrics','certificates','checks','spectra','allSlices','ge','gs');
cf=@(f)cellfun(@(c)c.(f),certificates);
cert=table(string(labels(:)),string(cellfun(@(c)c.status,certificates,'UniformOutput',false)),cf('certified_min_eig_M_lower_bound'),cf('certified_min_eig_10I_minus_M_lower_bound'),cf('certified_max_eig_contraction_matrix_upper_bound'),cf('epsilon_margin_holds_in_all_boxes'),cf('verified_boxes'),cf('unresolved_boxes'),cf('counterexample_boxes'),cf('min_certified_current_margin_pu'),cf('min_certified_voltage_margin_pu'),'VariableNames',{'branch','status','min_eig_M_lb','min_eig_10I_minus_M_lb','max_eig_Q_ub','eps_margin_all_boxes','verified_boxes','unresolved_boxes','counterexample_boxes','current_margin_pu','voltage_margin_pu'});
writetable(cert,fullfile(out,'interval_certificates.csv'));disp(cert);
function J=jac(x,p,d)
v=x(2);J=[0 p.omega;-d.totalPg/(p.mass*v),(-(1-d.P)/v^2-p.D)/p.mass];
end
function f=field(x,p,n)
d=coupled(x(1),1,p,n,0);f=[p.omega*(x(2)-1),((1-d.P)/x(2)-p.D*(x(2)-1))/p.mass];
end
function v=verify(C,grid,lam,cen,sc,p,n)
 x=cen+grid.*sc;d=coupled(x(:,1),1,p,n,0);v=struct('minM',zeros(size(grid,1),1),'maxM',zeros(size(grid,1),1),'maxQ',zeros(size(grid,1),1),'eigJ',zeros(size(grid,1),1),'rawI',d.rawI,'u',abs(d.u));
 for k=1:size(grid,1)
  nu=x(k,2);f=[p.omega*(nu-1),((1-d.P(k))/nu-p.D*(nu-1))/p.mass];J=[0 p.omega;-d.totalPg(k)/(p.mass*nu),(-(1-d.P(k))/nu^2-p.D)/p.mass];
  [M,Mdot]=metric(C,grid(k,:),f./sc);Q=J'*M+M*J+Mdot+2*lam*M;
  v.minM(k)=min(eig(M));v.maxM(k)=max(eig(M));v.maxQ(k)=max(eig(Q));v.eigJ(k)=max(real(eig(J)));
 end
end
function [M,Mdot]=metric(C,z,dz)
 b=[1 z z(1)^2 z(1)*z(2) z(2)^2];bd=[0 dz 2*z(1)*dz(1) z(2)*dz(1)+z(1)*dz(2) 2*z(2)*dz(2)];M=zeros(2);Mdot=M;for j=1:6,M=M+b(j)*C(:,:,j);Mdot=Mdot+bd(j)*C(:,:,j);end
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








function exportpng(fig,destination,varargin)
% Export through ASCII path to avoid PNG library issues with Unicode folders.
[~,id]=fileparts(tempname);tmp=fullfile('C:\Windows\Temp',[id '.png']);
cleanup=onCleanup(@()delete_if_exists(tmp));
exportgraphics(fig,tmp,varargin{:});copyfile(tmp,destination,'f');
end
function delete_if_exists(f)
if isfile(f),delete(f);end
end


