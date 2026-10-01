%% SCR5/XR10: hypothetical current ratings 1.2, 1.3 and 1.4 pu.
% Same network, H=3.4, K=20, Kd=5, D=25 and ORIGINAL metric construction.
% Saved Simulink limiter remains 1.2 pu. No controller/metric optimization.
clearvars;clc;close all;format long g
root=fileparts(mfilename('fullpath'));out=fullfile(root,'current_limit_sensitivity_SCR5_XR10');
if ~isfolder(out),mkdir(out);end
f=fopen(fullfile(out,'analysis_log.txt'),'w');fclose(f);diary(fullfile(out,'analysis_log.txt'));guard=onCleanup(@()diary('off'));
a=load(fullfile(root,'baseline_SCR5_XR10','SCR5_XR10_baseline_results.mat'),'results');base=a.results;p=base.parameters;n=base.network;p.t=n.Vth/p.U;p.Z=n.Zeq;p.phi=angle(n.Vth);p.Vmax=1.3;p.AD=.2;
assert(p.H==3.4&&p.K==20&&p.Kd==5&&p.D==25);
ratings=[1.2 1.3 1.4];gammaGrid=linspace(-180,180,180001).';runs=cell(3,1);summary=table();allIntervals=table();allBounds=table();nominals=table();
for z=1:numel(ratings)
 p.Imax=ratings(z);d=electrical(deg2rad(gammaGrid),p.U,p,n,0);sw=crossings(gammaGrid,d.Ireq-p.Imax,@(x)electrical(deg2rad(x),p.U,p,n,0).Ireq-p.Imax);cuts=[-180;sw;180];scan=table();I=struct('stable',zeros(0,2),'unstable',zeros(0,2),'certified',zeros(0,2),'mathCertified',zeros(0,2),'voltageInvalid',zeros(0,2),'inactive',zeros(0,2),'active',zeros(0,2));kinds=fieldnames(I);modes=[];
 for k=1:numel(cuts)-1
  lo=cuts(k);hi=cuts(k+1);mode=1+double(electrical(deg2rad((lo+hi)/2),p.U,p,n,0).Ireq>p.Imax);gg=unique([lo;gammaGrid(gammaGrid>lo&gammaGrid<hi);hi]);t=point_rows(gg,p,n,mode);scan=[scan;t(2:end-1,:)]; %#ok<AGROW>
  for j=1:numel(kinds),kind=kinds{j};ii=interval_runs(gg,t,p,n,mode,kind);I.(kind)=[I.(kind);ii];if strcmp(kind,'certified'),modes=[modes;repmat(mode,size(ii,1),1)];end;end %#ok<AGROW>
 end
 vb=crossings(gammaGrid,abs(d.u)-p.Vmax,@(x)abs(electrical(deg2rad(x),p.U,p,n,0).u)-p.Vmax);v=validate_all(p,n,base,sw,I);bt=boundary_table(sw,vb,p,n);bt.Imax=repmat(p.Imax,height(bt),1);allBounds=[allBounds;bt];reg=region_table(I.certified,modes,sw,vb,p,n);nom=point_rows(rad2deg(p.gammaNominal),p,n,0);nom.Imax=p.Imax;nominals=[nominals;nom]; %#ok<AGROW>
 for j=1:numel(kinds),ii=I.(kinds{j});allIntervals=[allIntervals;table(repmat(p.Imax,size(ii,1),1),repmat(string(kinds{j}),size(ii,1),1),ii(:,1),ii(:,2),'VariableNames',{'Imax','classification','lower_deg','upper_deg'})];end %#ok<AGROW>
 positiveSwitch=min(sw(sw>0));stablePositive=intersect_width([0 180],I.stable);certPositive=intersect_width([0 180],I.certified);mathPositive=intersect_width([0 180],I.mathCertified);activeStable=0;activeCert=0;
 for k=1:size(I.active,1),aa=[max(0,I.active(k,1)),I.active(k,2)];if aa(2)>aa(1),activeStable=activeStable+intersect_width(aa,I.stable);activeCert=activeCert+intersect_width(aa,I.certified);end;end
 at=electrical(deg2rad(positiveSwitch),p.U,p,n,0);left=electrical(deg2rad(positiveSwitch),p.U,p,n,1);right=electrical(deg2rad(positiveSwitch),p.U,p,n,2);
 summary=[summary;table(p.Imax,positiveSwitch,stablePositive,activeStable,certPositive,mathPositive,activeCert,size(I.certified,1),max([abs(d.u);abs(electrical(deg2rad(sw),p.U,p,n,0).u)]),abs(at.u),left.Pg,right.Pg,p.Imax-nom.Iconv,'VariableNames',{'Imax','positiveSwitch_deg','positiveLocalStableWidth_deg','positiveActiveStableWidth_deg','positiveCertifiedWidth_deg','positiveMathematicalCertificateWidth_deg','positiveActiveCertifiedWidth_deg','certifiedIntervalsFullPeriod','maximumRequiredVoltage_pu','voltageAtSwitch_pu','slopeInactiveAtSwitch','slopeActiveAtSwitch','nominalCurrentMargin_pu'})]; %#ok<AGROW>
 runs{z}=struct('Imax',p.Imax,'parameters',p,'scan',scan,'intervals',I,'switches_deg',sw,'voltageBoundaries_deg',vb,'regions',reg,'nominal',nom,'validation',v);
 fprintf('\nImax=%.1f pu | current switches:\n',p.Imax);disp(sw);fprintf('Local stable intervals:\n');disp(I.stable);fprintf('Certified intervals (voltage valid):\n');disp(I.certified);fprintf('Mathematical certificates ignoring voltage capability:\n');disp(I.mathCertified);fprintf('Voltage-invalid intervals:\n');disp(I.voltageInvalid);
end
writetable(summary,fullfile(out,'Imax_comparison.csv'));writetable(allIntervals,fullfile(out,'all_intervals.csv'));writetable(allBounds,fullfile(out,'switching_and_voltage_boundaries.csv'));writetable(nominals,fullfile(out,'nominal_points.csv'));
fig=figure('Color','w','Position',[50 40 1400 1100]);tl=tiledlayout(fig,4,1,'TileSpacing','compact');colors=lines(3);
for z=1:3
 s=runs{z}.scan;sel=find(s.gamma0_deg>=0 & s.gamma0_deg<=40);sel=sel(1:5:end);s=s(sel,:);ys={s.Pe,s.maxRealEigJ,s.lambda_max_SymF,s.u_required};for k=1:4,ax=nexttile(tl,k);hold(ax,'on');plot(ax,s.gamma0_deg,ys{k},'Color',colors(z,:),'LineWidth',1.4,'DisplayName',sprintf('Imax=%.1f',ratings(z)));xline(ax,runs{z}.switches_deg(end),'--','Color',colors(z,:),'HandleVisibility','off');grid(ax,'on');xlim(ax,[0 40]);end
end
labels={'Grid Pe (pu)','max Re eig(J)','max eig(Sym(F))','Required command voltage (pu)'};for k=1:4,ax=nexttile(tl,k);ylabel(ax,labels{k});xlabel(ax,'Voltage-target gamma (deg)');if k==2||k==3,yline(ax,0,'k:','HandleVisibility','off');end;end
ax=nexttile(tl,1);legend(ax,'Location','best');ax=nexttile(tl,3);ylim(ax,[-.5 10]);ax=nexttile(tl,4);yline(ax,p.Vmax,'r:','1.3 pu capability','HandleVisibility','off');title(tl,{'SCR5/XR10 ideal radial-current benchmark: current-rating sensitivity','H=3.4, K=20, Kd=5, D=25; unchanged metric | dashed: current switches | Sym(F) range clipped'});save_plot(fig,out,'Imax_positive_angle_comparison');
fig=figure('Color','w','Position',[50 50 1350 800]);ax=axes(fig);hold(ax,'on');yt=[];yl={};
for z=1:3
 sets={runs{z}.intervals.stable,runs{z}.intervals.certified,runs{z}.intervals.voltageInvalid};labs={'Locally stable','Contraction-certified','Voltage-invalid'};
 for k=1:3,y=(z-1)*3+k;yt(end+1)=y;yl{end+1}=sprintf('Imax %.1f: %s',ratings(z),labs{k});ii=sets{k};for j=1:size(ii,1),plot(ax,ii(j,:),[y y],'LineWidth',9,'Color',colors(z,:));end;if isempty(ii),text(ax,0,y,'None','HorizontalAlignment','center');end;end
end
yticks(ax,yt);yticklabels(ax,yl);ylim(ax,[.5 9.5]);xlim(ax,[-40 40]);grid(ax,'on');xlabel(ax,'Target gamma (deg); display zoom -40..40, full-period data in MAT/CSV');title(ax,{'SCR5/XR10 current-rating comparison | fixed H=3.4, D=25 and metric','Smooth branches only; no global switched-system certificate'});save_plot(fig,out,'Imax_interval_comparison');
results=struct('summary',summary,'runs',{runs},'allIntervals',allIntervals,'allBoundaries',allBounds,'nominals',nominals,'baselineRegions',base.regions_deg);save(fullfile(out,'Imax_sensitivity_results.mat'),'results','-v7.3');
fid=fopen(fullfile(out,'report.md'),'w','n','UTF-8');fprintf(fid,'# Imax sensitivity: SCR5/XR10\n\nNetwork-consistent two-state ideal radial-current benchmark. Hypothetical analytical ratings; saved Simulink setting is not changed. H=3.4, K=20, Kd=5, D=25, E/network/metric unchanged. No metric optimization, sag or trajectory study is included. Full period -180..180 deg scanned at 0.002 deg followed by root refinement. Exact switches are excluded from ordinary smooth certificates.\n\n');
fprintf(fid,'| Imax pu | Positive switch deg | Positive local-stable width deg | Positive certified width deg | Full-period certified intervals | Max required voltage pu |\n|---|---:|---:|---:|---:|---:|\n');for k=1:height(summary),s=summary(k,:);fprintf(fid,'| %.1f | %.9f | %.9f | %.9f | %d | %.9f |\n',s.Imax,s.positiveSwitch_deg,s.positiveLocalStableWidth_deg,s.positiveCertifiedWidth_deg,s.certifiedIntervalsFullPeriod,s.maximumRequiredVoltage_pu);end
fprintf(fid,'\nIncreasing Imax delays saturation and enlarges the positive-angle locally stable interval in these tests. Certification is evaluated separately, both with and without the 1.3-pu voltage-capability check. See all_intervals.csv for every connected interval and switching_and_voltage_boundaries.csv for refined switching/voltage boundaries. The voltage limit is not simulated; any exceedance invalidates the current-only realization assumption.\n\nThe unchanged baseline R1 starts at %.9f deg and requires %.9f pu converter current there. Raising Imax to 1.3/1.4 alone need not expose this R1: the unchanged metric must still be tested on each branch. No current rating is recommended as a hardware change by this analytical experiment.\n\n',base.regions_deg(1,1),electrical(deg2rad(base.regions_deg(1,1)),p.U,p,n,1).Ireq);
fprintf(fid,'## Validation\n\n');for k=1:3,fn=fieldnames(runs{k}.validation);fprintf(fid,'Imax=%.1f\n\n',ratings(k));for j=1:numel(fn),fprintf(fid,'- %s: %.9g\n',fn{j},runs{k}.validation.(fn{j}));end;fprintf(fid,'\n');end;fclose(fid);
disp(summary);fprintf('Completed separate Imax sensitivity. No existing model or analysis modified.\n');diary off
function w=intersect_width(a,I)
w=0;for k=1:size(I,1),w=w+max(0,min(a(2),I(k,2))-max(a(1),I(k,1)));end
end
function save_plot(fig,out,name)
axs=findall(fig,'Type','axes');for k=1:numel(axs),axs(k).Toolbar.Visible='off';end;drawnow;exportgraphics(fig,fullfile(out,[name,'.png']),'Resolution',180);savefig(fig,fullfile(out,[name,'.fig']));close(fig);
end
function d=electrical(g,U,p,n,mode)
% mode 0 automatic; 1 force inactive; 2 force active (one-sided derivatives).
v=p.E*exp(1i*(g+p.phi));q=(v-p.t*U)/p.Z;qg=1i*v/p.Z;qU=-p.t/p.Z+zeros(size(q));r=abs(q);
active=r>p.Imax;if mode==1,active=false(size(q));elseif mode==2,active=true(size(q));end
ic=q;icg=qg;icU=qU;icgU=zeros(size(q));
if any(active(:))
 a=active;qa=q(a);ra=r(a);ga=qg(a);ua=qU(a);rg=real(conj(qa).*ga)./ra;rU=real(conj(qa).*ua)./ra;rgU=real(conj(ua).*ga)./ra-rg.*rU./ra;
 ic(a)=p.Imax*qa./ra;icg(a)=p.Imax*(ga./ra-qa.*rg./ra.^2);icU(a)=p.Imax*(ua./ra-qa.*rU./ra.^2);
 icgU(a)=p.Imax*(-ga.*rU./ra.^2-ua.*rg./ra.^2-qa.*rgU./ra.^2+2*qa.*rg.*rU./ra.^3);
end
vc=p.t*U+p.Z*ic;vcg=p.Z*icg;vcU=p.t+p.Z*icU;vcgU=p.Z*icgU;
% On the inactive branch set exact target to avoid cancellation roundoff.
vc(~active)=v(~active);vcg(~active)=1i*v(~active);vcU(~active)=0;vcgU(~active)=0;
vg=n.vgE*vc+n.vgU*U;ig=n.igE*vc+n.igU*U;
vgg=n.vgE*vcg;igg=n.igE*vcg;vgU=n.vgE*vcU+n.vgU;igU=n.igE*vcU+n.igU;
vgm=n.vgE*vcgU;igm=n.igE*vcgU;
Pe=real(vg.*conj(ig));Pg=real(vgg.*conj(ig)+vg.*conj(igg));PU=real(vgU.*conj(ig)+vg.*conj(igU));PgU=real(vgm.*conj(ig)+vgg.*conj(igU)+vgU.*conj(igg)+vg.*conj(igm));
vb=vc+n.Zconverter*ic;u=vb+p.AD*ic;
Uvec=U+zeros(size(vc));
d=struct('Pe',Pe,'Pg',Pg,'PU',PU,'PgU',PgU,'q',q,'Ireq',r,'Iconv',ic,'Vc',vc,'Vg',vg,'Ig',ig,'u',u,'Vbridge',vb,'active',active,'Ic_g',icg,'Ic_U',icU,'nodes',n.nodeCoefficients*[vc(:).';Uvec(:).']);
end
function t=point_rows(gd,p,n,mode)
g=deg2rad(gd(:));d=electrical(g,p.U,p,n,mode);beta=-d.PU/p.mass;c=-d.PgU/p.mass;r=c./beta;a=-d.Pg/p.mass;b=p.D/p.mass;off=(1+a+b*r-r.^2)/2;
lm=-b/2+hypot(b/2-r,off);mn=-b/2-hypot(b/2-r,off);mm=(2*beta.^2+c.^2+hypot(c.^2,2*beta.*c))/2;mi=beta.^4./mm;
e1=(-b+sqrt(complex(b*b+4*a)))/2;e2=(-b-sqrt(complex(b*b+4*a)))/2;local=max(real(e1),real(e2));singular=abs(beta)<=1e-12;lm(singular)=NaN;mn(singular)=NaN;
valid=abs(d.u)<=p.Vmax;onSwitch=abs(d.Ireq-p.Imax)<1e-10;math=local<-p.tol & mi>p.tol & lm<-p.tol & ~singular;cert=math & valid & ~onSwitch;
t=table(gd(:),d.Pe,d.Pe,d.Ireq,abs(d.Iconv),abs(d.Vc),rad2deg(angle(d.Vc)),abs(d.u),d.active,valid,d.Pg,d.PU,d.PgU,real(e1),imag(e1),real(e2),imag(e2),local,local<-p.tol,beta,c,beta.^2,mi,mm,lm,mn,math,cert,onSwitch,zeros(size(g)), ...
 'VariableNames',{'gamma0_deg','Pstar','Pe','Ireq','Iconv','Vc_abs','Vc_angle_deg','u_required','limiter_active','voltage_capability_valid','Pe_gamma','Pe_U','Pe_gammaU','eig1_real','eig1_imag','eig2_real','eig2_imag','maxRealEigJ','local_stable','beta','c','detTheta','lambdaMinM','lambdaMaxM','lambda_max_SymF','lambda_min_SymF','mathematical_branch_certificate','contraction_certified','on_switch','equilibriumResidual'});
end
function s=matrix_point(g,p,n,mode)
d=electrical(g,p.U,p,n,mode);beta=-d.PU/p.mass;c=-d.PgU/p.mass;A=[0 1;-d.Pg/p.mass -p.D/p.mass];T=[-beta 0;-c -beta];M=T'*T;F=nan(2);S=F;if abs(beta)>1e-12,F=T*A/T;S=(F+F')/2;end
s=struct('gamma',g,'Pstar',d.Pe,'A',A,'Theta',T,'M',M,'F',F,'SymF',S,'beta',beta,'c',c);
end
function dx=rhs(~,x,p,n,ps,U,mode)
d=electrical(x(1),U,p,n,mode);dx=[x(2)-1;((ps-d.Pe)/x(2)-p.D*(x(2)-1))/p.mass];
end
function A=jac(g,nu,p,n,ps,U,mode)
d=electrical(g,U,p,n,mode);A=[0 1;-d.Pg/(p.mass*nu),(-(ps-d.Pe)/nu^2-p.D)/p.mass];
end
function roots=crossings(g,y,fun)
idx=find(y(1:end-1).*y(2:end)<0);roots=g(y==0);for k=idx(:).',roots(end+1,1)=fzero(fun,g(k+[0 1]),optimset('TolX',1e-11));end;roots=unique(sort(roots));
end
function m=margin(t,p,kind)
switch kind
 case 'stable',m=-t.maxRealEigJ-p.tol;
 case 'unstable',m=t.maxRealEigJ-p.tol;
 case 'certified',m=min([t.lambdaMinM-p.tol,-t.lambda_max_SymF-p.tol,-t.maxRealEigJ-p.tol,p.Vmax-t.u_required],[],2);m(abs(t.beta)<=1e-12)=-Inf;
 case 'mathCertified',m=min([t.lambdaMinM-p.tol,-t.lambda_max_SymF-p.tol,-t.maxRealEigJ-p.tol],[],2);m(abs(t.beta)<=1e-12)=-Inf;
 case 'voltageInvalid',m=t.u_required-p.Vmax;
 case 'active',m=t.Ireq-p.Imax;
 case 'inactive',m=p.Imax-t.Ireq;
end
end
function I=interval_runs(g,t,p,n,mode,kind)
if strcmp(kind,'active') || strcmp(kind,'inactive')
 I=zeros(0,2);if (strcmp(kind,'active') && mode==2)||(strcmp(kind,'inactive') && mode==1),I=[g(1) g(end)];end;return
end
y=margin(t,p,kind);ok=y>0;beg=find(diff([false;ok])==1);en=find(diff([ok;false])==-1);I=zeros(numel(beg),2);
for k=1:numel(beg)
 a=g(beg(k));b=g(en(k));
 if beg(k)>1,a=fzero(@(x)margin(point_rows(x,p,n,mode),p,kind),g(beg(k)+[-1 0]),optimset('TolX',1e-11));end
 if en(k)<numel(g),b=fzero(@(x)margin(point_rows(x,p,n,mode),p,kind),g(en(k)+[0 1]),optimset('TolX',1e-11));end
 I(k,:)=[a b];
end
I=I(diff(I,1,2)>1e-8,:);
end
function t=boundary_table(sw,vb,p,n)
t=table();for j=1:2,if j==1,rr=sw;kind="current-limiter switching";else,rr=vb;kind="voltage-capability limit";end
 for g=rr(:).'
  a=electrical(deg2rad(g-1e-5),p.U,p,n,0);b=electrical(deg2rad(g+1e-5),p.U,p,n,0);d=electrical(deg2rad(g),p.U,p,n,0);
  before="inactive";after="inactive";if a.active,before="active";end;if b.active,after="active";end
  l=electrical(deg2rad(g),p.U,p,n,1);r=electrical(deg2rad(g),p.U,p,n,2);
  t=[t;table(g,kind,before,after,d.Ireq,abs(d.u),l.Pg,r.Pg,'VariableNames',{'gamma_deg','boundary_type','branch_before','branch_after','Ireq','u_required','inactive_side_Pe_gamma','active_side_Pe_gamma'})]; %#ok<AGROW>
 end
end
end
function cause=boundary_cause(g,sw,vb,p,n,mode)
if abs(abs(g)-180)<1e-7,cause="domain boundary";elseif any(abs(sw-g)<1e-6),cause="current-limiter switching";elseif any(abs(vb-g)<1e-5),cause="voltage-capability limit";else
 t=point_rows(g,p,n,mode);if t.lambdaMinM<2*p.tol,cause="metric positive-definiteness threshold / Theta singularity vicinity";elseif abs(t.maxRealEigJ)<2e-8,cause="local-stability loss";else,cause="contraction condition lambda_max(SymF) = -tol";end
end
end
function t=region_table(I,modes,sw,vb,p,n)
t=table('Size',[0 9],'VariableTypes',{'double','double','double','logical','string','logical','logical','string','string'},'VariableNames',{'regionID','lower_gamma_deg','upper_gamma_deg','local_stable','limiter_mode','contraction_certified','voltage_capability_valid','lower_cause','upper_cause'});
for k=1:size(I,1),mode="inactive";if modes(k)==2,mode="active";end;t(end+1,:)={k,I(k,1),I(k,2),true,mode,true,true,boundary_cause(I(k,1),sw,vb,p,n,modes(k)),boundary_cause(I(k,2),sw,vb,p,n,modes(k))};end
end
function v=validate_all(p,n,base,sw,I)
v=struct('inactiveBaseline',0,'powerContinuity',0,'voltageContinuity',0,'currentLimit',0,'Pe_gamma',0,'Pe_U',0,'Pe_gammaU',0,'projectionDerivative',0,'Jacobian',0,'matrixFormula',0,'matrixSymmetry',0,'networkResidual',0,'equilibriumResidual',0,'nominalMatch',0);
h=2e-5;gs=linspace(-179,179,121);gs=gs(all(abs(gs(:)-sw(:).')>.05,2));
for gd=gs(:).'
 g=deg2rad(gd);d=electrical(g,p.U,p,n,0);mode=1+double(d.active);a=electrical(g+h,p.U,p,n,mode);b=electrical(g-h,p.U,p,n,mode);c=electrical(g,p.U+h,p,n,mode);e=electrical(g,p.U-h,p,n,mode);
 v.Pe_gamma=max(v.Pe_gamma,abs((a.Pe-b.Pe)/(2*h)-d.Pg));v.Pe_U=max(v.Pe_U,abs((c.Pe-e.Pe)/(2*h)-d.PU));v.Pe_gammaU=max(v.Pe_gammaU,abs((a.PU-b.PU)/(2*h)-d.PgU));
 if ~d.active
  co=p.pc;pe=co(1)*p.E^2+co(2)*p.U^2+p.E*p.U*(co(3)*cos(g)+co(4)*sin(g));v.inactiveBaseline=max(v.inactiveBaseline,abs(d.Pe-pe));
 else
  v.currentLimit=max(v.currentLimit,abs(abs(d.Iconv)-p.Imax));q=[real(d.q);imag(d.q)];D=p.Imax/norm(q)*(eye(2)-q*q'/norm(q)^2);Dn=zeros(2);for k=1:2,z=zeros(2,1);z(k)=h;Dn(:,k)=p.Imax*((q+z)/norm(q+z)-(q-z)/norm(q-z))/(2*h);end;v.projectionDerivative=max(v.projectionDerivative,norm(D-Dn,inf));
 end
 v.networkResidual=max(v.networkResidual,abs(n.icE*d.Vc+n.icU*p.U-d.Iconv));
 for nu=[.85 1 1.1]
  x=[g;nu];J=jac(g,nu,p,n,d.Pe,p.U,mode);Jn=zeros(2);for k=1:2,z=zeros(2,1);z(k)=h;Jn(:,k)=(rhs(0,x+z,p,n,d.Pe,p.U,mode)-rhs(0,x-z,p,n,d.Pe,p.U,mode))/(2*h);end;v.Jacobian=max(v.Jacobian,norm(J-Jn,inf));
 end
 m=matrix_point(g,p,n,mode);t=point_rows(gd,p,n,mode);if abs(m.beta)>1e-8,v.matrixFormula=max(v.matrixFormula,abs(max(eig(m.SymF))-t.lambda_max_SymF)/max(1,abs(t.lambda_max_SymF)));v.matrixSymmetry=max(v.matrixSymmetry,norm(m.SymF-m.SymF'));assert(min(eig(m.M))>0);end
 v.equilibriumResidual=max(v.equilibriumResidual,norm(rhs(0,[g;1],p,n,d.Pe,p.U,mode),inf));
end
for gd=sw(:).',a=electrical(deg2rad(gd),p.U,p,n,1);b=electrical(deg2rad(gd),p.U,p,n,2);v.powerContinuity=max(v.powerContinuity,abs(a.Pe-b.Pe));v.voltageContinuity=max(v.voltageContinuity,abs(a.Vc-b.Vc));end
for k=1:size(I.certified,1),gg=linspace(I.certified(k,1)+1e-6,I.certified(k,2)-1e-6,101);t=point_rows(gg,p,n,0);assert(all(t.contraction_certified));end
nom=point_rows(rad2deg(p.gammaNominal),p,n,0);bn=base.nominalEquilibria(base.nominalEquilibria.isSelected,:);v.nominalMatch=max(abs([nom.Pe-1,nom.eig1_real-bn.eig1_real,nom.eig2_real-bn.eig2_real,nom.lambda_max_SymF-bn.lambdaMaxSymF]));assert(~nom.limiter_active&&nom.voltage_capability_valid);
assert(v.inactiveBaseline<1e-10&&v.powerContinuity<1e-10&&v.voltageContinuity<1e-10&&v.currentLimit<1e-12&&v.networkResidual<1e-11);
assert(v.Pe_gamma<1e-6&&v.Pe_U<1e-6&&v.Pe_gammaU<1e-5&&v.Jacobian<1e-6&&v.matrixFormula<1e-10&&v.projectionDerivative<1e-6&&v.nominalMatch<1e-9);
fprintf('\nNUMERICAL VALIDATION PASSED\n');disp(v);
end

