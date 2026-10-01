%% Network-consistent two-state ideal radial-current benchmark
% New SCR5/XR10 network; not an exact PI/freeze controller model.
% All angles are voltage TARGET angles; actual Vc differs in saturation.
clearvars;clc;close all;format long g
root=fileparts(mfilename('fullpath'));out=fullfile(root,'current_limiter_SCR5_XR10');if ~isfolder(out),mkdir(out);end
logfile=fullfile(out,'analysis_log.txt');f=fopen(logfile,'w');fclose(f);diary(logfile);guard=onCleanup(@()diary('off'));
r=load(fullfile(root,'baseline_SCR5_XR10','SCR5_XR10_baseline_results.mat'),'results');base=r.results;
p=base.parameters;n=base.network;p.t=n.Vth/p.U;p.Z=n.Zeq;p.phi=angle(n.Vth);
addpath(genpath(p.libraryRoot));load_system(p.modelPath);g='GFMVSG1/Grid_Forming_Converter';
p.Imax=str2double(get_param(g,'InnerCtrl_i_conv_max'));p.Vmax=str2double(get_param(g,'InnerCtrl_v_conv_max'));p.AD=str2double(get_param(g,'InnerCtrl_AD_R_virt_pu'));
assert(p.H==3.4&&p.K==20&&p.Kd==5&&p.D==25&&p.mass==6.8);
assert(p.Imax==1.2&&p.Vmax==1.3&&p.AD==.2);assert(strcmp(get_param(g,'InnerCtrl_i_lim_type'),'dq coordinate system'));
assert(str2double(get_param(g,'OuterCtrl_VSM_H'))==p.H&&str2double(get_param(g,'OuterCtrl_VSM_Kd'))==p.Kd);
frozenNames={'baseline_gfm_contraction.m','gfm_current_limiter_analysis.m','baseline_gfm_contraction_SCR5_XR10.m'};frozen=cellfun(@(x)fileread(fullfile(root,x)),frozenNames,'UniformOutput',false);
grid=linspace(-180,180,360001).';d=electrical(deg2rad(grid),p.U,p,n,0);
switches=crossings(grid,d.Ireq-p.Imax,@(x)electrical(deg2rad(x),p.U,p,n,0).Ireq-p.Imax);
cuts=[-180;switches;180];scan=table();intervals=struct('stable',zeros(0,2),'unstable',zeros(0,2),'certified',zeros(0,2),'mathCertified',zeros(0,2),'voltageInvalid',zeros(0,2),'inactive',zeros(0,2),'active',zeros(0,2));
kinds=fieldnames(intervals);regionModes=[];
for k=1:numel(cuts)-1
 a=cuts(k);b=cuts(k+1);mode=1+double(electrical(deg2rad((a+b)/2),p.U,p,n,0).Ireq>p.Imax);
 gg=unique([a;grid(grid>a&grid<b);b]);t=point_rows(gg,p,n,mode);scan=[scan;t(2:end-1,:)]; %#ok<AGROW>
 for j=1:numel(kinds)
  kind=kinds{j};ii=interval_runs(gg,t,p,n,mode,kind);intervals.(kind)=[intervals.(kind);ii]; %#ok<AGROW>
  if strcmp(kind,'certified'),regionModes=[regionModes;repmat(mode,size(ii,1),1)];end %#ok<AGROW>
 end
end
scan=sortrows(scan,'gamma0_deg');voltageBounds=crossings(grid,abs(d.u)-p.Vmax,@(x)abs(electrical(deg2rad(x),p.U,p,n,0).u)-p.Vmax);
nominal=point_rows(rad2deg(p.gammaNominal),p,n,0);
validation=validate_all(p,n,base,switches,intervals);
boundaries=boundary_table(switches,voltageBounds,p,n);regions=region_table(intervals.certified,regionModes,switches,voltageBounds,p,n);
[rep,labels]=representatives(intervals,p,n);[rep,labels]=add_active_stable(rep,labels,intervals);points=unique([rep;rad2deg(p.gammaNominal);mean(base.regions_deg,2);reshape(base.regions_deg,[],1);switches-.001;switches+.001]);points=point_rows(points,p,n,0);
matrices=arrayfun(@(v)matrix_point(deg2rad(v),p,n,0),points.gamma0_deg);
coverage=coverage_table(base.regions_deg,intervals);comparison=comparison_table(base,nominal,intervals.certified);
network=table(p.E,real(n.Vth),imag(n.Vth),real(p.Z),imag(p.Z),abs(p.Z),rad2deg(angle(p.Z)),p.U,p.Imax,p.Vmax,p.H,p.K,p.Kd,p.D,'VariableNames',{'E','Vth_real','Vth_imag','Req','Xeq','absZeq','angle_Zeq_deg','physicalBusU','Imax','voltage_limit','H','K','Kd','D_eff'});
writetable(network,fullfile(out,'current_limiter_network_parameters.csv'));writetable(boundaries,fullfile(out,'current_limiter_switching_boundaries.csv'));writetable(points,fullfile(out,'current_limiter_operating_points.csv'));writetable(regions,fullfile(out,'current_limiter_regions.csv'));writetable(comparison,fullfile(out,'baseline_vs_current_limiter_comparison.csv'));writetable(coverage,fullfile(out,'baseline_region_coverage.csv'));writetable(scan,fullfile(out,'current_limiter_dense_scan.csv'));
allIntervals=table();for j=1:numel(kinds),ii=intervals.(kinds{j});allIntervals=[allIntervals;table(repmat(string(kinds{j}),size(ii,1),1),ii(:,1),ii(:,2),'VariableNames',{'classification','lower_deg','upper_deg'})];end %#ok<AGROW>
writetable(allIntervals,fullfile(out,'all_classification_intervals.csv'));
figureData=comparison_figures(scan,base,switches,voltageBounds,intervals,p,n,out);
[trajectories,trajectorySummary]=trajectory_study(rep,labels,p,n,out);
% Only nominal is a shared feasible equilibrium. Baseline R1/R2 Pstar cannot
% be used as initial current-limited equilibria. Sag comparison is nominal only.
[sagRuns,sagSummary]=nominal_sags(p,n,base,out);
results=struct('label','Network-consistent two-state ideal radial-current benchmark','parameters',p,'network',n,'baselineRegions',base.regions_deg,'baselineNominal',base.nominalEquilibria,'scan',scan,'switches_deg',switches,'voltageBoundaries_deg',voltageBounds,'intervals',intervals,'regions',regions,'boundaries',boundaries,'nominal',nominal,'operatingPoints',points,'matrices',matrices,'coverage',coverage,'comparison',comparison,'validation',validation,'figureData',figureData,'trajectories',{trajectories},'trajectorySummary',trajectorySummary,'sagRuns',{sagRuns},'sagSummary',sagSummary);
save(fullfile(out,'current_limiter_SCR5_XR10_results.mat'),'results','-v7.3');report(results,base,out);
for k=1:numel(frozenNames),assert(strcmp(frozen{k},fileread(fullfile(root,frozenNames{k}))));end
fprintf('\nFinished analytical benchmark. No EMT validation run. Outputs: %s\n',out);diary off

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
function [g,labels]=representatives(I,p,n)
g=rad2deg(p.gammaNominal);labels="Nominal";
sets={I.certified,I.active,I.stable,I.unstable};names=["Certified","Limiter active","Stable uncertified","Locally unstable"];
for k=1:4
 for j=1:size(sets{k},1)
  ii=sets{k}(j,:);if ii(2)<=0,continue;end;gg=mean([max(0,ii(1)),ii(2)]);if k==4,gg=min(max(0,ii(1))+15,gg);end;t=point_rows(gg,p,n,0);
  if k==3 && t.contraction_certified,continue;end
  if any(abs(g-gg)<.1),continue;end
  g(end+1,1)=gg;labels(end+1,1)=names(k);if ~t.voltage_capability_valid,labels(end)=labels(end)+" (voltage invalid)";end %#ok<AGROW>
  if k>1,break;end
 end
end
end
function [g,labels]=add_active_stable(g,labels,I)
for j=1:size(I.stable,1)
 for k=1:size(I.active,1)
 a=max(I.stable(j,1),I.active(k,1));b=min(I.stable(j,2),I.active(k,2));
 if b-a>1e-5,g(end+1,1)=(a+b)/2;labels(end+1,1)="Stable active (import)";return;end
 end
end
end
function w=overlap_width(a,I)
w=0;for k=1:size(I,1),w=w+max(0,min(a(2),I(k,2))-max(a(1),I(k,1)));end
end
function t=coverage_table(B,I)
t=table();for k=1:size(B,1),a=B(k,:);width=diff(a);cf=overlap_width(a,I.inactive);act=overlap_width(a,I.active);invalid=overlap_width(a,I.voltageInvalid);st=overlap_width(a,I.stable);cert=overlap_width(a,I.certified);
t=[t;table(k,a(1),a(2),width,cf,act,width-invalid,invalid,st,cert,width-cert,'VariableNames',{'baselineRegion','lower_deg','upper_deg','width_deg','currentFeasibleWidth_deg','currentLimitedWidth_deg','limitedVoltageValidWidth_deg','limitedVoltageInvalidWidth_deg','limitedLocalStableWidth_deg','limitedCertifiedWidth_deg','lostCertifiedWidth_deg'})];end
end
function t=comparison_table(base,nom,I)
B=base.regions_deg;t=table("nominal angle",rad2deg(base.parameters.gammaNominal),nom.gamma0_deg,'VariableNames',{'quantity','baseline','limited'});t=[t;table("nominal Pstar",1,nom.Pstar,'VariableNames',t.Properties.VariableNames)];
for k=1:max(size(B,1),size(I,1)),for j=1:2,a=NaN;b=NaN;if k<=size(B,1),a=B(k,j);end;if k<=size(I,1),b=I(k,j);end;t=[t;table("region "+k+" endpoint "+j,a,b,'VariableNames',t.Properties.VariableNames)];end;end
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
function save_plot(fig,out,name)
axs=findall(fig,'Type','axes');for aa=1:numel(axs),axs(aa).Toolbar.Visible='off';end;drawnow;
exportgraphics(fig,fullfile(out,[name,'.png']),'Resolution',180);savefig(fig,fullfile(out,[name,'.fig']));close(fig);
end
function label_plot(ax,txt)
grid(ax,'on');xlabel(ax,'Voltage-target gamma (deg)');title(ax,{'Network-consistent two-state ideal radial-current benchmark','SCR=5, X/R=10 | H=3.4, K=20, Kd=5 | '+string(txt)},'FontSize',11);
end
function mark_boundaries(ax,sw,vb,B)
for x=sw(:).',xline(ax,x,'k--','HandleVisibility','off');end
for x=vb(:).',xline(ax,x,'m:','HandleVisibility','off');end
for x=B(:).',xline(ax,x,'b:','HandleVisibility','off');end
end
function figureData=comparison_figures(scan,base,sw,vb,I,p,n,out)
g=unique([linspace(-180,180,12001).';sw]);d=electrical(deg2rad(g),p.U,p,n,0);u=electrical(deg2rad(g),p.U,p,n,1);t=point_rows(g,p,n,0);bt=point_rows(g,p,n,1);near=false(size(g));for x=sw(:).',near=near|abs(g-x)<1e-8;end
vals={ [u.Pe,d.Pe],[d.Ireq,abs(d.Iconv)],[abs(u.u),abs(d.u)],[u.Pg,d.Pg],[bt.maxRealEigJ,t.maxRealEigJ],[bt.lambda_max_SymF,t.lambda_max_SymF]};
figureData=struct('gamma_deg',g,'baselineElectrical',u,'limitedElectrical',d,'baselineRows',bt,'limitedRows',t,'curves',{vals});
names={'A_power_angle','B_converter_current','C_voltage_capability','D_power_derivative','E_local_eigenvalue','F_contraction_diagnostic'};
ylabs={'Grid-measured Pe (pu)','Converter current magnitude (pu)','Pre-active-damping voltage magnitude (pu)','dPe/dgamma (pu/rad)','max Re eig(J)','max eig(Sym(F))'};
for k=1:6
 fig=figure('Color','w','Position',[60 60 1400 650]);ax=axes(fig);hold(ax,'on');y=vals{k};if k>=4,y(near,2)=NaN;end
 plot(ax,g,y(:,1),'Color',[.4 .4 .4],'LineWidth',1.5,'DisplayName','Limiter-free baseline');plot(ax,g,y(:,2),'b','LineWidth',1.6,'DisplayName','Ideal current-limited');
 if k==2,ax.Children(2).DisplayName='Unlimited current demand';ax.Children(1).DisplayName='Ideal delivered current';yline(ax,p.Imax,'r--','Imax=1.2','HandleVisibility','off');end
 if k==3,yline(ax,p.Vmax,'r--','Voltage capability=1.3','HandleVisibility','off');bad=abs(d.u)>p.Vmax;plot(ax,g(bad),abs(d.u(bad)),'r.','DisplayName','Outside current-only validity');end
 if k==4,for gg=sw(:).',dd1=electrical(deg2rad(gg),p.U,p,n,1);dd2=electrical(deg2rad(gg),p.U,p,n,2);plot(ax,[gg gg],[dd1.Pg dd2.Pg],'ko','HandleVisibility','off');end;end
 if k>=4,yline(ax,0,'k:','HandleVisibility','off');end
 if k==6,ylim(ax,[-2 10]);for j=1:size(I.certified,1),plot(ax,I.certified(j,:),[-1 -1],'g','LineWidth',6,'HandleVisibility','off');end;for j=1:size(base.regions_deg,1),plot(ax,base.regions_deg(j,:),[-1.5 -1.5],'Color',[.4 .4 .4],'LineWidth',6,'HandleVisibility','off');end;end
 mark_boundaries(ax,sw,vb,base.regions_deg);xline(ax,rad2deg(p.gammaNominal),'Color',[0 .5 0],'Label','Nominal','HandleVisibility','off');xlim(ax,[-180 180]);ylabel(ax,ylabs{k});label_plot(ax,ylabs{k});legend(ax,'Location','best');save_plot(fig,out,names{k});
end
fig=figure('Color','w','Position',[50 50 1400 650]);ax=axes(fig);hold(ax,'on');sets={base.regions_deg,I.certified,I.active,I.voltageInvalid};colors=[.4 .4 .4;0 .5 0;.9 .5 0;.7 0 .7];labels={'Baseline certified','Limited branchwise certified','Limiter active','Voltage-invalid'};
for j=1:4,ii=sets{j};for k=1:size(ii,1),plot(ax,ii(k,:),[j j],'LineWidth',12,'Color',colors(j,:));end;if isempty(ii),text(ax,0,j,'None','HorizontalAlignment','center');end;end
yticks(ax,1:4);yticklabels(ax,labels);ylim(ax,[.5 4.5]);xlim(ax,[-180 180]);label_plot(ax,'Open smooth-branch intervals; no global switched-system certificate');save_plot(fig,out,'G_region_comparison');
end
function [value,terminal,direction]=stop_event(~,x)
value=x(2)-.5;terminal=1;direction=-1;
end
function tr=integrate(g,ps,p,n,mode,offset,U,T)
opt=odeset('RelTol',1e-9,'AbsTol',1e-11,'MaxStep',max(.05,T/5000),'Events',@stop_event);sol=ode15s(@(t,x)rhs(t,x,p,n,ps,U,mode),[0 T],[g+deg2rad(offset);1],opt);tt=unique([linspace(0,min(20,sol.x(end)),1001),linspace(0,sol.x(end),5001)]).';xx=deval(sol,tt).';dd=electrical(xx(:,1),U,p,n,mode);
tr=struct('time',tt,'state',xx,'electrical',dd,'mode',mode,'Pstar',ps,'equilibriumGamma',g,'terminated',sol.x(end)<T-1e-7);
end
function [runs,summary]=trajectory_study(gd,labels,p,n,out)
runs=cell(numel(gd),1);summary=table();figH=figure('Color','w','Position',[30 40 1700 800]);th=tiledlayout(figH,2,numel(gd));figI=figure('Color','w','Position',[30 40 1700 800]);ti=tiledlayout(figI,2,numel(gd));
for k=1:numel(gd)
 g=deg2rad(gd(k));m=matrix_point(g,p,n,0);b=matrix_point(g,p,n,1);T=max(200,20/max(abs(max(real(eig(m.A)))),1e-5));offset=1;if contains(labels(k),'Stable active'),offset=.1;end;l=integrate(g,m.Pstar,p,n,0,offset,p.U,T);un=integrate(g,b.Pstar,p,n,1,offset,p.U,T);nb=integrate(g,m.Pstar,p,n,0,offset+.01,p.U,T);
 % Same LIMITED reference Theta is applied to both curves. Each model uses
 % its own equilibrium Pstar at the SAME target angle (not same power).
 l.distance=vecnorm(m.Theta*(l.state-[g 1]).',2,1).';un.distance=vecnorm(m.Theta*(un.state-[g 1]).',2,1).';nt=min(l.time(end),nb.time(end));ok=l.time<=nt;xx=interp1(nb.time,nb.state,l.time(ok),'pchip');l.neighborDistance=vecnorm(m.Theta*(xx-l.state(ok,:)).',2,1).';l.neighborTime=l.time(ok);
 valid=all(abs(l.electrical.u)<=p.Vmax+1e-9);tail=l.time>=.9*l.time(end);recovered=all(abs(rad2deg(l.state(tail,1)-g))<.1 & abs(l.state(tail,2)-1)<1e-4)&&~l.terminated;pr=l.electrical.Pe(end)-m.Pstar;
 row=table(labels(k),gd(k),offset,m.Pstar,b.Pstar,rad2deg(l.state(end,1)),l.state(end,2),pr,recovered,any(l.electrical.active),valid,'VariableNames',{'caseLabel','gamma0_deg','perturbation_deg','limitedPstar','baselinePstar','gammaFinal_deg','nuFinal','finalPowerResidual','recoveredOriginal','limiterActivated','voltageValidWholeTrajectory'});summary=[summary;row];runs{k}=struct('limited',l,'baseline',un,'neighbor',nb,'metric',m,'summary',row); %#ok<AGROW>
 ax=nexttile(th,k);plot(ax,un.time,rad2deg(un.state(:,1)),'--',l.time,rad2deg(l.state(:,1)),'LineWidth',1.2);hold(ax,'on');yline(ax,gd(k),'r:');grid(ax,'on');title(ax,{char(labels(k)),sprintf('gamma0=%.4f; valid path=%d',gd(k),valid)},'FontSize',9);ylabel(ax,'Unwrapped gamma (deg)');xlabel(ax,'Normalized time');
 ax=nexttile(th,numel(gd)+k);plot(ax,rad2deg(l.state(:,1)),l.state(:,2)-1);hold(ax,'on');plot(ax,gd(k),0,'ro');plot(ax,gd(k)+offset,0,'ks');plot(ax,rad2deg(l.state(end,1)),l.state(end,2)-1,'gd');grid(ax,'on');xlabel(ax,'Unwrapped target gamma (deg)');ylabel(ax,'nu-1');
 ax=nexttile(ti,k);semilogy(ax,un.time,max(un.distance,1e-16),'--',l.time,max(l.distance,1e-16),'LineWidth',1.2);grid(ax,'on');title(ax,labels(k),'FontSize',9);xlabel(ax,'Normalized time');ylabel(ax,'L: SAME limited-reference metric');
 ax=nexttile(ti,numel(gd)+k);semilogy(ax,l.neighborTime,max(l.neighborDistance,1e-16));grid(ax,'on');xlabel(ax,'Normalized time');ylabel(ax,'Neighbor distance (0.01 deg separation)');
end
title(th,{'Network-consistent two-state ideal radial-current benchmark | SCR5/XR10 | H=3.4, D=25','+1 deg test; narrow stable import branch uses +0.1 deg; dashed=baseline, solid=limited; each uses own equilibrium Pstar'});title(ti,{'SCR5/XR10 ideal radial benchmark: fixed metric distances','Same limited-reference Theta for each pair; no certificate implied by trajectory recovery'});save_plot(figH,out,'H_representative_trajectories');save_plot(figI,out,'I_metric_distance');writetable(summary,fullfile(out,'trajectory_summary.csv'));
end
function [runs,summary]=nominal_sags(p,n,base,out)
% Optional fair same-Pstar comparison only at common nominal equilibrium.
Us=[.5 .3 .1];runs=cell(numel(Us),2);summary=table();T=200;g=p.gammaNominal;ps=1;
for j=1:numel(Us)
 for k=1:2
 mode=1;if k==2,mode=0;end
 opt=odeset('RelTol',1e-9,'AbsTol',1e-11,'MaxStep',.05,'Events',@stop_event);f=ode15s(@(t,x)rhs(t,x,p,n,ps,p.U*Us(j),mode),[0 1],[g;1],opt);tf=linspace(0,f.x(end),501).';xf=deval(f,tf).';z=ode15s(@(t,x)rhs(t,x,p,n,ps,p.U,mode),[1 T],xf(end,:).',opt);tz=linspace(1,z.x(end),6001).';xz=deval(z,tz).';tt=[tf;tz];xx=[xf;xz];UU=[p.U*Us(j)*ones(size(tf));p.U*ones(size(tz))];d=electrical(xx(:,1),UU,p,n,mode);tail=tt>=.9*tt(end);recovered=all(abs(rad2deg(xx(tail,1)-g))<.1 & abs(xx(tail,2)-1)<1e-4);valid=all(abs(d.u)<=p.Vmax+1e-9);model="baseline";if k==2,model="limited";end
 row=table(model,Us(j),1,recovered,any(d.active),max(abs(d.u)),valid,max(abs(rad2deg(xx(:,1)-g))),xx(end,2),'VariableNames',{'model','Usag','Tfault','recoveredOriginal','limiterActivated','maxUrequired','voltageValidWholeTrajectory','maxAngleExcursion_deg','nuFinal'});summary=[summary;row];runs{j,k}=struct('time',tt,'state',xx,'voltage',UU,'electrical',d,'summary',row); %#ok<AGROW>
 end
 fig=figure('Color','w','Position',[30 30 1300 1150]);tl=tiledlayout(fig,4,2);labs={'gamma (deg)','nu','Unlimited demand (pu)','Delivered current (pu)','Limiter active','Required command voltage (pu)','Grid Pe (pu)'};
 for k=1:2,tr=runs{j,k};ys={rad2deg(tr.state(:,1)),tr.state(:,2),tr.electrical.Ireq,abs(tr.electrical.Iconv),double(tr.electrical.active),abs(tr.electrical.u),tr.electrical.Pe};for a=1:7,ax=nexttile(tl,a);hold(ax,'on');plot(ax,tr.time,ys{a},'DisplayName',char(tr.summary.model));grid(ax,'on');ylabel(ax,labs{a});xlabel(ax,'Normalized time');xlim(ax,[0 40]);end;end
 ax=nexttile(tl,1);legend(ax,'Location','best');ax=nexttile(tl,6);yline(ax,p.Vmax,'r--');title(tl,{'Network-consistent two-state ideal radial-current benchmark | SCR5/XR10',sprintf('Shared nominal Pstar=1; Usag=%.1f; Tfault=1; path voltage-valid=%d',Us(j),runs{j,2}.summary.voltageValidWholeTrajectory)});save_plot(fig,out,sprintf('J_nominal_sag_%02d',round(Us(j)*10)));
end
writetable(summary,fullfile(out,'nominal_sag_comparison.csv'));
end

function report(r,base,out)
p=r.parameters;n=r.network;fprintf('\nSCR5/XR10 LIMITER-FREE BASELINE\n');fprintf('E=%.12g; Vth=%.12g%+.12gj; Zeq=%.12g%+.12gj\nH=3.4 K=20 Kd=5 D_eff=25; nominal=%.12g deg\n',p.E,real(n.Vth),imag(n.Vth),real(n.Zeq),imag(n.Zeq),rad2deg(p.gammaNominal));fprintf('Local stability intervals:\n');disp(base.locallyStableIntervals_deg);fprintf('Baseline contraction-certified intervals:\n');disp(base.regions_deg);
fprintf('\nSCR5/XR10 + RADIAL CURRENT LIMITER\nImax=%.9g; voltage capability=%.9g\nSwitching angles:\n',p.Imax,p.Vmax);disp(r.switches_deg);fprintf('Voltage boundaries:\n');disp(r.voltageBoundaries_deg);fprintf('Local stable intervals (branchwise):\n');disp(r.intervals.stable);fprintf('Local unstable intervals:\n');disp(r.intervals.unstable);fprintf('Contraction-certified intervals, voltage valid, open smooth branches:\n');disp(r.regions);fprintf('NOMINAL:\n');disp(r.nominal);fprintf('EFFECT ON BASELINE REGIONS:\n');disp(r.coverage);
fprintf('Branchwise contraction negativity is NOT a global switched-system contraction certificate.\n');
fid=fopen(fullfile(out,'current_limiter_report.md'),'w','n','UTF-8');guard=onCleanup(@()fclose(fid));
fprintf(fid,'# SCR5/XR10 current-limiter comparison\n\n**Network-consistent two-state ideal radial-current benchmark.** Not an exact PI/freeze controller model. All electrical coefficients are read from the saved SCR5/XR10 baseline. H=3.4, K=20, Kd=5, D_eff=25; normalized analytical time and baseline metric method unchanged.\n\n');
fprintf(fid,'## Electrical closure\n\nVtarget=E exp(j(gamma+phi)); Vth=tU; q=(Vtarget-tU)/Zeq. Inactive: Iconv=q, Vc=Vtarget. Active: Iconv=1.2 q/abs(q), Vc=tU+Zeq Iconv. Grid voltage/current and transformer nodes use baseline nodal coefficients; Pe=real(Vg conj(Ig)). Gamma is the target angle during saturation, not actual capacitor angle. No rho-times-unlimited-power shortcut is used.\n\nZf is excluded from downstream Zeq and used only in Vbridge=Vc+Zf Iconv; u_required=Vbridge+0.2 Iconv. Points exceeding abs(u_required)=1.3 are marked voltage-invalid, with no added voltage saturation. The primary certificates require voltage validity.\n\n');
fprintf(fid,'## Derivatives and metric\n\nAll branches use analytic Pe_gamma, Pe_U and Pe_gammaU from voltage/current reconstruction. For r=abs(q), r_g=real(conj(q)q_g)/r, r_U=real(conj(q)q_U)/r, r_gU=real(conj(q_U)q_g)/r-r_g r_U/r. The active mixed current derivative is 1.2[-q_g r_U/r^2-q_U r_g/r^2-q r_gU/r^2+2q r_g r_U/r^3]. q_g=j Vtarget/Zeq and q_U=-t/Zeq.\n\nJ=[0,1;-Pe_gamma/(6.8nu),(-(Pstar-Pe)/nu^2-25)/6.8]. At each new equilibrium Pstar=Pe(gamma0), nu0=1. beta=-Pe_U/6.8; c=-Pe_gammaU/6.8; Theta=[-beta,0;-c,-beta]; M=Theta transpose times Theta; F=Theta J/Theta with fixed reference Theta_dot=0. No metric optimization.\n\nScan covers -180..180 deg at 0.001 deg spacing, partitioned at all refined current switches. Each interval endpoint is refined by a scalar root solve or is an exact switching/domain endpoint. Numerical certificate tolerances are 1e-9, matching baseline. Exact switches have only one-sided derivatives and no ordinary smooth certificate. Adjacent branches are not joined into a switched-system proof.\n\n');
fprintf(fid,'## Results\n\nNominal gamma=%.12g deg, Pe=%.12g, current=%.12g pu, current margin=%.12g pu, u_required=%.12g pu. Inactive=%d, voltage-valid=%d, local-stable=%d, contraction-certified=%d. Eigenvalues %.12g and %.12g; max eig SymF=%.12g. Nominal remains unchanged.\n\n',r.nominal.gamma0_deg,r.nominal.Pe,r.nominal.Iconv,p.Imax-r.nominal.Iconv,r.nominal.u_required,~r.nominal.limiter_active,r.nominal.voltage_capability_valid,r.nominal.local_stable,r.nominal.contraction_certified,r.nominal.eig1_real,r.nominal.eig2_real,r.nominal.lambda_max_SymF);
fprintf(fid,'Current-limit switching angles (deg): ');fprintf(fid,'%.12g ',r.switches_deg);fprintf(fid,'\n\nVoltage-capability boundaries (deg): ');if isempty(r.voltageBoundaries_deg),fprintf(fid,'none in scanned domain');else,fprintf(fid,'%.12g ',r.voltageBoundaries_deg);end;fprintf(fid,'\n\n');
fprintf(fid,'| Limited certified region | Lower deg | Upper deg | Mode | Lower cause | Upper cause |\n|---|---:|---:|---|---|---|\n');for k=1:height(r.regions),t=r.regions(k,:);fprintf(fid,'| %d | %.9f | %.9f | %s | %s | %s |\n',k,t.lower_gamma_deg,t.upper_gamma_deg,t.limiter_mode,t.lower_cause,t.upper_cause);end;if isempty(r.regions),fprintf(fid,'| None | | | | | |\n');end
fprintf(fid,'\nBaseline region coverage (widths in degrees; same target-angle domains, limited equilibria have new Pstar):\n\n| Region | Baseline interval | Current feasible | Current limited | Limited voltage valid | Limited locally stable | Limited certified | Lost certified width |\n|---|---|---:|---:|---:|---:|---:|---:|\n');for k=1:height(r.coverage),t=r.coverage(k,:);fprintf(fid,'| R%d | %.9f..%.9f | %.9f | %.9f | %.9f | %.9f | %.9f | %.9f |\n',k,t.lower_deg,t.upper_deg,t.currentFeasibleWidth_deg,t.currentLimitedWidth_deg,t.limitedVoltageValidWidth_deg,t.limitedLocalStableWidth_deg,t.limitedCertifiedWidth_deg,t.lostCertifiedWidth_deg);end
activeStable=0;activeCert=0;for k=1:size(r.intervals.active,1),activeStable=activeStable+overlap_width(r.intervals.active(k,:),r.intervals.stable);activeCert=activeCert+overlap_width(r.intervals.active(k,:),r.intervals.certified);end
fprintf(fid,'\nLimiter-active locally stable interval width=%.9g deg; limiter-active voltage-valid certified width=%.9g deg. Active does not automatically imply unstable. See all_classification_intervals.csv for exact intervals and current_limiter_regions.csv for boundary causes.\n\n',activeStable,activeCert);
for k=1:height(r.coverage),t=r.coverage(k,:);status='partially survives';if t.limitedCertifiedWidth_deg<1e-8,status='disappears';elseif abs(t.limitedCertifiedWidth_deg-t.width_deg)<1e-8,status='fully survives';end;fprintf(fid,'- Limiter-free R%d %s as a limited branchwise contraction-certified angle interval. Current-feasible here means the original prescribed-voltage equilibrium does not require current limitation.\n',k,status);end
fprintf(fid,'\nAll electrical network, E target, measurement mapping, tuning and metric methodology are identical; changes in power, equilibrium family, derivatives and certificates come from the ideal current projection and its necessary voltage reconstruction. Voltage-capability exclusion is separately reported, not voltage-saturation dynamics.\n\n');
fprintf(fid,'## Trajectories and optional sags\n\n+1 deg is our test choice; the narrow stable limiter-active import branch uses +0.1 deg to stay on that branch initially. The neighboring initial angle is always a further +0.01 deg. Baseline/limited trajectory pairs share target equilibrium angle but generally have DIFFERENT equilibrium Pstar; tables state both. For metric comparisons the SAME limited-reference Theta is used for both curves, with distance to their original common target-angle equilibrium. Invalid voltage paths remain mathematical diagnostics and are flagged.\n\nOptional nominal sag comparisons use identical Pstar=1 and initial equilibrium, retained voltages 0.5/0.3/0.1, duration 1, with zero angle jump. All requested electrical quantities are retained in MAT. Only nominal is used because old R1/R2 powers are not shared feasible equilibria. These representative sag comparisons are not a full duration/magnitude boundary sweep.\n\n');
fprintf(fid,'## Numerical verification\n\n');fn=fieldnames(r.validation);for k=1:numel(fn),fprintf(fid,'- %s: %.8g\n',fn{k},r.validation.(fn{k}));end
fprintf(fid,'\nAnalytical formulation is ready as an ideal benchmark for later EMT validation within its stated voltage-capability domain. It omits current/voltage PI states, delayed freeze and release behavior, measurement/Q-V/PLL dynamics and physical current overshoot. A smooth-branch test is not a global switched-system certificate. No EMT validation was started.\n');
end


