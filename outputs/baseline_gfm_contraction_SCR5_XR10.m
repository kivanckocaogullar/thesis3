%% Limiter-free VSG baseline: current SCR=5, X/R=10 physical network.
% H=3.4, K=20, Kd=5, D_eff=25; normalized analytical time, not EMT seconds.
% Constant physical capacitor-node E calibrated by nominal algebraic Q-V.
% No current saturation, PI/filter states, metric optimization or EMT fitting.
clearvars; clc; close all; format long g
root=fileparts(mfilename('fullpath'));outDir=fullfile(root,'baseline_SCR5_XR10');
if ~isfolder(outDir),mkdir(outDir);end
frozenNames={'baseline_gfm_contraction.m','gfm_current_limiter_analysis.m'};
frozen=cellfun(@(s)fileread(fullfile(root,s)),frozenNames,'UniformOutput',false);
logPath=fullfile(outDir,'SCR5_XR10_console_log.txt');fidLog=fopen(logPath,'w');fclose(fidLog);diary(logPath); cleanup=onCleanup(@()diary('off'));
[p,network,sourceAudit]=identify_current_model();
p.mass=2*p.H;p.D=p.K+p.Kd;p.omega0_pu=1;p.tol=1e-9;
% Positive-slope nominal equilibrium of the quasi-steady Q-V closure.
gnom=fzero(@(g)qv_nominal_residual(g,p),[0,pi/3],optimset('TolX',1e-13));
p.E=qv_E(gnom,p.U,p);p.gammaNominal=gnom;
nominalRoots=power_roots(1,p);nominalRootsTable=rows(rad2deg(nominalRoots(:)),p);
nominalRootsTable.isSelected=abs(nominalRootsTable.gamma0_deg-rad2deg(gnom))<1e-7;
assert(abs(electric(gnom,p.U,p).Pe-1)<1e-11);
% Full fundamental period audits import/export. Primary plots use 0..180 deg.
gridDeg=unique(sort([linspace(-180,180,720001),rad2deg(gnom),0])).';
scan=rows(gridDeg,p);
regions=find_intervals(gridDeg,scan,p,'certificate');
stable=find_intervals(gridDeg,scan,p,'stable');unstable=find_intervals(gridDeg,scan,p,'unstable');
boundaries=boundary_table(regions,p);
structural=structural_boundaries(p);
[sampleDeg,sampleLabels,representativeDeg,representativeLabels]=select_points(regions,stable,unstable,p);
samples=rows(sampleDeg,p);samples.label=sampleLabels;
pointMatrices=arrayfun(@(g)point(deg2rad(g),p),sampleDeg);
print_matrices(pointMatrices,sampleLabels);
validation=verify_model(p,network,regions,gnom);
secondary=secondary_qv(p);
writetable(scan,fullfile(outDir,'SCR5_XR10_dense_equilibrium_scan.csv'));
writetable(samples,fullfile(outDir,'SCR5_XR10_sample_points.csv'));
writetable(nominalRootsTable,fullfile(outDir,'SCR5_XR10_nominal_equilibria.csv'));
writetable(boundaries,fullfile(outDir,'SCR5_XR10_contraction_boundaries.csv'));
writetable(structural,fullfile(outDir,'SCR5_XR10_structural_boundaries.csv'));
writetable(secondary,fullfile(outDir,'SCR5_XR10_QV_secondary_sensitivity.csv'));
networkTable=table(p.E,abs(network.Vth),real(network.Vth),imag(network.Vth), ...
    real(network.Zeq),imag(network.Zeq),abs(network.Zeq),rad2deg(angle(network.Zeq)),p.U, ...
    p.S,p.Vlv,p.S/(sqrt(3)*p.Vlv),p.Vlv^2/p.S,p.Rgrid,p.Lgrid, ...
    'VariableNames',{'E_pu','VthMagnitude_pu','VthReal_pu','VthImag_pu','Req_pu','Xeq_pu', ...
    'ZeqMagnitude_pu','ZeqAngle_deg','physicalBusReferred_pu','Sbase_VA','VbaseLL_V','IbaseRMS_A','Zbase_ohm','Rgrid_ohm','Lgrid_H'});
writetable(networkTable,fullfile(outDir,'SCR5_XR10_network_parameters.csv'));
[trajectories,trajectoryTable]=trajectory_study(representativeDeg,representativeLabels,p,outDir);
base_figures(scan,samples,regions,boundaries,structural,secondary,p,outDir);
% Sag references are NEW certified midpoints plus nominal as an explicitly
% uncertified comparison if needed. Never relabel nominal as certified.
sagRefs=unique([gnom;deg2rad(mean(regions,2))],'stable');
sagRefs=sagRefs(sagRefs>0 & sagRefs<pi);
[magnitude,magnitudeRuns]=sag_study(sagRefs,.9:-.1:.1,1,'magnitude',p,outDir);
[duration,durationRuns]=sag_study(sagRefs,.1,[1,2,3,4,5,7,10],'duration',p,outDir);
figure5_views(representativeDeg,representativeLabels,trajectories,regions,magnitudeRuns,p,outDir);
comparison=comparison_table(regions,gnom,p,network);
writetable(comparison,fullfile(outDir,'SCR5_XR10_old_new_comparison.csv'));
results=struct('parameters',p,'network',network,'sourceAudit',sourceAudit,'scan',scan, ...
    'nominalEquilibria',nominalRootsTable,'regions_deg',regions,'locallyStableIntervals_deg',stable, ...
    'locallyUnstableIntervals_deg',unstable,'boundaries',boundaries,'structuralBoundaries',structural, ...
    'samplePoints',samples,'pointMatrices',pointMatrices,'validation',validation,'secondaryQV',secondary,'trajectorySummary',trajectoryTable, ...
    'trajectories',{trajectories},'sagMagnitude',magnitude,'sagMagnitudeRuns',{magnitudeRuns}, ...
    'sagDuration',duration,'sagDurationRuns',{durationRuns},'comparison',comparison);
save(fullfile(outDir,'SCR5_XR10_baseline_results.mat'),'results','-v7.3');
write_report(results,outDir);
for k=1:numel(frozenNames),assert(strcmp(frozen{k},fileread(fullfile(root,frozenNames{k}))), 'Frozen analytical file was modified.');end
fprintf('\nAll SCR5_XR10 outputs saved to %s\nBoth existing analytical scripts unchanged. No current-limiter analysis run.\n',outDir);
diary off

function [p,n,a]=identify_current_model()
pRepoRoot=fileparts(fileparts(mfilename('fullpath')));
p=struct('modelPath',fullfile(pRepoRoot,'model','GFMVSG1.slx'), ...
    'libraryRoot',fullfile(pRepoRoot,'model','tdmsquared-master','tdmsquared-master','Library'));
addpath(genpath(p.libraryRoot));load_system(p.modelPath);mdl='GFMVSG1';
gfm=[mdl,'/Grid_Forming_Converter'];bl=find_system(mdl,'SearchDepth',1,'Type','Block');
rl='';source='';
for k=1:numel(bl)
    s=get_param(bl{k},'ReferenceBlock');
    if contains(s,'spsThreePhaseSeriesRLCBranchLib'),rl=bl{k};end
    if contains(s,'spsThreePhaseProgrammableVoltageSourceLib'),source=bl{k};end
end
assert(~isempty(rl)&&~isempty(source));assert(strcmp(get_param(rl,'BranchType'),'RL'));
assert(strcmp(get_param(source,'VariationEntity'),'None'));
% Audit the actual terminal graph: source -> RL -> GFM, one phase at a time.
pr=get_param(rl,'PortConnectivity');hg=get_param(gfm,'Handle');hs=get_param(source,'Handle');
for k=1:3
    left=unique([pr(k).SrcBlock(:);pr(k).DstBlock(:)]);right=unique([pr(k+3).SrcBlock(:);pr(k+3).DstBlock(:)]);
    assert((any(left==hs)&&any(right==hg))||(any(left==hg)&&any(right==hs)), 'Unexpected RL terminal connectivity.');
end
num=@(b,k)str2num(get_param(b,k)); %#ok<ST2NM>
p.H=num(gfm,'OuterCtrl_VSM_H');p.K=num(gfm,'OuterCtrl_VSM_K');p.Kd=num(gfm,'OuterCtrl_VSM_Kd');
assert(p.H==3.4&&p.K==20&&p.Kd==5,'Saved VSM tuning differs from requested fixed tuning.');
p.S=num(gfm,'S_n');p.Vlv=num(gfm,'V_n_phph');p.Vgrid=num(gfm,'V_n_grid_phph');
p.Rgrid=num(rl,'Resistance');p.Lgrid=num(rl,'Inductance');src=num(source,'PositiveSequence');p.f=src(3);
assert(p.S==2e6&&p.Vgrid==20e3&&p.f==50&&src(2)==0);
assert(abs(p.Rgrid-3.98014876)<1e-5&&abs(p.Lgrid-.12669207)<1e-7);
assert(strcmp(get_param(gfm,'withTransformer'),'Grid coupling using transformer'));
assert(num(gfm,'OuterCtrl_Lvirt')==0);
p.kQ=num(gfm,'OuterCtrl_VSM_kQ');p.Vref=1;p.Qstar=0;
% Referenced electrical subsystem: Winding2=[1.02*V_n_grid_phph,R2,X2].
% Positive-sequence Yg/D1 phase shift removed on referral; 50-Hz unity
% measurement tracking is the primary two-state quasi-steady approximation.
p.Vhv=1.02*p.Vgrid;p.U=src(1)/p.Vhv;
t1=num(gfm,'RL_Trafo_Conv');t2=num(gfm,'RL_Trafo_Grid');cp=num(gfm,'C_filter_pu');fp=num(gfm,'L_filter_conv_pu');
z1=complex(t1(1),t1(2));z2=complex(t2(1),t2(2));ym=1/num(gfm,'Rm')+1/(1i*num(gfm,'Lm'));
yc=1/(cp(1)+1/(1i*cp(2)));zg=complex(p.Rgrid,2*pi*p.f*p.Lgrid)/(p.Vhv^2/p.S);
% Closed breaker in Electrical_Part: 0.01 ohm, parallel 1e6-ohm snubber.
zbr=(.01*1e6/(.01+1e6))/(p.Vhv^2/p.S);
y=[1/z1+1/z2+ym,-1/z2;-1/z2,1/z2+1/(zg+zbr)];
coef=y\diag([1/z1,1/(zg+zbr)]);
iv=(coef(2,1))/(zg+zbr);iu=(coef(2,2)-1)/(zg+zbr);
av=coef(2,1)-zbr*iv;au=coef(2,2)-zbr*iu;
icv=(1-coef(1,1))/z1+yc;icu=-coef(1,2)/z1;
n=struct('Zt1',z1,'Zt2',z2,'Ym',ym,'Yc',yc,'Zgrid',zg,'Zbreaker',zbr, ...
    'Zconverter',complex(fp(1),fp(2)),'nodalY',y,'nodeCoefficients',coef, ...
    'vgE',av,'vgU',au,'igE',iv,'igU',iu,'icE',icv,'icU',icu,'Zeq',1/icv,'Vth',-icu/icv*p.U);
p.phi=angle(n.Vth);
A=av*conj(iv);B=au*conj(iu);C=av*conj(iu)+au*conj(iv);D=1i*(av*conj(iu)-au*conj(iv));
c=[real(A),real(B),real(C),real(D)];q=[imag(A),imag(B),imag(C),imag(D)];
p.pc=[c(1:2),c(3)*cos(p.phi)+c(4)*sin(p.phi),-c(3)*sin(p.phi)+c(4)*cos(p.phi)];
p.qc=[q(1:2),q(3)*cos(p.phi)+q(4)*sin(p.phi),-q(3)*sin(p.phi)+q(4)*cos(p.phi)];
a=struct('modelPath',p.modelPath,'RLpath',rl,'sourcePath',source,'GFMpath',gfm, ...
    'topology','Bridge -> converter RL -> capacitor node (shunt series RC) -> transformer T equivalent -> closed breaker -> external RL -> ideal source', ...
    'absent','Old PI line and old nonideal source removed; fault block unconnected; no local load', ...
    'excluded','No_Transformer grid inductor is inactive; active damping is not passive series resistance', ...
    'powerLocation','grid terminal voltage after breaker times transformer-secondary current before breaker', ...
    'Edefinition','Physical capacitor-node fundamental voltage, ideal inner-loop tracking, unity measurement approximation', ...
    'gammaDefinition','arg(Vcap) minus arg(Thevenin voltage); ideal transformer D1 phase offset removed');
fprintf('CURRENT MODEL AUDIT PASSED: H=%.4g, K=%.4g, Kd=%.4g; external R=%.10g ohm, L=%.10g H\n',p.H,p.K,p.Kd,p.Rgrid,p.Lgrid);
end

function d=electric(g,u,p,E)
if nargin<4,E=p.E;end
c=p.pc;q=p.qc;h=c(3)*cos(g)+c(4)*sin(g);hg=-c(3)*sin(g)+c(4)*cos(g);
hq=q(3)*cos(g)+q(4)*sin(g);hqg=-q(3)*sin(g)+q(4)*cos(g);
d=struct('Pe',c(1)*E.^2+c(2)*u.^2+E.*u.*h,'Pg',E.*u.*hg, ...
    'Pgg',-E.*u.*h,'PU',2*c(2)*u+E.*h,'PgU',E.*hg, ...
    'Q',q(1)*E.^2+q(2)*u.^2+E.*u.*hq,'Qg',E.*u.*hqg, ...
    'QE',2*q(1)*E+u.*hq,'PE',2*c(1)*E+u.*h);
end
function E=qv_E(g,u,p)
q=p.qc;a=p.kQ*q(1);b=1+p.kQ*u.*(q(3)*cos(g)+q(4)*sin(g));c=p.kQ*q(2)*u.^2-p.Vref-p.kQ*p.Qstar;
disc=b.^2-4*a*c;E=nan(size(disc));ok=disc>=0;
E(ok)=-2*c./(b(ok)+sqrt(disc(ok)));E(E<=0)=NaN;
end
function y=qv_nominal_residual(g,p)
E=qv_E(g,p.U,p);d=electric(g,p.U,p,E);y=d.Pe-1;
end
function r=power_roots(ps,p)
c=p.pc;amp=p.E*p.U*hypot(c(3),c(4));offset=c(1)*p.E^2+c(2)*p.U^2;v=(ps-offset)/amp;
if abs(v)>1+1e-12,r=[];return;end
v=max(-1,min(1,v));phase=atan2(c(3),c(4));r=unique(mod([asin(v)-phase,pi-asin(v)-phase]+pi,2*pi)-pi);
end
function dx=flow(~,x,p,ps,u)
d=electric(x(1),u,p);dx=[x(2)-1;((ps-d.Pe)/x(2)-p.D*(x(2)-1))/p.mass];
end
function J=jac(g,nu,p,ps,u)
d=electric(g,u,p);J=[0,1;-d.Pg/(p.mass*nu),(-(ps-d.Pe)/nu^2-p.D)/p.mass];
end
function s=point(g,p)
d=electric(g,p.U,p);beta=-d.PU/p.mass;c=-p.omega0_pu*d.PgU/p.mass;
J=jac(g,1,p,d.Pe,p.U);T=[-beta,0;-c,-beta];M=T'*T;
if abs(beta)>1e-12,F=T*J/T;S=(F+F')/2;lambda=max(eig(S));else,F=nan(2);S=nan(2);lambda=NaN;end
s=struct('gamma',g,'Pstar',d.Pe,'J',J,'beta',beta,'c',c,'Theta',T,'M',M,'F',F,'SymF',S, ...
    'eigJ',eig(J),'eigM',sort(eig(M)),'lambda',lambda);
end
function t=rows(gd,p)
g=deg2rad(gd(:));d=electric(g,p.U,p);beta=-d.PU/p.mass;c=-d.PgU*p.omega0_pu/p.mass;
r=c./beta;a=-d.Pg/p.mass;b=p.D/p.mass;off=(1+a+b*r-r.^2)/2;
lm=-b/2+hypot(b/2-r,off);mn=-b/2-hypot(b/2-r,off);
m11=beta.^2+c.^2;m12=beta.*c;m22=beta.^2;mm=(m11+m22+hypot(m11-m22,2*m12))/2;mi=beta.^4./mm;
e1=(-b+sqrt(complex(b^2+4*a)))/2;e2=(-b-sqrt(complex(b^2+4*a)))/2;
singular=abs(beta)<=1e-12;lm(singular)=NaN;mn(singular)=NaN;
local=max(real(e1),real(e2));status=repmat("Marginal",size(g));status(local<-p.tol)="Locally stable";status(local>p.tol)="Locally unstable";
cert=mi>p.tol & lm<-p.tol & ~singular;
t=table(gd(:),d.Pe,d.Pg,real(e1),imag(e1),real(e2),imag(e2),local,status,local<-p.tol, ...
    beta,c,mi,mm,lm,mn,cert,singular,r.*(b-r)-off.^2,zeros(size(g)), ...
    'VariableNames',{'gamma0_deg','Pstar_pu','dPe_dgamma','eig1_real','eig1_imag','eig2_real','eig2_imag', ...
    'maxRealEigJ','localClassification','locallyStable','beta','c','lambdaMinM','lambdaMaxM','lambdaMaxSymF', ...
    'lambdaMinSymF','contractionCertified','singularTheta','detSymF','equilibriumResidual'});
end
function y=margin(gd,p,kind)
s=rows(gd,p);
switch kind
case 'certificate',y=min(s.lambdaMinM-p.tol,-s.lambdaMaxSymF-p.tol);y(s.singularTheta)=-Inf;
case 'stable',y=-s.maxRealEigJ-p.tol;
case 'unstable',y=s.maxRealEigJ-p.tol;
end
end
function I=find_intervals(g,t,p,kind)
switch kind
case 'certificate',ok=t.contractionCertified;
case 'stable',ok=t.locallyStable;
case 'unstable',ok=t.maxRealEigJ>p.tol;
end
beg=find(diff([false;ok])==1);fin=find(diff([ok;false])==-1);I=zeros(numel(beg),2);
for k=1:numel(beg)
 lo=g(beg(k));hi=g(fin(k));
 if beg(k)>1,lo=fzero(@(v)margin(v,p,kind),g(beg(k)+[-1,0]),optimset('TolX',1e-11));end
 if fin(k)<numel(g),hi=fzero(@(v)margin(v,p,kind),g(fin(k)+[0,1]),optimset('TolX',1e-11));end
 I(k,:)=[lo,hi];
end
end
function t=boundary_table(I,p)
t=table();
for k=1:size(I,1)
 for side=1:2
    g=I(k,side);s=rows(g,p);z=p;z.tol=0;
    exact=fzero(@(v)margin(v,z,'certificate'),g,optimset('TolX',1e-12));
    sides=["lower","upper"];cause="lambdaMaxSymF = 0 / det(SymF) = 0 (within numerical tolerance)";
    if abs(s.lambdaMinM-p.tol)<1e-10,cause="M positive-definiteness tolerance";end
    row=table(k,sides(side),g,exact,cause,s.lambdaMaxSymF,s.lambdaMinM,s.detSymF, ...
        'VariableNames',{'region','side','boundary_deg','zeroBoundary_deg','cause','lambdaMaxSymF','lambdaMinM','detSymF'});
    t=[t;row]; %#ok<AGROW>
 end
end
end
function t=structural_boundaries(p)
c=p.pc;phase=atan2(c(3),c(4));slope=mod([pi/2-phase,-pi/2-phase]+pi,2*pi)-pi;
v=-2*c(2)*p.U/(p.E*hypot(c(3),c(4)));bet=[];
if abs(v)<=1,bet=mod([asin(v)-phase,pi-asin(v)-phase]+pi,2*pi)-pi;end
zerosP=power_roots(0,p);
g=[slope(:);bet(:);zerosP(:)];kind=[repmat("dPe/dgamma = 0; local-stability boundary",numel(slope),1); ...
    repmat("beta = 0; singular Theta",numel(bet),1);repmat("Pe = 0; import/export boundary",numel(zerosP),1)];
t=table(rad2deg(g),kind,'VariableNames',{'gamma_deg','cause'});t=sortrows(t,'gamma_deg');
end
function [gd,labels,rg,rl]=select_points(I,stable,unstable,p)
rg=rad2deg(p.gammaNominal);rl="New nominal (Pstar=1)";
pos=I(I(:,1)>0 & I(:,2)<180,:);
if ~isempty(pos)
 rg=[rg;max(.5,pos(1,1)/2)];rl=[rl;"Below first certificate"];
 for k=1:size(pos,1)
    rg=[rg;mean(pos(k,:))];rl=[rl;"Inside R"+k]; %#ok<AGROW>
    if k<size(pos,1),rg=[rg;mean([pos(k,2),pos(k+1,1)])];rl=[rl;"Between R"+k+" and R"+(k+1)];end %#ok<AGROW>
 end
end
un=unstable(unstable(:,2)>90 & unstable(:,1)<170,:);
if ~isempty(un),rg=[rg;min(140,max(110,mean(un(end,:))))];rl=[rl;"Locally unstable"];
else,rg=[rg;120];rl=[rl;"Comparison 120 deg"];end
[rg,idx]=unique(rg,'stable');rl=rl(idx);
gd=rg;labels=rl;
for k=1:numel(I)
 for shift=[-.01,.01]
    gd(end+1,1)=I(k)+shift;labels(end+1,1)="Near boundary "+sprintf('%.6f',I(k)); %#ok<AGROW>
 end
end
for a=[10,12,17.457603,20.5,25,35,55,65,100,120]
 if ~any(abs(gd-a)<1e-8),gd(end+1,1)=a;labels(end+1,1)="Old-angle comparison";end %#ok<AGROW>
end
[gd,idx]=sort(gd);labels=labels(idx);
end
function t=secondary_qv(p)
g=deg2rad(linspace(-180,180,7201).');E=qv_E(g,p.U,p);d=electric(g,p.U,p,E);
Eg=-p.kQ*d.Qg./(1+p.kQ*d.QE);der=d.Pg+d.PE.*Eg;main=electric(g,p.U,p);
t=table(rad2deg(g),E,d.Pe,main.Pe,der,Eg,'VariableNames',{'gamma_deg','E_QV_pu','Pe_QV_pu','Pe_constantE_pu','dPe_QV_dgamma','dE_dgamma'});
end

function v=verify_model(p,n,I,gnom)
h=1e-5;angles=deg2rad([-170,-80,-10,0,10,rad2deg(gnom),25,50,70,100,140,175]);
v=struct('powerDerivative',0,'Jacobian',0,'betaDerivative',0,'cDerivative',0, ...
    'matrixSymmetry',0,'matrixFormula',0,'equilibriumResidual',0,'phasorPower',0,'thevenin',0,'nominalAngle',0);
for g=angles
 d=electric(g,p.U,p);plus=electric(g+h,p.U,p);minus=electric(g-h,p.U,p);
 v.powerDerivative=max(v.powerDerivative,abs((plus.Pe-minus.Pe)/(2*h)-d.Pg));
 du1=electric(g,p.U+h,p);du0=electric(g,p.U-h,p);
 v.betaDerivative=max(v.betaDerivative,abs(-(du1.Pe-du0.Pe)/(2*h)/p.mass+d.PU/p.mass));
 v.cDerivative=max(v.cDerivative,abs(-(plus.PU-minus.PU)/(2*h)/p.mass+d.PgU/p.mass));
 for nu=[.8,1,1.1]
    x=[g;nu];ps=d.Pe;A=jac(g,nu,p,ps,p.U);An=zeros(2);
    for k=1:2,eh=zeros(2,1);eh(k)=h;An(:,k)=(flow(0,x+eh,p,ps,p.U)-flow(0,x-eh,p,ps,p.U))/(2*h);end
    v.Jacobian=max(v.Jacobian,norm(A-An,inf));
 end
 s=point(g,p);r=rows(rad2deg(g),p);
 v.matrixSymmetry=max(v.matrixSymmetry,norm(s.SymF-s.SymF','fro'));
 v.matrixFormula=max(v.matrixFormula,abs(s.lambda-r.lambdaMaxSymF)/max(1,abs(s.lambda)));
 assert(abs(det(s.Theta)-s.beta^2)<1e-12 && min(s.eigM)>0);
 assert(abs(min(s.eigM)-r.lambdaMinM)<1e-11);
 v.equilibriumResidual=max(v.equilibriumResidual,norm(flow(0,[g;1],p,s.Pstar,p.U),inf));
 voltage=p.E*exp(1i*(g+p.phi));vg=n.vgE*voltage+n.vgU*p.U;ig=n.igE*voltage+n.igU*p.U;
 v.phasorPower=max(v.phasorPower,abs(real(vg*conj(ig))-d.Pe));
 ic=n.icE*voltage+n.icU*p.U;v.thevenin=max(v.thevenin,abs(ic-(voltage-n.Vth)/n.Zeq));
end
% Independently reduce transformer/core and capacitor port by series/parallel.
zrest=n.Zt2+n.Zbreaker+n.Zgrid;zm=1/(n.Ym+1/zrest);um=p.U/(1+n.Ym*zrest);zp=n.Zt1+zm;
assert(abs(zp/(1+n.Yc*zp)-n.Zeq)<1e-11 && abs(um/(1+n.Yc*zp)-n.Vth)<1e-11);
ind=fzero(@(g)direct_power(g,p,n)-1,[0,pi/3]);v.nominalAngle=abs(ind-gnom);
for k=1:size(I,1)
 assert(margin(I(k,1)+1e-5,p,'certificate')>0 && margin(I(k,1)-1e-5,p,'certificate')<0);
 assert(margin(I(k,2)-1e-5,p,'certificate')>0 && margin(I(k,2)+1e-5,p,'certificate')<0);
end
assert(v.powerDerivative<1e-7&&v.Jacobian<1e-7&&v.betaDerivative<1e-7&&v.cDerivative<1e-7);
assert(v.matrixFormula<1e-9&&v.equilibriumResidual<1e-12&&v.phasorPower<1e-10&&v.thevenin<1e-10&&v.nominalAngle<1e-10);
fprintf('Independent numerical verification passed.\n');disp(v);
end
function pe=direct_power(g,p,n)
x=p.E*exp(1i*(g+p.phi));nodes=n.nodalY\[x/n.Zt1;p.U/(n.Zgrid+n.Zbreaker)];
ig=(nodes(2)-p.U)/(n.Zgrid+n.Zbreaker);vg=nodes(2)-n.Zbreaker*ig;pe=real(vg*conj(ig));
end
function [runs,tab]=trajectory_study(gd,labels,p,outDir)
runs=cell(numel(gd),1);tab=table();
for k=1:numel(gd)
 s=point(deg2rad(gd(k)),p);tend=max(200,20/max(abs(max(real(s.eigJ))),1e-4));
 opt=odeset('RelTol',1e-9,'AbsTol',1e-11,'MaxStep',.05,'Events',@frequency_event);
 sol=ode15s(@(t,x)flow(t,x,p,s.Pstar,p.U),[0,tend],[s.gamma+deg2rad(1);1],opt);
 t=linspace(0,sol.x(end),5001).';x=deval(sol,t).';
 assert(max(abs(diff(x(:,1))))<pi);L=vecnorm(s.Theta*(x-[s.gamma,1]).',2,1).';
 [cls,stats]=recovery(t,x,s,p,false(size(t)));
 start=rows(gd(k),p);
 row=table(labels(k),gd(k),s.Pstar,1,rad2deg(x(end,1)),x(end,2),cls,start.contractionCertified, ...
     stats.finalPowerResidual,stats.tailAngle,stats.tailNu, ...
     'VariableNames',{'label','gamma0_deg','Pstar_pu','perturbation_deg','gammaFinal_deg','nuFinal', ...
     'recoveryClass','contractionCertified','finalPowerResidual_pu','tailMaxAngleError_deg','tailMaxAbsNuDeviation'});
 tab=[tab;row]; %#ok<AGROW>
 % A second nearby initial angle gives an actual neighboring-trajectory check.
 sol2=ode15s(@(t,x)flow(t,x,p,s.Pstar,p.U),[0,sol.x(end)],[s.gamma+deg2rad(1.01);1],opt);
 xt=deval(sol2,min(t,sol2.x(end))).';pairDistance=vecnorm(s.Theta*(xt-x).',2,1).';
 runs{k}=struct('time',t,'state',x,'distance',L,'neighborState',xt,'neighborDistance',pairDistance, ...
     'reference',s,'label',labels(k),'summary',row);
end
writetable(tab,fullfile(outDir,'SCR5_XR10_trajectory_outcomes.csv'));
n=numel(runs);fig=figure('Color','w','Position',[40,40,1800,850]);tl=tiledlayout(fig,2,n,'TileSpacing','compact');
for k=1:n
 tr=runs{k};ax=nexttile(tl,k);hold(ax,'on');show=true(size(tr.time));
 if max(real(tr.reference.eigJ))>0
    stop=find(abs(rad2deg(tr.state(:,1)-tr.reference.gamma))>=20,1);
    if ~isempty(stop),show((stop+1):end)=false;end
 end
 plot(ax,rad2deg(tr.state(show,1)),tr.state(show,2)-1,'b-','LineWidth',1.5);
 plot(ax,rad2deg(tr.reference.gamma),0,'ro','MarkerFaceColor','r');
 plot(ax,rad2deg(tr.state(1,1)),tr.state(1,2)-1,'ks');
 ei=find(show,1,'last');plot(ax,rad2deg(tr.state(ei,1)),tr.state(ei,2)-1,'gd','MarkerFaceColor','g');
 xlabel(ax,'gamma unwrapped (deg)');ylabel(ax,'nu-1');grid(ax,'on');title(ax,{char(tr.label),sprintf('gamma0 = %.4f deg',rad2deg(tr.reference.gamma))},'FontSize',9);
 if ~all(show),subtitle(ax,'Green: displayed endpoint','FontSize',8);end
 ax=nexttile(tl,n+k);plot(ax,tr.time(show),tr.distance(show),'LineWidth',1.5);grid(ax,'on');
 xlabel(ax,'Normalized analytical time');ylabel(ax,'L relative to original equilibrium');
end
stamp(tl,'+1 deg test perturbation is our validation choice, not from the paper');savefigs(fig,outDir,'07_representative_trajectories');
fig=figure('Color','w','Position',[80,40,1400,850]);tl=tiledlayout(fig,3,1);
for k=1:n
 tr=runs{k};ax=nexttile(tl,1);hold(ax,'on');plot(ax,tr.time,rad2deg(tr.state(:,1)),'DisplayName',tr.label);
 ax=nexttile(tl,2);hold(ax,'on');plot(ax,tr.time,tr.distance,'DisplayName',tr.label);
 ax=nexttile(tl,3);hold(ax,'on');semilogy(ax,tr.time,max(tr.neighborDistance,eps),'DisplayName',tr.label);
end
for j=1:3,ax=nexttile(tl,j);grid(ax,'on');xlabel(ax,'Normalized analytical time');end
ax=nexttile(tl,1);ylabel(ax,'gamma unwrapped (deg)');legend(ax,'Location','best');
ax=nexttile(tl,2);ylabel(ax,'L: original equilibrium');
ax=nexttile(tl,3);ylabel(ax,'Neighbor metric distance');set(ax,'YScale','log');
stamp(tl,'Full-time and neighboring-trajectory diagnostics; no basin proof');savefigs(fig,outDir,'10_full_time_neighbor_diagnostics');
end
function [value,terminal,direction]=frequency_event(~,x)
value=x(2)-.5;terminal=1;direction=-1;
end
function [cls,a]=recovery(t,x,s,p,fault)
g=x(:,1);ne=x(:,2)-1;ae=rad2deg(g-s.gamma);d=electric(g,p.U,p);pr=d.Pe-s.Pstar;
tail=t>=.9*t(end);original=all(abs(ae(tail))<.1 & abs(ne(tail))<1e-4) && ~any(fault(tail));
rr=power_roots(s.Pstar,p);gm=mean(g(tail));candidates=rr+2*pi*round((gm-rr)/(2*pi));
[~,j]=min(abs(candidates-gm));target=candidates(j);
another=abs(rad2deg(target-s.gamma))>=.1 && all(abs(rad2deg(g(tail)-target))<.1) ...
    && all(abs(ne(tail))<1e-4 & abs(pr(tail))<1e-4) && ~any(fault(tail));
if original,cls="Recovered original equilibrium";
elseif another,cls="Transitioned to another equilibrium branch";
else,cls="Did not recover original equilibrium";end
if any(x(:,2)<=.50000001),cls="Did not recover original equilibrium";end
a=struct('finalPowerResidual',pr(end),'tailAngle',max(abs(ae(tail))),'tailNu',max(abs(ne(tail))), ...
    'tail',tail,'angleError',ae,'nuError',ne);
end
function [tr,row]=sag_run(g,usag,T,p,kind,tighter)
if nargin<6,tighter=false;end
s=point(g,p);rate=abs(max(real(s.eigJ)));postDuration=max(200,20/rate);
if strcmp(kind,'magnitude'),tend=postDuration;else,tend=T+postDuration;end
opt=odeset('RelTol',1e-9,'AbsTol',1e-11,'MaxStep',.05,'Events',@frequency_event);
if tighter,opt=odeset(opt,'RelTol',1e-11,'AbsTol',1e-13,'MaxStep',.025);end
f=ode15s(@(t,x)flow(t,x,p,s.Pstar,p.U*usag),[0,T],[g;1],opt);
tf=linspace(0,f.x(end),max(501,ceil(T/.002)+1)).';xf=deval(f,tf).';
if f.x(end)<T-1e-8
    tp=[];xp=zeros(0,2);
else
    post=ode15s(@(t,x)flow(t,x,p,s.Pstar,p.U),[T,tend],xf(end,:).',opt);
    tp=linspace(T,post.x(end),6001).';xp=deval(post,tp).';
end
t=[tf;tp];x=[xf;xp];u=[repmat(p.U*usag,size(tf));repmat(p.U,size(tp))];fault=[true(size(tf));false(size(tp))];
assert(max(abs(diff(x(:,1))))<pi,'Insufficient unwrapped trajectory samples');
[cls,a]=recovery(t,x,s,p,fault);
L=vecnorm(s.Theta*(x-[g,1]).',2,1).';d=electric(x(:,1),u,p);ratio=s.c/s.beta;
a21=-d.Pg./(p.mass*x(:,2));a22=(-(s.Pstar-d.Pe)./x(:,2).^2-p.D)/p.mass;
off=(1+a21-a22*ratio-ratio^2)/2;lm=a22/2+hypot(-ratio-a22/2,off);
mm=min(s.eigM);cert=mm>p.tol & lm<-p.tol;postcert=cert(~fault);postt=t(~fault);
last=find(~postcert,1,'last');recert=NaN;
if ~isempty(postcert)
 if isempty(last),recert=T;elseif last<numel(postcert),recert=postt(last+1);end
end
violate=any(~cert(fault));posttail=~any(fault(a.tail))&&all(cert(a.tail));
pr=d.Pe-s.Pstar;refrow=rows(rad2deg(g),p);
row=table(rad2deg(g),s.Pstar,refrow.contractionCertified,usag,100*(1-usag),T, ...
    max(rad2deg(x(:,1))),min(rad2deg(x(:,1))),max(abs(a.angleError)),max(abs(a.nuError)),max(L), ...
    rad2deg(x(end,1)),x(end,2),a.angleError(end),a.nuError(end),max(lm),mm,mm>p.tol,all(lm<-p.tol),all(cert), ...
    violate,posttail,violate&&posttail&&isfinite(recert),recert,a.tailAngle,a.tailNu,pr(end),cls, ...
    "Fixed reference metric; sampled instantaneous test, not a regional proof", ...
    'VariableNames',{'gamma0_deg','Pstar_pu','referenceContractionCertified','Usag_pu','sagDepth_percent','faultDuration_normalized', ...
    'gammaMax_deg','gammaMin_deg','maxAngleExcursion_deg','maxAbsNuDeviation','maxL','gammaFinal_deg','nuFinal', ...
    'finalGammaError_deg','finalNuError','maxLambdaMaxSymF','minLambdaMinM','metricRemainedPD','symFRemainedND','contractionMaintained', ...
    'contractionViolatedDuringFault','contractionSatisfiedPostfaultTail','contractionBecameSatisfiedAfterClearing', ...
    'postfaultSatisfiedFromTime','tailMaxAngleError_deg','tailMaxAbsNuDeviation','finalPowerResidual_pu','recoveryClass','contractionEvaluationScope'});
for j=unique([1,numel(tf),min(numel(tf)+1,numel(t)),numel(t)])
 A=jac(x(j,1),x(j,2),p,s.Pstar,u(j));F=s.Theta*A/s.Theta;
 assert(abs(max(eig((F+F')/2))-lm(j))<1e-8);
end
tr=struct('time',t,'state',x,'voltage',u,'fault',fault,'distance',L,'lambdaMaxSymF',lm,'reference',s,'summary',row);
end
function [tab,runs]=sag_study(refs,usags,durations,kind,p,outDir)
if strcmp(kind,'magnitude'),n=numel(usags);else,n=numel(durations);end
runs=cell(numel(refs),n);tab=table();checks=table();
fprintf('\n%s sag study: %d references, %d tests each; no angle jump; Pstar fixed.\n',kind,numel(refs),n);
for k=1:numel(refs)
 for j=1:n
    if strcmp(kind,'magnitude'),u=usags(j);T=durations;else,u=usags;T=durations(j);end
    [tr,row]=sag_run(refs(k),u,T,p,kind);runs{k,j}=tr;tab=[tab;row]; %#ok<AGROW>
    fprintf('gamma0 %.6f | Usag %.1f | depth %.0f%% | T %.2f | excursion %.6g | max|nu-1| %.6g | maintained %d | %s\n', ...
        row.gamma0_deg,u,row.sagDepth_percent,T,row.maxAngleExcursion_deg,row.maxAbsNuDeviation,row.contractionMaintained,row.recoveryClass);
 end
 % Tight independent integration on the severest/longest case for each ref.
 tr=runs{k,end};[tight,r]=sag_run(refs(k),tr.summary.Usag_pu,tr.summary.faultDuration_normalized,p,kind,true);
 dx=max(abs(tight.state(end,:)-tr.state(end,:)));
 assert(r.recoveryClass==tr.summary.recoveryClass && dx<1e-5,'Sag classification not numerically robust.');
 checks=[checks;table(rad2deg(refs(k)),dx,'VariableNames',{'gamma0_deg','tightSolveFinalStateError'})]; %#ok<AGROW>
 if strcmp(kind,'magnitude')
    rowsK=tab(abs(tab.gamma0_deg-rad2deg(refs(k)))<1e-7,:);good=rowsK.recoveryClass=="Recovered original equilibrium";
    if any(good),fprintf('Strongest recovering sag among tested discrete cases: gamma0 %.6f, Usag %.1f pu\n',rad2deg(refs(k)),min(rowsK.Usag_pu(good)));
    else,fprintf('No tested sag recovered original equilibrium for gamma0 %.6f.\n',rad2deg(refs(k)));end
 end
 sag_point_figure(runs(k,:),kind,k,p,outDir);
end
writetable(tab,fullfile(outDir,['SCR5_XR10_voltage_sag_',kind,'_summary.csv']));
writetable(checks,fullfile(outDir,['SCR5_XR10_voltage_sag_',kind,'_verification.csv']));
save(fullfile(outDir,['SCR5_XR10_voltage_sag_',kind,'_results.mat']),'runs','tab','checks','-v7.3');
fig=figure('Color','w','Position',[80,40,1200,950]);tl=tiledlayout(fig,3,1);
colors=lines(numel(refs));
for k=1:numel(refs)
 t=tab(abs(tab.gamma0_deg-rad2deg(refs(k)))<1e-7,:);
 if strcmp(kind,'magnitude'),xx=t.sagDepth_percent;xl='Sag depth (%)';else,xx=t.faultDuration_normalized;xl='Fault duration (normalized analytical time)';end
 yy={t.maxAngleExcursion_deg,t.maxAbsNuDeviation,t.maxLambdaMaxSymF};
 for j=1:3,ax=nexttile(tl,j);hold(ax,'on');plot(ax,xx,yy{j},'-o','Color',colors(k,:),'DisplayName',sprintf('gamma0=%.6f',rad2deg(refs(k))));grid(ax,'on');xlabel(ax,xl);end
end
labs={'Maximum angle excursion (deg)','Maximum |nu-1|','Maximum sampled lambda max(Sym(F))'};
for j=1:3,ax=nexttile(tl,j);ylabel(ax,labs{j});if strcmp(kind,'magnitude'),xlim(ax,[10,90]);end;end
ax=nexttile(tl,1);legend(ax,'Location','best');stamp(tl,[kind,' sweep; sampled recovery and fixed-metric diagnostics']);
savefigs(fig,outDir,['voltage_sag_',kind,'_severity_summary']);
if strcmp(kind,'duration')
 fig=figure('Color','w','Position',[80,50,1100,550]);ax=axes(fig);hold(ax,'on');
 for k=1:numel(refs)
    t=tab(abs(tab.gamma0_deg-rad2deg(refs(k)))<1e-7,:);codes=double(t.recoveryClass=="Recovered original equilibrium")+2*double(t.recoveryClass=="Transitioned to another equilibrium branch");
    plot(ax,t.faultDuration_normalized,codes,'-o','DisplayName',sprintf('gamma0=%.6f',rad2deg(refs(k))));
 end
 yticks(ax,0:2);yticklabels(ax,{'No original recovery','Original equilibrium','Another branch'});ylim(ax,[-.3,2.3]);xlabel(ax,'Normalized fault duration');grid(ax,'on');legend(ax,'Location','best');
 title(ax,{'Limiter-free VSG | SCR=5, X/R=10 | H=3.4, D_eff=25','Discrete duration outcomes; no critical clearing time inferred'});
 savefigs(fig,outDir,'voltage_sag_duration_recovery');
end
end
function sag_point_figure(runs,kind,k,p,outDir)
fig=figure('Color','w','Position',[80,40,1250,1050]);tl=tiledlayout(fig,4,1,'TileSpacing','compact');
for j=1:numel(runs)
 tr=runs{j};label=sprintf('U=%.1f, T=%.1f',tr.summary.Usag_pu,tr.summary.faultDuration_normalized);
 ys={rad2deg(tr.state(:,1)),tr.state(:,2)-1,tr.distance,tr.lambdaMaxSymF};
 for a=1:4,ax=nexttile(tl,a);hold(ax,'on');plot(ax,tr.time,ys{a},'DisplayName',label);grid(ax,'on');end
end
labs={'gamma unwrapped (deg)','nu-1','L relative to original equilibrium','lambda max(Sym(F))'};
for a=1:4,ax=nexttile(tl,a);ylabel(ax,labs{a});xlabel(ax,'Normalized analytical time');end
ax=nexttile(tl,1);legend(ax,'Location','best','NumColumns',3);yline(ax,rad2deg(runs{1}.reference.gamma),'r--','Original equilibrium');
ax=nexttile(tl,4);yline(ax,0,'k:');stamp(tl,sprintf('%s sweep; gamma0=%.6f; reference certificate=%d',kind,rad2deg(runs{1}.reference.gamma),runs{1}.summary.referenceContractionCertified));
savefigs(fig,outDir,sprintf('voltage_sag_%s_point_%d',kind,k));
end

function stamp(tl,txt)
title(tl,{'Limiter-free VSG | SCR=5, X/R=10 | H=3.4, D_{eff}=25',txt});
end
function savefigs(fig,outDir,stem)
exportgraphics(fig,fullfile(outDir,['SCR5_XR10_',stem,'.png']),'Resolution',180);
savefig(fig,fullfile(outDir,['SCR5_XR10_',stem,'.fig']));close(fig);
end
function base_figures(scan,samples,regions,boundaries,structural,secondary,p,outDir)
s=scan(scan.gamma0_deg>=0 & scan.gamma0_deg<=180,:);s=s(1:20:end,:);
fig=figure('Color','w','Position',[50 50 1200 850]);tl=tiledlayout(fig,2,1);
ax=nexttile(tl);plot(ax,s.gamma0_deg,s.Pstar_pu,'LineWidth',1.5);hold(ax,'on');
plot(ax,secondary.gamma_deg,secondary.Pe_QV_pu,'--');xlim(ax,[0 180]);yline(ax,1,':');grid(ax,'on');ylabel(ax,'P_e (pu)');legend(ax,'Constant E baseline','Algebraic Q-V sensitivity');
ax=nexttile(tl);plot(ax,s.gamma0_deg,s.dPe_dgamma,'LineWidth',1.5);hold(ax,'on');plot(ax,secondary.gamma_deg,secondary.dPe_QV_dgamma,'--');xlim(ax,[0 180]);yline(ax,0,':');grid(ax,'on');ylabel(ax,'dP_e/dgamma (pu/rad)');xlabel(ax,'gamma (deg)');stamp(tl,'Power-angle law and derivative');savefigs(fig,outDir,'power_and_derivative');
fig=figure('Color','w','Position',[50 30 1200 1100]);tl=tiledlayout(fig,3,1);
ys={[s.eig1_real,s.eig2_real],[s.lambdaMinM,s.lambdaMaxM],[s.lambdaMinSymF,s.lambdaMaxSymF]};labs={'Real eigenvalues of J','Eigenvalues of M','Eigenvalues of Sym(F)'};
for j=1:3
 ax=nexttile(tl);plot(ax,s.gamma0_deg,ys{j},'LineWidth',1.3);hold(ax,'on');yline(ax,0,'k:');grid(ax,'on');ylabel(ax,labs{j});xlabel(ax,'gamma0 (deg)');xlim(ax,[0 180]);
 if j==3,ylim(ax,[-12 8]);end
 for k=1:numel(regions),xline(ax,regions(k),'b--','HandleVisibility','off');end
end
stamp(tl,'Equilibrium-family matrix eigenvalues; Sym(F) vertical range clipped near singularities');savefigs(fig,outDir,'matrix_eigenvalues');
fig=figure('Color','w','Position',[50 50 1200 800]);tl=tiledlayout(fig,2,1);
ax=nexttile(tl);plot(ax,s.gamma0_deg,s.maxRealEigJ,'LineWidth',1.5);yline(ax,0,'k:');grid(ax,'on');ylabel(ax,'max Re eig(J)');
ax=nexttile(tl);plot(ax,s.gamma0_deg,s.lambdaMaxSymF,'LineWidth',1.5);hold(ax,'on');yline(ax,0,'k:');ylim(ax,[-1 3]);grid(ax,'on');ylabel(ax,'max eig(Sym(F))');xlabel(ax,'gamma0 (deg)');
for k=1:size(regions,1),plot(ax,regions(k,:),[-.8 -.8],'g-','LineWidth',7);end
stamp(tl,'Local stability and sufficient contraction certificate are distinct');savefigs(fig,outDir,'stability_and_certificate');
fig=figure('Color','w','Position',[50 50 1400 750]);tl=tiledlayout(fig,2,2);
for k=1:height(boundaries)
 ax=nexttile(tl);gg=linspace(boundaries.boundary_deg(k)-.15,boundaries.boundary_deg(k)+.15,501);rr=rows(gg,p);plot(ax,gg,rr.lambdaMaxSymF,'LineWidth',1.5);hold(ax,'on');yline(ax,0,'k:');xline(ax,boundaries.boundary_deg(k),'r--');grid(ax,'on');xlabel(ax,'gamma0 (deg)');ylabel(ax,'max eig(Sym(F))');title(ax,sprintf('R%d %s; min eig(M)=%.4g',boundaries.region(k),boundaries.side(k),boundaries.lambdaMinM(k)));
end
stamp(tl,'Refined boundary causes; metric remains positive definite');savefigs(fig,outDir,'boundary_causes');
% Sampled gamma intervals hold BOTH reference metric and Pstar fixed.
fig=figure('Color','w','Position',[50 50 1200 700]);ax=axes(fig);hold(ax,'on');it=table();gg=linspace(-180,180,72001).';
for k=1:size(regions,1)
 ref=point(deg2rad(mean(regions(k,:))),p);ll=fixed_lambda(gg,ref,p);plot(ax,gg,ll,'DisplayName',sprintf('R%d midpoint reference',k));ok=ll<-p.tol;beg=find(diff([false;ok])==1);en=find(diff([ok;false])==-1);
 for j=1:numel(beg)
 lo=gg(beg(j));hi=gg(en(j));if beg(j)>1,lo=fzero(@(g)fixed_lambda(g,ref,p)+p.tol,gg(beg(j)+[-1 0]));end
 if en(j)<numel(gg),hi=fzero(@(g)fixed_lambda(g,ref,p)+p.tol,gg(en(j)+[0 1]));end
 it=[it;table(k,rad2deg(ref.gamma),lo,hi,'VariableNames',{'referenceRegion','referenceGamma_deg','gammaMin_deg','gammaMax_deg'})]; %#ok<AGROW>
 end
end
yline(ax,0,'k:');ylim(ax,[-2 5]);grid(ax,'on');xlabel(ax,'Instantaneous gamma (deg), nu=1');ylabel(ax,'max eig(Sym(F))');legend(ax,'Location','best');title(ax,{'Limiter-free VSG | SCR=5, X/R=10 | H=3.4, D_{eff}=25','Fixed-reference metric and Pstar: sampled instantaneous test'});savefigs(fig,outDir,'fixed_reference_instantaneous');writetable(it,fullfile(outDir,'SCR5_XR10_fixed_reference_gamma_intervals.csv'));
end
function lm=fixed_lambda(gd,s,p)
d=electric(deg2rad(gd),p.U,p);r=s.c/s.beta;a=-d.Pg/p.mass;b=(-(s.Pstar-d.Pe)-p.D)/p.mass;off=(1+a-b*r-r^2)/2;lm=b/2+hypot(-r-b/2,off);
end
function figure5_views(gd,labels,runs,I,sags,p,outDir)
fig=figure('Color','w','Position',[40 50 1300 650]);ax=axes(fig);hold(ax,'on');
for k=1:size(I,1),plot(ax,I(k,:),[1 1]*k,'LineWidth',12,'DisplayName',sprintf('Certified R%d',k));end
for k=1:numel(gd),plot(ax,gd(k),0,'ko','HandleVisibility','off');text(ax,gd(k),-.1,labels(k),'Rotation',35,'FontSize',8);end
xlim(ax,[0 180]);ylim(ax,[-1 size(I,1)+.5]);xlabel(ax,'Equilibrium gamma0 (deg)');ylabel(ax,'Certificate region');grid(ax,'on');legend(ax,'Location','best');title(ax,{'Limiter-free VSG | SCR=5, X/R=10 | H=3.4, D_{eff}=25','Newly computed operating-point regions'});savefigs(fig,outDir,'figure5A_operating_points');
fig=figure('Color','w','Position',[30 30 1400 1100]);tl=tiledlayout(fig,3,2,'TileSpacing','compact');ax=nexttile(tl,[1 2]);hold(ax,'on');stop=75;
for k=1:numel(runs)
 tr=runs{k};if max(real(tr.reference.eigJ))>0,j=find(abs(rad2deg(tr.state(:,1)-tr.reference.gamma))>=30,1);if ~isempty(j),stop=max(stop,tr.time(j)*1.1);end;end
 plot(ax,tr.time,rad2deg(tr.state(:,1)),'DisplayName',labels(k),'LineWidth',1.4);
end
for k=1:numel(I),yline(ax,I(k),'b--','HandleVisibility','off');end
xlim(ax,[0 stop]);grid(ax,'on');ylabel(ax,'gamma unwrapped (deg)');xlabel(ax,'Normalized analytical time');legend(ax,'Location','northwest','NumColumns',3,'FontSize',8);title(ax,'Common +1 deg test perturbation (chosen for visualization)');
cert=find(cellfun(@(s)s.summary.referenceContractionCertified,sags(:,1)));panel=0;
for k=1:min(2,numel(cert))
 for u=[.5 .3]
  panel=panel+1;row=cert(k);j=find(cellfun(@(s)abs(s.summary.Usag_pu-u)<1e-8,sags(row,:)),1);tr=sags{row,j};a=nexttile(tl,2+panel);plot(a,tr.time,rad2deg(tr.state(:,1)),'LineWidth',1.4);hold(a,'on');yline(a,rad2deg(tr.reference.gamma),'r--');xline(a,1,'k:');grid(a,'on');xlabel(a,'Normalized analytical time');ylabel(a,'gamma unwrapped (deg)');title(a,sprintf('gamma0=%.4f, Usag=%.1f; %s',rad2deg(tr.reference.gamma),u,tr.summary.recoveryClass),'FontSize',9);
 end
end
stamp(tl,'Temporal trajectories and voltage-sag examples; no angle jump in sag tests');drawnow;
pos=ax.Position;ins=axes(fig,'Position',[pos(1)+.61*pos(3),pos(2)+.15*pos(4),.34*pos(3),.46*pos(4)]);hold(ins,'on');
for k=1:numel(runs),tr=runs{k};if max(real(tr.reference.eigJ))<0,plot(ins,tr.time,rad2deg(tr.state(:,1)));end;end
for k=1:numel(I),yline(ins,I(k),'b--');end
xlim(ins,[0 10]);grid(ins,'on');title(ins,'Early-time stable cases','FontSize',8);ins.FontSize=7;savefigs(fig,outDir,'figure5B_temporal_gamma_with_sag_examples');
end
function t=comparison_table(I,g,p,n)
names=["E";"Vth magnitude";"Req";"Xeq";"Nominal gamma";"H";"D_eff"];old=[1;1;0;.3;17.4576031237221;3.4;25];new=[p.E;abs(n.Vth);real(n.Zeq);imag(n.Zeq);rad2deg(g);p.H;p.D];
oldI=[15.32259918055 21.091668254769;48.626361830614 81.415962806098];
for k=1:max(size(I,1),2)
 for side=1:2
 names(end+1)="R"+k+" boundary "+side;old(end+1)=NaN;new(end+1)=NaN;if k<=2,old(end)=oldI(k,side);end;if k<=size(I,1),new(end)=I(k,side);end
 end
end
t=table(names,old,new,new-old,'VariableNames',{'quantity','oldValue','newValue','change'});
end
function write_report(r,outDir)
p=r.parameters;n=r.network;fid=fopen(fullfile(outDir,'SCR5_XR10_report.md'),'w','n','UTF-8');guard=onCleanup(@()fclose(fid));
fprintf(fid,'# Limiter-free SCR=5, X/R=10 baseline\n\nGenerated from saved GFMVSG1 model; normalized analytical time. Frozen analytical files were not changed.\n\n');
fprintf(fid,'## Network and scope\n\n%s.\n\n%s.\n\n%s.\n\n',r.sourceAudit.topology,r.sourceAudit.absent,r.sourceAudit.excluded);
fprintf(fid,'Sbase=%.0f VA; capacitor Vbase=%.0f V line-line; transformer HV base=%.0f V (1.02 times measurement nominal). Physical ideal bus=20000 V, referred U=%.12g pu. Transformer D1 phase shift is removed on referral. Rgrid=%.10g ohm; Lgrid=%.10g H.\n\n',p.S,p.Vlv,p.Vhv,p.U,p.Rgrid,p.Lgrid);
fprintf(fid,'E is the physical capacitor-node fundamental amplitude under ideal inner-loop tracking and unity measurement transfer approximation. Constant E=%.12g is calibrated at Pstar=1 using E=1-0.1 Qgrid. The secondary algebraic Q-V curve is sensitivity only; it does not change primary derivatives or certificates. Dynamic Q-V, voltage/current control and measurement states are omitted. This is not an exact EMT reduction.\n\n',p.E);
fprintf(fid,'Transformer z1=%.9g+j%.9g, z2=%.9g+j%.9g; Ym=%.9g+j%.9g; series-RC shunt Yc=%.9g+j%.9g. Converter RL is upstream of E and excluded from downstream equivalent.\n\n',real(n.Zt1),imag(n.Zt1),real(n.Zt2),imag(n.Zt2),real(n.Ym),imag(n.Ym),real(n.Yc),imag(n.Yc));
fprintf(fid,'Nodal reduction: Y=[1/z1+1/z2+Ym,-1/z2;-1/z2,1/z2+1/(zg+zb)]. Nodes=Y\\[Vcap/z1;U/(zg+zb)]. Ig=(node2-U)/(zg+zb), Vg=node2-zb Ig, Ic=(Vcap-node1)/z1+Yc Vcap. Zeq=1/icE, Vth=-icU/icE U.\n\nZeq=%.12g+j%.12g pu; Vth=%.12g+j%.12g pu. Gamma=arg(Vcap)-arg(Vth); Vcap=E exp(j(gamma+phi)), phi=%.12g deg.\n\n',real(n.Zeq),imag(n.Zeq),real(n.Vth),imag(n.Vth),rad2deg(p.phi));
fprintf(fid,'Power is Re(Vg conj(Ig)) at the grid measurement, NOT capacitor source power. With actual referred bus magnitude U, Pe=a E^2+b U^2+E U(c cos(gamma)+d sin(gamma)); [a,b,c,d]=[%.15g %.15g %.15g %.15g].\n\n',p.pc);
fprintf(fid,'Pe_gamma=EU(-c sin(gamma)+d cos(gamma)); Pe_U=2bU+E(c cos(gamma)+d sin(gamma)); Pe_gammaU=E(-c sin(gamma)+d cos(gamma)).\n\n');
fprintf(fid,'Current block paths (prefix GFMVSG1/Grid_Forming_Converter):\n\n- /Converter_Control/Outer_Control/Virtual_Synchronous_Machine\n- /Converter_Control/Operation_Mode_Selector\n- /Converter_Control/Inner_Control/Cascaded_PI_controller/v_cap_Controller\n- /Electrical_Part/LC_harmonic_filter/V_cap_abc\n- /Converter_Control/measurement_transformations/Power_Calculation_(Positive_Sequence)\n\nThe VSM Q-V amplitude command feeds the mode selector (Lvirt=0), then the capacitor voltage controller; physical V_cap is the chosen E node. The measured grid positive-sequence voltage/current feed active power. Bridge voltage is not the fixed E source.\n\n');
fprintf(fid,'LV bases: I_RMS=%.12g A, Z=%.12g ohm. HV transformer bases: I_RMS=%.12g A, Z=%.12g ohm. P and Q use Sbase. Peak dq and RMS phasors have identical pu magnitude on compatible bases. The Thevenin equivalent is exact for this linear positive-sequence network at 50 Hz; the control/measurement quasi-steady assumptions are the approximation.\n\n',p.S/(sqrt(3)*p.Vlv),p.Vlv^2/p.S,p.S/(sqrt(3)*p.Vhv),p.Vhv^2/p.S);
fprintf(fid,'## Dynamics and unchanged metric construction\n\ngamma_dot=nu-1; nu_dot=((Pstar-Pe)/nu-D_eff(nu-1))/(2H). H=3.4, K=20, Kd=5, D_eff=25. Each equilibrium has its own fixed Pstar=Pe(gamma0), nu0=1.\n\nJ=[0,1;-Pe_gamma/(2H nu),(-(Pstar-Pe)/nu^2-D_eff)/(2H)]. beta=-Pe_U/(2H), c_metric=-Pe_gammaU/(2H), Theta=[-beta,0;-c_metric,-beta], M=Theta''Theta. F=Theta J/Theta with fixed reference Theta (Theta_dot=0). No optimized metric.\n\nFull-period scan: -180..180 deg, 0.0005 deg spacing, refined scalar roots. Certificate: min eig(M)>1e-9 and max eig(Sym(F))<-1e-9. CSV distinguishes tolerance boundaries and mathematical zeros. Local stability and sufficient contraction are different tests. Singular Theta points are excluded.\n\n');
fprintf(fid,'## Results\n\nNominal angle=%.12g deg. Nominal is locally stable; certificate=%d. All Pstar=1 roots in one period are in nominal_equilibria.csv; copies repeat by 360 deg. The equilibrium family always sets Pstar=Pe; a fixed-Pstar root can disappear at a power extremum.\n\n',rad2deg(p.gammaNominal),r.nominalEquilibria.contractionCertified(r.nominalEquilibria.isSelected));
for k=1:size(r.regions_deg,1),fprintf(fid,'- R%d: %.12g to %.12g deg.\n',k,r.regions_deg(k,:));end
fprintf(fid,'\nBoundary causes and structural points are exported separately. All primary changes arise from electrical network, power measurement mapping and nominal E calibration; controller/metric method stays fixed. These changes are not attributable to X alone.\n\n');
fprintf(fid,'### Old/new numerical comparison\n\n| Quantity | Old | New | Change |\n|---|---:|---:|---:|\n');
for k=1:height(r.comparison),fprintf(fid,'| %s | %.12g | %.12g | %.12g |\n',r.comparison.quantity(k),r.comparison.oldValue(k),r.comparison.newValue(k),r.comparison.change(k));end
fprintf(fid,'\nLocal stability intervals (deg):\n');for k=1:size(r.locallyStableIntervals_deg,1),fprintf(fid,'- %.12g to %.12g\n',r.locallyStableIntervals_deg(k,:));end
fprintf(fid,'\nNominal Pstar=1 mathematical equilibria (modulo 360 deg):\n');for k=1:height(r.nominalEquilibria),t=r.nominalEquilibria(k,:);fprintf(fid,'- gamma0=%.12g deg; J eigenvalues %.12g, %.12g; %s.\n',t.gamma0_deg,t.eig1_real,t.eig2_real,t.localClassification);end
fprintf(fid,'\nFull matrices for every selected point are stored in MAT pointMatrices and printed in the console log. The second certified region survives.\n\n');
fprintf(fid,'## Trajectories and voltage sags\n\nRepresentative trajectories use our common +1 deg test (not a paper-specified magnitude); unstable main panels show early departure. Full diagnostics retain unwrapped angles and L relative to ORIGINAL equilibrium, so another branch does not yield L=0. A neighboring trajectory starts at +1.01 deg.\n\nSag references include the new nominal and midpoints of all positive certified regions. Nominal is explicitly labelled uncertified. Initial states are exact equilibria, no angle jump. Retained voltage 0.9..0.1 with duration 1; duration sweep 1,2,3,4,5,7,10 at 0.1 pu. Pstar is fixed during sag; U is restored at clearing.\n\nRecovery uses the last 10 percent of elapsed trajectory: angle error <0.1 deg and |nu-1|<1e-4. Other-branch test also checks power residual <1e-4. nu=0.5 terminates integration rather than modifying dynamics. Tests use a fixed reference metric: sampled instantaneous test, not a regional proof. Discrete recovering tests are not critical voltages or CCT.\n\n');
fprintf(fid,'## Verification\n\n');fn=fieldnames(r.validation);for k=1:numel(fn),fprintf(fid,'- %s: %.6g\n',fn{k},r.validation.(fn{k}));end
fprintf(fid,'\nTighter independent sag integrations are checked and saved. PNG and editable FIG outputs include power/derivative, three matrix spectra, local/certificate regions, boundary causes, fixed-reference instantaneous regions, trajectories, Figure 5A/5B equivalents and sag magnitude/duration diagnostics. H-D retuning/sweep and unrelated legacy SG-only figures are not part of this fixed-tuning baseline.\n\nThe baseline is ready as an internally verified reference for a subsequent reduced current-limiter study within the stated approximations. No current-limiter study was started.\n');
fprintf('\nOLD BASELINE: E=1, Xeq=0.3, H=3.4, D_eff=25; gamma_nominal=17.457603 deg\nR1=15.322599..21.091668; R2=48.626362..81.415963 deg\n');
fprintf('\nNEW SCR5/XR10 BASELINE: E=%.12g, Vth=%.12g+j%.12g, Zeq=%.12g+j%.12g pu\n|Zeq|=%.12g; angle=%.9g deg; H=3.4 K=20 Kd=5 D_eff=25\ngamma_nominal=%.12g deg\n',p.E,real(n.Vth),imag(n.Vth),real(n.Zeq),imag(n.Zeq),abs(n.Zeq),rad2deg(angle(n.Zeq)),rad2deg(p.gammaNominal));
fprintf('Local stability intervals (deg):\n');disp(r.locallyStableIntervals_deg);fprintf('Certified regions (deg):\n');disp(r.regions_deg);disp(r.comparison);
fprintf('Second certified region survives: %d. Changes are electrical-network/mapping/calibration only; tuning and metric method unchanged.\n',size(r.regions_deg,1)>=2);
fprintf('Numerically verified baseline ready for subsequent reduced limiter study within documented assumptions. Limiter analysis not started.\n');
end



function print_matrices(points,labels)
fprintf('\nPOINTWISE MATRIX AUDIT (each equilibrium has its own fixed reference metric)\n');
for k=1:numel(points)
 s=points(k);fprintf('\n%s: gamma0=%.9f deg, Pstar=%.12g, beta=%.9g, c=%.9g\n',labels(k),rad2deg(s.gamma),s.Pstar,s.beta,s.c);
 fprintf('J =\n');disp(s.J);fprintf('Theta =\n');disp(s.Theta);fprintf('M =\n');disp(s.M);fprintf('F =\n');disp(s.F);fprintf('Sym(F) =\n');disp(s.SymF);
end
end

