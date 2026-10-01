% Isolated numerical diagnostics reuse prior sweeps and preserve final figures.
if strcmp(getenv('GFM_HD_SANITY_ONLY'),'1')
    hdOut=fileparts(mfilename('fullpath'));
    hdSaved=load(fullfile(hdOut,'HD_eff_sweep_contraction_boundaries.mat'),'sweep');
    hd_boundary_sanity(hdSaved.sweep,hdOut);
    return
end
% Legacy plotting-only request intentionally creates no figures.
if strcmp(getenv('GFM_HD_PLOT_ONLY'),'1')
    fprintf('H-D_eff visualization disabled; numerical data retained.\n');
    return
end
% Optional isolated execution of the additive sweep preserves existing exports.
if strcmp(getenv('GFM_HD_SWEEP_ONLY'),'1')
    hd_parameter_sweep(fileparts(mfilename('fullpath')));
    return
end
%% CURRENT-LIMITER-FREE RETUNED VSM-GFM BASELINE
% Candidate contraction-oriented retuning; no optimization is claimed.
% Revised IN PLACE. No resistance, limiter, sag, inner control or LCL states.
% Pourazad synchronism workflow, with VSM dynamics in the requested normalized analytical convention.
% Main scan: gamma=gamma0. gamma_rep=0 is a SEPARATE legacy-SG diagnostic.
% Method A uses c_paper=omega0_pu*d(beta)/d(gamma), omega0_pu=1 to DESIGN a metric at gamma0.
% c_paper is the paper's beta_dot convention, NOT the true VSM beta derivative.
% Once chosen for a reference, Theta_ref is fixed during its trajectories.
% Eq15c/2 equals Sym(Theta_ref*A(gamma)/Theta_ref). The complete direct
% definition includes Theta_dot_ref=0 because the reference is fixed.
% Also audit (i) Method B's actual derivatives and its different metric,
% and (ii) a moving-angle literal derivative gamma_dot=omega0_pu. The latter
% is NOT a synchronous VSM equilibrium and is never used to claim a certificate.
clearvars; clc; close all; format long g
outDir=fileparts(mfilename('fullpath'));if isempty(outDir),outDir=pwd;end
logFile=fullfile(outDir,'baseline_gfm_diagnostics.txt');
if isfile(logFile),delete(logFile);end
diary(logFile);

%% 1. Normalization audit -- BEFORE optional SG reproduction
f_base=50;omega_base=2*pi*f_base;omega0_pu=1.0;
H=3.4;K=20;Kd=5;D=K+Kd;
E=1.0;U1=1.0;Xeq=0.3;Pstar_nominal=1.0;
% MAIN analytical state is x=[gamma;nu], NOT [gamma;Delta_omega_SI].
% Requested paper numerical convention: gamma_dot=nu-1.
% Physical EMT time would require gamma_dot=omega_base*(nu-1).
% Changing ONLY the angle equation is not an equivalent time/coordinate
% transformation of the physical EMT system. A common time rescaling would
% also rescale the frequency equation. We simulate the specified normalized
% analytical model, not claim identical physical EMT trajectories.
par=struct('f_base',f_base,'omega0_pu',omega0_pu,'H',H,'K',K,'Kd',Kd,'D',D, ...
    'mass',2*H,'E',E,'U1',U1,'Xeq',Xeq,'Pstar_nominal',Pstar_nominal, ...
    'idealPLL',true,'frequencyReference_pu',1,'powerBase_VA',2e6);
% omega_base is deliberately NOT a field of the main-model parameter struct.
C=E/(2*H*Xeq);b=D/(2*H);Pmax=E*U1/Xeq;
assert(omega0_pu==1 && H==3.4 && K==20 && Kd==5 && D==25);
fprintf('====================================================\nNORMALIZATION AUDIT\n====================================================\n');
fprintf('Physical EMT base:\nf_base = %.0f Hz\nomega_base = %.9f rad/s\n',f_base,omega_base);
fprintf('Pourazad-style analytical convention:\nomega0_pu = %.1f pu\n',omega0_pu);
fprintf('VSM parameters:\nH = %.1f s\nK = %.0f pu\nKd = %.0f pu\nD_eff = %.0f pu\n',H,K,Kd,D);
fprintf(['Main contraction state: x=[gamma;nu]\ngamma_dot=nu-1; nu0=1 pu\n' ...
    'Physical mapping retained only for reference: omega_phys=omega_base*nu\n' ...
    'NO 314.159 factor is used inside beta2_dot, Theta or Eq15c of the main branch.\n' ...
    'Analytical time is the specified normalized-model time, NOT physical EMT seconds.\n' ...
    'This is NOT a coordinate-only re-expression of the physical EMT dynamics.\n']);
fprintf('Candidate contraction-oriented retuning: D/(2H)=%.12g; D/H=%.12g; Pmax=%.12g pu\n',b,D/H,Pmax);
fprintf('SI equivalents for reference ONLY: Mv_SI=%.12g, Dv_SI=%.12g\n',2*H/omega_base,D/omega_base);
fprintf('====================================================\n');
sgValidation=legacy_sg_verification(outDir);
assert(sgValidation.pass,'Legacy SG algebra/classification reproduction failed.');
fprintf('\nCORRECTED PU VSM-GFM CONTRACTION RESULTS\n');disp(par);
fprintf(['Pe=(E*U1/Xeq)*sin(gamma)\n' ...
    'f=[nu-1;((Pstar-Pe)/nu-D*(nu-1))/(2*H)]\n' ...
    'A=[0,1;-E*U1*cos(gamma)/(2*H*Xeq*nu),(-(Pstar-Pe)/nu^2-D)/(2*H)]\n' ...
    'B=[0;-E*sin(gamma)/(2*H*Xeq*nu)]\n' ...
    'At equilibrium: alpha1=-E*U1*cos(gamma0)/(2*H*Xeq); alpha2=-D/(2*H) (printed numerically below)\n' ...
    'beta2=-E*sin(gamma)/(2*H*Xeq*nu)\n' ...
    'beta2_dot_paper=-E*omega0_pu*cos(gamma0)/(2*H*Xeq), omega0_pu=1\n' ...
    'beta2_dot_paper/beta2=cot(gamma0)\n' ...
    'beta2_dot_coordinate=-C*cos(gamma)*(nu-1)/nu+C*sin(gamma)*nu_dot/nu^2\n' ...
    'The full coordinate derivative is zero at a true equilibrium.\n' ...
    'Main branch: paper-literal PU contraction branch, SymF=Eq15c/2.\n' ...
    'Same frozen-reference-metric direct definition is checked independently.\n']);
fprintf('alpha2=%.12g\n',-b);
fprintf('Nominal Pstar=1 pu angle=%.12f deg\n',asind(Pstar_nominal/Pmax));

%% 3. Dense equilibrium scan; all Pstar vary with gamma0
% Step <=0.0005 deg resolves even the narrow second component near 90 deg.
% All mathematical powers are retained; the <=1 pu sector is marked in plots.
gamma0_deg=unique(sort([linspace(0.001,119.999,240001),90,asind(Pstar_nominal/Pmax)]));
gamma0=deg2rad(gamma0_deg);N=numel(gamma0);
scan=scan_model(gamma0,par);
intervals=refine_intervals(gamma0,scan.lambdaPaper,@(g) lambda_at(g,par,'paper'));
intervalsDirect=refine_intervals(gamma0,scan.lambdaDirect,@(g) lambda_at(g,par,'direct'));
intervalsCoordinate=zeros(0,2);intervalsMoving=zeros(0,2);
% Both other constructions have SymF(1,1)=0 analytically => cannot be strictly ND.
assert(all(scan.lambdaCoordinate>=-1e-8) && all(scan.lambdaMoving>=-1e-7));
assert(isequal(size(intervals),size(intervalsDirect)) && ...
    max(abs(intervals-intervalsDirect),[],'all')<1e-8);
assert(max(scan.equilibriumResidual)<1e-10);
assert(max(scan.directResidual)<1e-8);
assert(all(scan.minEigM>0) && all(scan.detM>0));
print_intervals('Paper-style VSM-GFM contraction interval (gamma=gamma0)',intervals);
print_intervals('Direct-definition interval, same frozen paper metric',intervalsDirect);
boundaryTable=boundary_causes(intervals,par);
fprintf('Boundary cause diagnostics (half-symmetric convention):\n');disp(boundaryTable);
writetable(boundaryTable,fullfile(outDir,'boundary_causes.csv'));
fullSylvester=scan.s11<0 & scan.s22<0 & scan.detSymF>0 & scan.minEigM>0;
assert(isequal(fullSylvester,scan.paperCertified),'Full determinant/eigenvalue tests disagree.');
fprintf('Method B actual-derivative metric interval: EMPTY (SymF11=0).\n');
fprintf('Literal moving-angle derivative interval: EMPTY (SymF11=0); not a VSM equilibrium.\n');
fprintf('Necessary diagonal sector only: %.10f < gamma0 < 90 deg\n',rad2deg(atan(omega0_pu/b)));
fprintf('0 <= Pstar <= 1 pu in scanned sector: 0 <= gamma0 <= %.9f deg\n',asind(1/Pmax));
fprintf('Contraction interval powers at boundaries [pu]:\n');disp(Pmax*sin(intervals));
fprintf('Maximum full direct/paper residual over scan: %.3e\n',max(scan.directResidual));
fprintf('Maximum equilibrium residual over scan: %.3e\n',max(scan.equilibriumResidual));

%% 4. Select true cases by eigenvalue properties, not assigned angle labels
nominal=asin(Pstar_nominal/Pmax);snom=point_model(nominal,par);
nominalStable=max(real(snom.eigA))<0;
nominalCertified=min(snom.eigM)>0 && max(snom.eigPaper)<0;
fprintf('\nNominal 1-pu operating point: gamma0=%.12f deg; residual=%.3e\n',rad2deg(nominal),snom.equilibriumResidual);
for fld={'A','eigA','M','eigM','SymF_paper','eigPaper'}
    fprintf('%s=\n',fld{1});disp(snom.(fld{1}));
end
fprintf('lambda_max(SymF)=%.12g; det(SymF)=%.12g\n',max(snom.eigPaper),det(snom.SymF_paper));
yesno={'NO','YES'};
fprintf('LOCALLY STABLE = %s\nCONTRACTION CERTIFIED = %s\n',yesno{1+nominalStable},yesno{1+nominalCertified});
% These are 12 different steady-state operating points, NOT disturbance magnitudes.
% Each has Pstar=Pe(gamma0). Five representatives avoid crowded figures;
% all twelve remain stored. A fixed-point severity study is deferred.
tol=1e-9;
testAngles=[10,12,15.8,17.457603,19,20.5,25,35,55,65,100,120];
regionLabels=[repmat("Below Region 1",1,2),repmat("Inside Region 1",1,4), ...
    repmat("Between Regions",1,2),repmat("Inside Region 2",1,2),repmat("Above Region 2",1,2)];
pointResults=table();pointModels=cell(1,numel(testAngles));
for k=1:numel(testAngles)
    q=point_model(deg2rad(testAngles(k)),par);ma=max(real(q.eigA));mm=min(q.eigM);mf=max(q.eigPaper);
    if ma < -tol,st="Locally stable";elseif ma>tol,st="Locally unstable";else,st="Marginal / boundary";end
    q.locallyStable=ma < -tol;q.contractionCertified=mm>tol && mf < -tol;
    q.stabilityClass=st;q.regionLabel=regionLabels(k);pointModels{k}=q;
    row=table(testAngles(k),q.Pstar,q.equilibriumResidual,ma,mm,mf,det(q.SymF_paper), ...
        q.locallyStable,q.contractionCertified,regionLabels(k),st,{q.A},{q.eigA},{q.M},{q.eigM},{q.SymF_paper},{q.eigPaper}, ...
        'VariableNames',{'gamma0_deg','Pstar_pu','eq_residual','maxRealEigA','lambdaMinM','lambdaMaxSymF','detSymF', ...
        'locallyStable','contractionCertified','regionLabel','stabilityClass','A','eigA','M','eigM','SymF','eigSymF'});
    pointResults=[pointResults;row]; %#ok<AGROW>
end
assert(height(pointResults)==12 && all(pointResults.eq_residual<1e-10));
fprintf('\n========================================================\n12-POINT VSM-GFM OPERATING-POINT CLASSIFICATION\n========================================================\n');
fprintf('H=%.1f s; K=%.0f; Kd=%.0f; D_eff=%.0f; classification tol=%.3g\n',H,K,Kd,D,tol);
print_intervals('Verified contraction intervals',intervals);disp(pointResults(:,1:11));
writetable(pointResults(:,1:11),fullfile(outDir,'12_point_gfm_classification.csv'));
save(fullfile(outDir,'12_point_gfm_classification.mat'),'pointResults','pointModels','tol');
representativeIndices=[2,4,8,10,12];selected=deg2rad(testAngles(representativeIndices));
caseNames={'Below R1 - 12 deg','Inside R1 - nominal 1 pu','Between R1/R2 - 35 deg','Inside R2 - 65 deg','Above R2 - 120 deg'};
caseLabels=cellstr(pointResults.stabilityClass(representativeIndices));
fprintf('Five representative trajectory angles (deg):\n');disp(testAngles(representativeIndices));
fprintf(['COMMON +1 deg test perturbation: gamma(0)=gamma0+1 deg; nu(0)=1 for all five representatives.\n' ...
    'This is our own visualization/test choice, NOT a magnitude taken from the reference paper.\n' ...
    'Equilibrium and metric reference remain xd=[gamma0;1]. No Figure-5 severity analysis is implemented.\n']);
% Previously verified OLD PU baseline summary, never used in the new scan.
old=struct('H',.1,'K',20,'Kd',40,'D',60,'intervals_deg',[89.6730509036269,89.8749234883390], ...
    'nominalLambda',337.117533861893,'nominalCertified',false,'source','Previously verified OLD PU baseline; not recalculated');
comparison=struct('old',old,'newIntervals_deg',rad2deg(intervals),'newNominalLambda',max(snom.eigPaper), ...
    'newNominalCertified',nominalCertified,'widthChange_deg',sum(rad2deg(diff(intervals,1,2)))-diff(old.intervals_deg), ...
    'nominalLambdaChange',max(snom.eigPaper)-old.nominalLambda);
fprintf('\n====================================================\nCONTROLLER TUNING COMPARISON\n====================================================\n');
fprintf('OLD: H=0.1 K=20 Kd=40 D=60; nominal certified=NO\n%s\n',old.source);
print_intervals('OLD interval',deg2rad(old.intervals_deg));
fprintf('NEW: H=%.1f K=%.0f Kd=%.0f D=%.0f; nominal certified=%s\n',H,K,Kd,D,yesno{1+nominalCertified});
print_intervals('NEW intervals',intervals);
fprintf('Total width: OLD %.12g deg; NEW %.12g deg; change %+.12g deg\n',diff(old.intervals_deg),sum(rad2deg(diff(intervals,1,2))),comparison.widthChange_deg);
fprintf('Nominal lambda_max: OLD %.12g; NEW %.12g; change %+.12g\n',old.nominalLambda,max(snom.eigPaper),comparison.nominalLambdaChange);
intervalTable=table();
for k=1:size(intervals,1)
    gd=rad2deg(intervals(k,:));pw=power_smib(intervals(k,:),par);
    row=table(k,gd(1),gd(2),pw(1),pw(2),diff(gd),boundaryTable.Cause(2*k-1),boundaryTable.Cause(2*k), ...
        'VariableNames',{'Interval','gamma_min_deg','gamma_max_deg','Pstar_min_edge','Pstar_max_edge','width_deg','lowerCause','upperCause'});
    intervalTable=[intervalTable;row]; %#ok<AGROW>
end
fprintf('Full interval summaries:\n');disp(intervalTable);
writetable(intervalTable,fullfile(outDir,'contraction_intervals.csv'));
available=isfinite(selected);cases=cell(1,5);trajectories=cell(1,5);
T=table();
for j=find(available)
    s=pointModels{representativeIndices(j)};s.name=caseNames{j};s.classification=caseLabels{j};cases{j}=s;
    fprintf('\n%s: %s\n',s.name,s.classification);
    fprintf('gamma0=%.12f deg; Pstar=%.12g pu; above 1 pu=%d; residual=%.3e\n', ...
        rad2deg(s.gamma0),s.Pstar,s.Pstar>1,s.equilibriumResidual);
    fprintf('A=\n');disp(s.A);fprintf('B=\n');disp(s.B);fprintf('T Eq12=\n');disp(s.T);
    fprintf('alpha1=%.12g; alpha2=%.12g; beta2=%.12g\n',s.alpha1,s.alpha2,s.beta);
    fprintf('beta2_dot_paper=%.12g; beta2_dot_direct=%.12g\n',s.betaDotPaper,s.betaDotCoordinate);
    fprintf('Theta_paper=\n');disp(s.Theta);fprintf('Theta_dot_reference=\n');disp(s.ThetaDotReference);
    fprintf('M=\n');disp(s.M);fprintf('det(M)=%.12g; symmetry residual=%.3e\n',det(s.M),norm(s.M-s.M','fro'));
    fprintf('F_direct=\n');disp(s.Fdirect);fprintf('SymF_paper=\n');disp(s.SymF_paper);
    fprintf('SymF_direct=\n');disp(s.SymF_direct);fprintf('SymF_coordinate (Method B)=\n');disp(s.SymF_coordinate);
    fprintf('SymF_literal_moving (not equilibrium)=\n');disp(s.SymF_literal_moving);
    fprintf('eig(A)=\n');disp(s.eigA);fprintf('eig(M)=\n');disp(s.eigM);
    fprintf('eig(SymF_paper)=\n');disp(s.eigPaper);fprintf('eig(SymF_direct)=\n');disp(s.eigDirect);
    fprintf('lambda_max paper=%+.12g; direct=%+.12g; coordinate=%+.12g; moving=%+.12g\n', ...
        max(s.eigPaper),max(s.eigDirect),max(s.eigCoordinate),max(s.eigMoving));
    fprintf('norm(SymF_paper-SymF_direct)=%.3e\n',norm(s.SymF_paper-s.SymF_direct,'fro'));
    fprintf('norm(SymF_paper-SymF_literal_moving)=%.12g\n',norm(s.SymF_paper-s.SymF_literal_moving,'fro'));
    fprintf('gamma_rep=0 diagnostic lambda_max=%+.12g (NOT main classification)\n',max(s.eigRepresentative));
    row=table(string(s.name),rad2deg(s.gamma0),s.Pstar,s.Pstar>1,s.equilibriumResidual, ...
        s.eigA(1),s.eigA(2),s.eigM(1),s.eigM(2),max(s.eigPaper),max(s.eigDirect), ...
        max(s.eigCoordinate),max(s.eigMoving),max(s.eigRepresentative),string(s.classification), ...
        'VariableNames',{'Case','gamma0_deg','Pstar_pu','AboveOnePu','EquilibriumResidual', ...
        'eigA_1','eigA_2','eigM_1','eigM_2','lambdaPaper','lambdaDirect','lambdaCoordinate', ...
        'lambdaLiteralMoving','lambdaRepresentative','Classification'});
    T=[T;row]; %#ok<AGROW>
    trajectories{j}=simulate_case(s,par,j);
end

%% 5. Analytical-model validation (finite differences are tests only)
validation=validate_model(par,cases,available);
fprintf('\nAnalytical Jacobian relative check: %.3e; input Jacobian check: %.3e\n',validation.jacobian,validation.input);
fprintf('Actual beta derivative check: %.3e; speed-coordinate equivalence: %.3e\n',validation.beta,validation.coordinates);

%% 6. Figures: requested eigenvalues, consistency, power, and trajectories
colors=[.5,.25,.75;0,.35,.9;.9,.5,.05;0,.55,.2;.85,.1,.1];markers={'s','o','^','d','x'};
fig=figure('Color','w','Position',[50,80,1550,570]);tl=tiledlayout(fig,1,3,'TileSpacing','compact');
fields={'eigA','eigM','eigPaper'};titles={'Eigenvalues of A','Eigenvalues of M','Eigenvalues of Sym(F)'};
for n=1:3
    ax=nexttile(tl);hold(ax,'on');
    for j=find(available)
        z=cases{j}.(fields{n});plot(ax,real(z),imag(z),markers{j},'Color',colors(j,:), ...
            'LineStyle','none','LineWidth',1.6,'MarkerSize',9,'DisplayName',caseNames{j});
    end
    xline(ax,0,'k--','HandleVisibility','off');grid(ax,'on');axis(ax,'padded');
    title(ax,titles{n});xlabel(ax,'Real part');ylabel(ax,'Imaginary part');legend(ax,'Location','best');
end
title(tl,'Pourazad-style per-unit VSM-GFM contraction model');
save_figure(fig,outDir,'05_case_eigenvalue_comparison');
fig=figure('Color','w','Position',[80,60,1300,950]);tl=tiledlayout(fig,2+size(intervals,1),1,'TileSpacing','compact');
ax=nexttile(tl);semilogy(ax,gamma0_deg,scan.minEigM,'LineWidth',1.5);grid(ax,'on');ylabel(ax,'lambda min (M)');
title(ax,'Positive reference metric; each gamma0 has a separate M');
xline(ax,rad2deg(nominal),'m-.','Nominal 1 pu');
ax=nexttile(tl);overview=scan;
for field={'lambdaPaper','lambdaDirect','lambdaCoordinate','lambdaMoving'}
    z=overview.(field{1});overview.(field{1})=sign(z).*log10(1+abs(z));
end
plot_consistency(ax,gamma0_deg,overview);xlim(ax,[0,120]);
shade_intervals(ax,rad2deg(intervals));xline(ax,rad2deg(nominal),'m-.','Nominal 1 pu','HandleVisibility','off');
ylabel(ax,'sign(lambda) log10(1+|lambda|)');title(ax,'Signed-log overview: zero and eigenvalue signs are preserved');
for n=1:size(intervals,1)
    ax=nexttile(tl);ii=min(n,size(intervals,1));bd=rad2deg(intervals(ii,:));pad=.3*diff(bd);
    gz=linspace(bd(1)-pad,bd(2)+pad,2001);sz=scan_model(deg2rad(gz),par);
    plot_consistency(ax,gz,sz);legend(ax,'off');
    ymin=min(sz.lambdaPaper);ymax=max(sz.lambdaPaper);yp=.12*(ymax-ymin);
    ylim(ax,[ymin-yp,ymax+yp]);xlim(ax,[gz(1),gz(end)]);
    shade_intervals(ax,rad2deg(intervals));
    if nominal>=intervals(n,1) && nominal<=intervals(n,2),xline(ax,rad2deg(nominal),'m-.','Nominal 1 pu');end
    for boundary=bd,xline(ax,boundary,'k:','HandleVisibility','off');end
    title(ax,sprintf('Component %d: shaded interval satisfies the full paper-metric test',n));
    xlabel(ax,'Equilibrium gamma0 (deg)');
end
title(tl,'Pourazad-style per-unit VSM-GFM contraction model: retuned candidate');
save_figure(fig,outDir,'02_contraction_interval');
% Required full Sylvester boundary diagnostic: no diagonal-only certificate.
fig=figure('Color','w','Position',[70,60,1350,1000]);tl=tiledlayout(fig,3,2,'TileSpacing','compact');
fields={'s11','s22','detSymF'};names={'s11','s22','det[Sym(F)]'};
zoomBounds=[min(rad2deg(intervals(:)))-.08,max(rad2deg(intervals(:)))+.08];
gz=linspace(zoomBounds(1),zoomBounds(2),3001);sz=scan_model(deg2rad(gz),par);
for n=1:3
    ax=nexttile(tl);z=scan.(fields{n});plot(ax,gamma0_deg,sign(z).*log10(1+abs(z)),'LineWidth',1.3);
    yline(ax,0,'k--');grid(ax,'on');xlabel(ax,'gamma0 (deg)');ylabel(ax,['signed log: ',names{n}]);
    title(ax,[names{n},': full scan, sign-preserving scale']);
    for g=rad2deg(intervals(:)).',xline(ax,g,'k:');end
    shade_intervals(ax,rad2deg(intervals));
    ax=nexttile(tl);plot(ax,gz,sz.(fields{n}),'b','LineWidth',1.5);yline(ax,0,'k--');grid(ax,'on');
    for g=rad2deg(intervals(:)).',xline(ax,g,'k:');end
    shade_intervals(ax,rad2deg(intervals));xlabel(ax,'gamma0 (deg)');ylabel(ax,names{n});
    title(ax,[names{n},': boundary detail (linear scale)']);
end
title(tl,'PU contraction boundary causes: both endpoints are tested, not assumed');
save_figure(fig,outDir,'09_boundary_causes');
fig=figure('Color','w','Position',[90,80,1150,650]);ax=axes(fig);hold(ax,'on');
low=scan.Pstar<=1;plo=scan.Pstar;plo(~low)=NaN;phi=scan.Pstar;phi(low)=NaN;
plot(ax,gamma0_deg,plo,'Color',[0,.5,.25],'LineWidth',2,'DisplayName','0 <= Pstar <= 1 pu');
plot(ax,gamma0_deg,phi,'Color',[.85,.4,.05],'LineWidth',2,'DisplayName','Pstar > 1 pu');
yline(ax,1,'k--','Pstar = 1 pu','HandleVisibility','off');
for g=intervals(:).',xline(ax,rad2deg(g),'b:','HandleVisibility','off');end
for j=find(available)
    s=cases{j};plot(ax,rad2deg(s.gamma0),s.Pstar,markers{j},'Color',colors(j,:), ...
        'LineWidth',1.8,'MarkerSize',10,'DisplayName',caseNames{j});
end
grid(ax,'on');xlabel(ax,'Equilibrium gamma0 (deg)');ylabel(ax,'Pstar = Pe(gamma0) (pu)');
title(ax,'Retuned VSM candidate: power and certified intervals');legend(ax,'Location','best');
shade_intervals(ax,rad2deg(intervals));xline(ax,rad2deg(nominal),'m-.','Nominal 1 pu','HandleVisibility','off');
save_figure(fig,outDir,'01_power_angle');
% Region for a FIXED chosen Case-I metric versus instantaneous gamma.
fig=figure('Color','w','Position',[100,80,1200,650]);tl=tiledlayout(fig,1,2);
ax=nexttile(tl);hold(ax,'on');
for j=find(available)
    s=cases{j};lam=frozen_lambda(gamma0,s,par);plot(ax,gamma0_deg,lam,'Color',colors(j,:), ...
        'LineWidth',1.4,'DisplayName',sprintf('%s, gamma0=%.4f deg',caseNames{j},rad2deg(s.gamma0)));
end
yline(ax,0,'k--','HandleVisibility','off');grid(ax,'on');xlabel(ax,'Instantaneous gamma (deg)');
ylabel(ax,'lambda max, fixed reference metric');legend(ax,'Location','best');title(ax,'gamma is distinct from gamma0');
ax=nexttile(tl);hold(ax,'on');
for j=[2,4]
    q=cases{j};gi=linspace(q.gamma0-deg2rad(2),q.gamma0+deg2rad(2),1001);
    plot(ax,rad2deg(gi),frozen_lambda(gi,q,par),'Color',colors(j,:),'DisplayName',caseNames{j});
end
yline(ax,0,'k--');grid(ax,'on');legend(ax,'Location','best');
xlabel(ax,'Instantaneous gamma (deg)');ylabel(ax,'lambda max (Sym F)');
title(ax,'Certified references: nu=1 slice, not a trajectory certificate');
save_figure(fig,outDir,'06_sg_operating_and_state_regions');
fig=figure('Color','w','Position',[20,60,2400,900]);ids=find(available);nc=numel(ids);
tl=tiledlayout(fig,2,nc,'TileSpacing','compact');
for col=1:nc
    j=ids(col);s=cases{j};tr=trajectories{j};ax=nexttile(tl,col);hold(ax,'on');
    show=true(size(tr.time));
    if j==5,show=tr.mainDisplayMask;end % Identical indices in phase and distance panels.
    phaseAngle=tr.gammaUnwrapped(show);phaseState=tr.state(show,:);
    shownTime=tr.time(show);shownDistance=tr.distance(show);
    plot(ax,rad2deg(phaseAngle),phaseState(:,2)-1,'b','LineWidth',1.6);
    plot(ax,rad2deg(s.gamma0),0,'ro','MarkerFaceColor','r');
    plot(ax,rad2deg(tr.state(1,1)),tr.state(1,2)-1,'ks');grid(ax,'on');axis(ax,'padded');
    xlabel(ax,'gamma (deg)');ylabel(ax,'frequency deviation (pu)');title(ax,[caseNames{j},': phase plane']);
    plot(ax,rad2deg(phaseAngle(end)),phaseState(end,2)-1,'gd','MarkerSize',7,'LineWidth',1.5);
    if true
        xlabel(ax,'gamma unwrapped (deg)');
        text(ax,rad2deg(s.gamma0),0,'  Original equilibrium','Color','r','VerticalAlignment','top');
    end
    subtitle(ax,sprintf('gamma0=%.5f deg, Pstar=%.5f pu',rad2deg(s.gamma0),s.Pstar));
    if j==5,subtitle(ax,sprintf('gamma0=120 deg; early time 0 to %.3g',shownTime(end)));end
    if j==5
        text(ax,rad2deg(phaseAngle(end)),phaseState(end,2)-1,'Displayed endpoint', ...
            'HorizontalAlignment','right','VerticalAlignment','bottom','Color',[0,.5,0],'FontSize',9);
    end
    ax.YAxis.Exponent=0;ytickformat(ax,'%.4f');
    ax=nexttile(tl,col+nc);plot(ax,shownTime,shownDistance,'b','LineWidth',1.6,'DisplayName','Metric distance');hold(ax,'on');
    if j==5
        plot(ax,shownTime(end),shownDistance(end),'gd','MarkerSize',7,'LineWidth',1.5);
        text(ax,shownTime(end),shownDistance(end),'Displayed endpoint','HorizontalAlignment','right', ...
            'VerticalAlignment','bottom','Color',[0,.5,0],'FontSize',9);
    end
    if tr.envelopeCertified
        plot(ax,shownTime,tr.envelope(show),'r--','LineWidth',1.4,'DisplayName','Certified regional envelope');legend(ax,'Location','best');
    end
    grid(ax,'on');xlabel(ax,'Normalized analytical time');ylabel(ax,'L(t) = ||Theta_ref (x-x0)||');title(ax,[caseNames{j},': metric distance']);
    ax.YAxis.Exponent=0;ytickformat(ax,'%.4f');
end
title(tl,'Five representative operating points: common +1 deg test perturbation; 120 deg: early-time segment');
save_figure(fig,outDir,'07_case_trajectories_figure4');
% Supplementary full-time data are never truncated or wrapped.
s=cases{5};tr=trajectories{5};gd=rad2deg(tr.gammaUnwrapped);q=tr.state(:,2)-1;
fig=figure('Color','w','Position',[70,80,1550,600]);tl=tiledlayout(fig,1,3,'TileSpacing','compact');
ax=nexttile(tl);plot(ax,gd,q,'b','LineWidth',1.5);hold(ax,'on');
plot(ax,rad2deg(s.gamma0),0,'ro','MarkerSize',10,'LineWidth',1.5,'DisplayName','Original equilibrium');
plot(ax,gd(1),q(1),'ks','MarkerSize',8,'DisplayName','Initial');
plot(ax,gd(end),q(end),'gd','MarkerSize',6,'DisplayName','Final');
xlabel(ax,'gamma unwrapped (deg)');ylabel(ax,'frequency deviation (pu)');title(ax,'Full-time phase plane');grid(ax,'on');
ax=nexttile(tl);plot(ax,tr.time,gd,'b','LineWidth',1.5);hold(ax,'on');
yline(ax,rad2deg(s.gamma0),'r--','Original equilibrium');
plot(ax,tr.time(1),gd(1),'ks');plot(ax,tr.time(end),gd(end),'gd');
xlabel(ax,'Normalized analytical time');ylabel(ax,'gamma unwrapped (deg)');title(ax,'Full unwrapped angle');grid(ax,'on');
ax=nexttile(tl);plot(ax,tr.time,tr.distance,'b','LineWidth',1.5);hold(ax,'on');
yline(ax,0,'r--');plot(ax,tr.time(1),tr.distance(1),'ks');plot(ax,tr.time(end),tr.distance(end),'gd');
xlabel(ax,'Normalized analytical time');ylabel(ax,'L(t), original reference');title(ax,'Full metric distance');grid(ax,'on');
title(tl,{'Supplementary: 120 deg, full simulation, common +1 deg test perturbation', ...
    'Red = original equilibrium, black square = initial, green diamond = final (markers may coincide)'});
save_figure(fig,outDir,'10_above_R2_full_time_diagnostic');


%% 7. Save reproducible numerical outputs and report
results=struct('status','RETUNED_CURRENT_LIMITER_FREE_VSM_GFM_COMPLETE','parameters',par,'physicalBase',struct('f_base',f_base,'omega_base',omega_base),'sgValidation',sgValidation, ...
    'scan',scan,'gamma0_rad',gamma0,'gamma0_deg',gamma0_deg,'paperIntervals_rad',intervals, ...
    'paperIntervals_deg',rad2deg(intervals),'directIntervals_deg',rad2deg(intervalsDirect), ...
    'coordinateIntervals_deg',intervalsCoordinate,'movingIntervals_deg',intervalsMoving, ...
    'caseAvailable',available,'cases',{cases},'caseTable',T,'trajectories',{trajectories},'validation',validation, ...
    'pointResults',pointResults,'pointModels',{pointModels},'representativeIndices',representativeIndices,'classificationTolerance',tol, ...
    'boundaryCauses',boundaryTable,'nominal',snom,'nominalCertified',nominalCertified, ...
    'nominalStable',nominalStable,'tuningComparison',comparison,'intervalTable',intervalTable);
results.method=['Paper-literal PU reference metric (c_paper=omega0_pu*beta_gamma, omega0_pu=1). ' ...
    'Full direct definition uses same frozen metric; Method B and moving-angle audits saved separately.'];
save(fullfile(outDir,'baseline_gfm_results.mat'),'results','-v7.3');
writetable(T,fullfile(outDir,'gfm_sg_case_comparison.csv'));
writetable(T,fullfile(outDir,'baseline_gfm_candidates.csv'));
write_report(outDir,results);
fprintf('\nFINAL CASE TABLE\n');disp(T);
fprintf('All selected equilibrium, derivative, matrix, trajectory and envelope checks passed.\n');
fprintf('====================================================\nRETUNED CURRENT-LIMITER-FREE VSM-GFM SUMMARY\n====================================================\n');
fprintf('Candidate contraction-oriented retuning\nH=%.1f s; K=%.0f; Kd=%.0f; D_eff=%.0f\n',H,K,Kd,D);
fprintf('Nominal Pstar=1 pu; gamma0=%.12f deg\n',rad2deg(nominal));disp(intervalTable);
fprintf('Nominal contraction certified: %s\nlambda_max(SymF): %.12g\nLocal stability at nominal: %s\n',yesno{1+nominalCertified},max(snom.eigPaper),yesno{1+nominalStable});
for j=1:5
    if available(j),fprintf('%s = %.12f deg\n',caseNames{j},rad2deg(selected(j)));
    else,fprintf('%s: category unavailable\n',caseNames{j});end
end
fprintf('\n========================================================\nFIVE REPRESENTATIVE NONLINEAR TRAJECTORIES\n========================================================\n');
for j=1:5,fprintf('%s: gamma0=%.9f deg; classification=%s\n',caseNames{j},rad2deg(selected(j)),trajectories{j}.finalDiagnostics.interpretation);end
fprintf('Same MATLAB file updated; outputs: %s\n',outDir);
%% Separate voltage-sag study: baseline variables/results above are unchanged.
voltageSagResults=voltage_sag_sweep(par,outDir); %#ok<NASGU>
%% Separate coarse duration sweep; previous analyses and exports are unchanged.
voltageSagDurationResults=voltage_sag_duration_sweep(par,outDir); %#ok<NASGU>
%% Figure-5-style post-processing; no prior results are modified.
figure5_gfm_postprocess(results,voltageSagResults,outDir);
%% Additive H-D_eff sweep; all earlier analysis definitions are unchanged.
hdSweepResults=hd_parameter_sweep(outDir);
hd_boundary_sanity(hdSweepResults,outDir);
diary off;

%% Local functions -- all implementation remains in THIS MATLAB file
function pe=power_smib(g,p)
pe=p.E*p.U1/p.Xeq*sin(g);
end
function dx=rhs(~,x,p,ps)
nu=x(2);assert(nu>0.5,'Normal-frequency domain nu>0.5 exceeded; no floor added.');
dx=[nu-1;((ps-power_smib(x(1),p))/nu-p.D*(nu-1))/p.mass];
end
function [A,B]=jacobians(g,p,nu,ps)
if nargin<3,nu=1;end
if nargin<4,ps=power_smib(g,p);end
C=p.E/(p.mass*p.Xeq);
A=[0,1;-C*p.U1*cos(g)/nu,(-(ps-power_smib(g,p))/nu^2-p.D)/p.mass];
B=[0;-C*sin(g)/nu];
end
function [beta,bd,bdd]=actual_beta(x,p,ps)
C=p.E/(p.mass*p.Xeq);g=x(1);nu=x(2);
f=rhs(0,x,p,ps);[A,~]=jacobians(g,p,nu,ps);
beta=-C*sin(g)/nu;
grad=[-C*cos(g)/nu,C*sin(g)/nu^2];
Hess=[C*sin(g)/nu,C*cos(g)/nu^2;C*cos(g)/nu^2,-2*C*sin(g)/nu^3];
bd=grad*f;bdd=f'*Hess*f+grad*A*f;
end
function S=paper_sym(a1,a2,r)
S=[-r,(1+a1-a2*r-r^2)/2;(1+a1-a2*r-r^2)/2,a2+r];
end
function s=point_model(g,p)
C=p.E/(p.mass*p.Xeq);ps=power_smib(g,p);flow=rhs(0,[g;1],p,ps);[A,B]=jacobians(g,p);
a1=A(1,1);a2=A(1,2);a3=A(2,1);a4=A(2,2);b1=B(1);b2=B(2);
den=a2*b2-a4*b1;assert(abs(den)>1e-12,'Singular beta: do not regularize.');
T=[a2,a2*b1/den;a2*(a3*b1-a1*b2)/den,a2*b2/den];
alpha1=a2*a3-a1*a4;alpha2=a1+a4;beta=den/a2;
assert(norm(T\A*T-[0,1;alpha1,alpha2],'fro')<1e-9 && norm(T\B-[0;beta])<1e-9);
c=-C*p.omega0_pu*cos(g);Theta=[-beta,0;-c,-beta];M=Theta'*Theta;
% Full derivative of the SAME reference-designed metric. Along the fixed
% reference it is zero. Even extending its coefficients as functions of g,
% DTheta(g)*f0 is zero because flow(1)=0. No beta_dot convention is hidden.
ThetaGamma=[C*cos(g),0;-C*p.omega0_pu*sin(g),C*cos(g)];
ThetaDotReference=ThetaGamma*flow(1);
Fdirect=(ThetaDotReference+Theta*A)/Theta;Sdirect=(Fdirect+Fdirect')/2;
Spaper=paper_sym(alpha1,alpha2,c/beta);
% Method B includes BOTH beta_gamma*gamma_dot and beta_w*w_dot.
[~,bd,bdd]=actual_beta([g;1],p,ps);
ThetaB=[-beta,0;-bd,-beta];ThetaDotB=[-bd,0;-bdd,-bd];
FB=(ThetaDotB+ThetaB*A)/ThetaB;SB=(FB+FB')/2;
% Additional audit: differentiating the paper coefficient along a path
% gamma_dot=omega0_pu gives nonzero Theta_dot; this path is NOT the equilibrium.
ThetaDotMoving=ThetaGamma*p.omega0_pu;Fm=(ThetaDotMoving+Theta*A)/Theta;Sm=(Fm+Fm')/2;
s=struct('gamma0',g,'Pstar',ps,'A',A,'B',B,'T',T,'alpha1',alpha1,'alpha2',alpha2, ...
    'beta',beta,'betaDotPaper',c,'betaDotCoordinate',bd,'Theta',Theta,'M',M, ...
    'ThetaDotReference',ThetaDotReference,'ThetaCoordinate',ThetaB,'ThetaDotCoordinate',ThetaDotB, ...
    'ThetaDotLiteralMoving',ThetaDotMoving,'Fdirect',Fdirect,'SymF_paper',Spaper, ...
    'SymF_direct',Sdirect,'SymF_coordinate',SB,'SymF_literal_moving',Sm, ...
    'eigA',eig(A),'eigM',sort(eig(M)),'eigPaper',sort(eig(Spaper)), ...
    'eigDirect',sort(eig(Sdirect)),'eigCoordinate',sort(eig(SB)),'eigMoving',sort(eig(Sm)), ...
    'eigRepresentative',sort(eig(paper_sym(-C*p.U1,alpha2-ps/p.mass,c/beta))), ...
    'equilibriumResidual',norm(flow,inf));
assert(min(s.eigM)>0 && det(M)>0 && norm(Spaper-Sdirect,'fro')<1e-7);
end
function s=scan_model(g,p)
C=p.E/(p.mass*p.Xeq);b=p.D/p.mass;beta=-C*sin(g);c=-C*p.omega0_pu*cos(g);r=c./beta;
a=-C*p.U1*cos(g);off=(1+a+b*r-r.^2)/2;
lmax=-b/2+hypot(b/2-r,off);lmin=-b/2-hypot(b/2-r,off);
% Stable metric eigenvalue computation avoids cancellation at ill-conditioned angles.
M11=beta.^2+c.^2;M12=beta.*c;M22=beta.^2;metricMax=(M11+M22+hypot(M11-M22,2*M12))/2;
metricMin=beta.^4./metricMax;
lamCoordinate=-b/2+hypot(b/2,(1+a)/2);
% Full moving-angle F11=0, F22=-b+2r,
% F21=a+b*r-2*r^2-omega0^2 (r_dot=-omega0^2-r^2).
f21=a+b*r-2*r.^2-p.omega0_pu^2;d=-b+2*r;
lamMoving=d/2+hypot(d/2,(1+f21)/2);
ps=power_smib(g,p);res=abs((ps-power_smib(g,p))/p.mass);
% Direct frozen matrix multiplication, independently vectorized by pages.
N=numel(g);A=zeros(2,2,N);A(1,2,:)=1;A(2,1,:)=a;A(2,2,:)=-b;
Theta=zeros(2,2,N);Theta(1,1,:)=-beta;Theta(2,1,:)=-c;Theta(2,2,:)=-beta;
ThetaInv=zeros(2,2,N);ThetaInv(1,1,:)=-1./beta;ThetaInv(2,1,:)=c./beta.^2;ThetaInv(2,2,:)=-1./beta;
Fd=pagemtimes(pagemtimes(Theta,A),ThetaInv);Sd=(Fd+permute(Fd,[2,1,3]))/2;
z11=reshape(Sd(1,1,:),1,[]);z12=reshape(Sd(1,2,:),1,[]);z22=reshape(Sd(2,2,:),1,[]);
ld=(z11+z22)/2+hypot((z11-z22)/2,z12);
residual=sqrt((z11+r).^2+2*(z12-off).^2+(z22+b-r).^2);
% Use relative residual across extremely ill-conditioned mathematical endpoints.
relativeResidual=residual./max(1,sqrt(r.^2+2*off.^2+(-b+r).^2));
M=zeros(2,2,N);M(1,1,:)=M11;M(1,2,:)=M12;M(2,1,:)=M12;M(2,2,:)=M22;
B=zeros(2,1,N);B(2,1,:)=beta;
s=struct('Pstar',ps,'equilibriumResidual',res,'A',A,'B',B,'Theta',Theta,'M',M, ...
    'alpha1',a,'alpha2',-b,'beta',beta,'betaDotPaper',c,'betaDotCoordinate',zeros(size(g)), ...
    'lambdaPaper',lmax,'lambdaPaperMin',lmin,'lambdaDirect',ld,'lambdaCoordinate',lamCoordinate, ...
    'lambdaMoving',lamMoving,'minEigM',metricMin,'maxEigM',metricMax,'detM',beta.^4, ...
    'maxRealEigA',real((-b+sqrt(complex(b^2+4*a)))/2),'directResidual',relativeResidual, ...
    'paperCertified',metricMin>0 & lmax<0,'aboveOnePu',ps>1, ...
    's11',-r,'s22',-b+r,'s12',off,'detSymF',r.*(b-r)-off.^2);
end
function tab=boundary_causes(I,p)
tab=table();
for k=1:size(I,1)
    for edge=1:2
        s=point_model(I(k,edge),p);S=s.SymF_paper;
        sc=max(1,norm(S,'fro'));dv=det(S);cause="unresolved";
        if min(s.eigM)<=0,cause="lambda_min(M)=0";
        elseif abs(S(1,1))<1e-10 && abs(dv)<1e-8,cause="s11=0 and det=0";
        elseif abs(S(2,2))<1e-10 && abs(dv)<1e-8,cause="s22=0 and det=0";
        elseif abs(dv)/sc^2<1e-10,cause="det(SymF)=0";
        end
        assert(cause~="unresolved" && abs(max(s.eigPaper))<1e-8);
        row=table(k,edge,rad2deg(I(k,edge)),S(1,1),S(2,2),dv,min(s.eigM),max(s.eigPaper),cause, ...
            'VariableNames',{'Interval','Edge','gamma0_deg','s11','s22','detSymF','minEigM','lambdaMax','Cause'});
        tab=[tab;row]; %#ok<AGROW>
    end
end
end
function shade_intervals(ax,I)
hold(ax,'on');yl=ylim(ax);
for k=1:size(I,1)
    q=patch(ax,[I(k,1),I(k,2),I(k,2),I(k,1)],[yl(1),yl(1),yl(2),yl(2)], ...
        [.6,.85,.65],'FaceAlpha',.18,'EdgeColor','none','HandleVisibility','off');
    uistack(q,'bottom');
end
ylim(ax,yl);
end
function z=lambda_at(g,p,method)
s=point_model(g,p);
if strcmp(method,'paper'),z=max(s.eigPaper);else,z=max(s.eigDirect);end
end
function intervals=refine_intervals(g,y,f)
idx=find(y(1:end-1).*y(2:end)<0);roots=zeros(size(idx));
for k=1:numel(idx),roots(k)=fzero(f,g(idx(k):idx(k)+1),optimset('TolX',1e-13));end
edges=[g(1),roots,g(end)];intervals=zeros(0,2);
for k=1:numel(edges)-1
    if f(mean(edges(k:k+1)))<0,intervals(end+1,:)=edges(k:k+1);end %#ok<AGROW>
end
end
function print_intervals(label,I)
fprintf('\n%s:\n',label);
if isempty(I),fprintf('EMPTY\n');else
    for k=1:size(I,1),fprintf('gamma0_min=%.12f deg; gamma0_max=%.12f deg\n',rad2deg(I(k,:)));end
end
end
function z=frozen_lambda(g,s,p)
% nu=1 slice of the FULL Jacobian, with Pstar held fixed at the reference.
r=s.betaDotPaper/s.beta;a=-p.E*p.U1/(p.mass*p.Xeq)*cos(g);
d=-p.D/p.mass-(s.Pstar-power_smib(g,p))/p.mass;
z=d/2+hypot(-r-d/2,(1+a-d*r-r^2)/2);
end
function [rate,R]=rectangle_certificate(s,p,dg,dw)
% Bound both variable Jacobian entries over a convex (gamma,w) rectangle.
% A21 interval from cos(g)/nu; A22 interval from (Pstar-Pe)/nu^2.
% lambda_max(Sym(Theta*A/Theta)) is convex in the two entries, so its
% maximum over their enclosing rectangle is bounded by the four corners.
gs=s.gamma0+[-dg,dg];nus=1+[-dw,dw];
assert(gs(1)>0 && gs(2)<pi/2 && nus(1)>0.5);
C=p.E*p.U1/(p.mass*p.Xeq);aVals=[];dVals=[];
for g=gs
    for nu=nus
        aVals(end+1)=-C*cos(g)/nu; %#ok<AGROW>
        dVals(end+1)=-p.D/p.mass-(s.Pstar-power_smib(g,p))/(p.mass*nu^2); %#ok<AGROW>
    end
end
rate=-inf;
for aa=[min(aVals),max(aVals)]
    for dd=[min(dVals),max(dVals)]
        F=s.Theta*[0,1;aa,dd]/s.Theta;
        rate=max(rate,max(eig((F+F')/2)));
    end
end
Ti=s.Theta\eye(2);R=min(dg/norm(Ti(1,:)),dw/norm(Ti(2,:)));
end
function tr=simulate_case(s,p,j)
% Common +1 degree perturbation is our own visualization/test choice,
% NOT a disturbance magnitude taken from the reference paper.
% Keep xd=[gamma0;1] and Pstar=Pe(gamma0), with identical offsets for all five.
% Figure-5-style disturbance-severity analysis remains a separate future task.
g0=s.gamma0;delta=deg2rad(1);cert=[NaN,NaN];inv=cert;rate=NaN;certified=false;R=NaN;dw=NaN;
if s.contractionCertified
    % Restore the previously verified regional-envelope check unchanged.
    dg=min(deg2rad(15),.9*min(g0,pi/2-g0));dw=.1;
    [rate,R]=rectangle_certificate(s,p,dg,dw);
    certified=rate<0 && norm(s.Theta*[delta;0])<R;
    cert=g0+[-dg,dg];inv=cert;
end
if ~certified,fprintf('%s: no validated envelope for common +1 deg displacement.\n',s.name);end
xinit=[g0+delta;1];
horizon=max(1,15/max(1e-4,abs(max(real(s.eigA)))));
% Resolve both the fast frequency transient and the slow angle evolution.
times=unique(sort([linspace(0,horizon,3001),logspace(-8,log10(min(1,horizon)),600)])).';
opt=odeset('RelTol',1e-9,'AbsTol',1e-11,'MaxStep',horizon/1000);
[t,x]=ode15s(@(tt,xx) rhs(tt,xx,p,s.Pstar),times,xinit,opt);
% Keep distance relative to each representative's ORIGINAL equilibrium.
% L(t) does not tend to zero if the trajectory settles on another branch.
L=vecnorm(s.Theta*(x-[g0,1]).',2,1).';
[~,xcheck]=ode15s(@(tt,xx) rhs(tt,xx,p,s.Pstar),times,xinit, ...
    odeset('RelTol',1e-10,'AbsTol',1e-12,'MaxStep',horizon/2000));
err=max(abs(x-xcheck),[],'all')/max(1,max(abs(xcheck),[],'all'));
nu=x(:,2);assert(all(nu>0.5));
envelope=nan(size(t));
if certified
    envelope=L(1)*exp(rate*t);
    assert(all(L<R) && min(x(:,1))>=cert(1)-1e-10 && max(x(:,1))<=cert(2)+1e-10 && max(abs(x(:,2)-1))<dw);
    assert(all(L<=envelope+1e-7*max(1,L(1))));
end
fprintf('%s trajectory: L0=%.8g, Lend=%.8g, maxL=%.8g, solver discrepancy=%.3e\n',s.name,L(1),L(end),max(L),err);
fprintf('  nu range=[%.10f,%.10f]; actual 1/nu torque conversion retained.\n',min(nu),max(nu));
if certified,fprintf('  Certified rectangle rate=%.12g; gamma [deg]=[%.9f,%.9f]; |nu-1|<=%.9g\n',rate,rad2deg(cert),dw);end
tr=struct('time',t,'state',x,'distance',L,'initialState',xinit,'solverRelativeError',err, ...
    'envelopeCertified',certified,'envelope',envelope,'regionalRate',rate, ...
    'certificateInterval',cert,'invariantInterval',inv,'metricBallRadius',R, ...
    'rectangleSpeedHalfWidth',dw,'nu',nu);
if true
    % The integrated angle is already continuous. Verify sampling before
    % unwrap so that a genuine rotation is never mistaken for a wrap jump.
    assert(all(abs(diff(x(:,1)))<pi),'Angle sampling too coarse for unwrap.');
    gu=unwrap(x(:,1));
    assert(max(abs(gu-x(:,1)))<1e-10,'Unwrap changed the integrated angular branch.');
    finalDeg=rad2deg(gu(end));finalMod=mod(finalDeg,360);
    finalNu=x(end,2);dnu=finalNu-1;dg=finalDeg-rad2deg(g0);
    pe=power_smib(gu(end),p);rp=pe-s.Pstar;
    assert(abs(pe-power_smib(deg2rad(finalMod),p))<1e-10);
    % Explicit absolute numerical thresholds; classification uses final
    % values, not the visual appearance of the phase-plane trajectory.
    angleTolDeg=1e-3;frequencyTol=1e-5;powerTol=1e-4;
    if ~all(isfinite([finalDeg,finalNu,rp,err])) || err>=1e-5
        interpretation='Numerically inconclusive';
    elseif abs(dnu)<=frequencyTol && abs(rp)<=powerTol
        if abs(dg)<=angleTolDeg
            interpretation='Returns to original equilibrium';
        else
            interpretation='Approaches another equilibrium branch';
        end
    else
        interpretation='Still drifting / divergent';
    end
    tr.gammaUnwrapped=gu;
    tr.mainDisplayMask=true(size(t));
    if j==5
        % Stop at first sampled 20-degree departure from the ORIGINAL equilibrium.
        % No simulation data or long-time diagnostics are truncated.
        stop=find(abs(rad2deg(gu-g0))>=20,1,'first');
        if isempty(stop),stop=numel(t);end
        tr.mainDisplayMask=(1:numel(t)).'<=stop;
        tr.displayedEndpointIndex=stop;tr.displayedEndpointTime=t(stop);
        assert(all(tr.mainDisplayMask(1:stop)) && ~any(tr.mainDisplayMask(stop+1:end)));
        assert(L(stop)>L(1),'Expected early-time growth from unstable reference.');
        fprintf('120 deg MAIN display ends at t=%.12g, gamma=%.12g deg, L=%.12g; full simulation ends at t=%.12g.\n',t(stop),rad2deg(gu(stop)),L(stop),t(end));
    end
    tr.finalDiagnostics=struct('originalGamma_deg',rad2deg(g0),'gammaUnwrapped_deg',finalDeg, ...
        'gammaMod360_deg',finalMod,'nu',finalNu,'frequencyDeviation',dnu, ...
        'initialGamma_deg',rad2deg(xinit(1)),'L_final',L(end),'deltaGamma_deg',dg,'Pe_pu',pe,'Pstar_pu',s.Pstar,'powerResidual_pu',rp, ...
        'angleTolerance_deg',angleTolDeg,'frequencyTolerance_pu',frequencyTol, ...
        'powerTolerance_pu',powerTol,'interpretation',interpretation);
    fprintf('\n====================================================\nREPRESENTATIVE FINAL-STATE DIAGNOSTICS\n====================================================\n');
    fprintf('%s\nInitial perturbed gamma = %.12f deg\nL_final = %.12g\n',s.name,rad2deg(xinit(1)),L(end));
    fprintf('Original equilibrium gamma0      = %.12f deg\n',rad2deg(g0));
    fprintf('Final unwrapped gamma            = %.12f deg\n',finalDeg);
    fprintf('Final gamma mod 360              = %.12f deg\n',finalMod);
    fprintf('Final nu                         = %.12g\n',finalNu);
    fprintf('Final frequency deviation nu-1   = %+.12g\n',dnu);
    fprintf('Delta gamma from original eq     = %+.12f deg\n\n',dg);
    fprintf('Pe(gamma_final)                  = %.12g pu\n',pe);
    fprintf('Pstar reference                   = %.12g pu\n',s.Pstar);
    fprintf('Pe(gamma_final)-Pstar residual   = %+.12g pu\n\n',rp);
    fprintf('Interpretation:\n- %s\n',interpretation);
    fprintf('Final-sample tolerances: angle %.3g deg; nu-1 %.3g pu; power %.3g pu\n',angleTolDeg,frequencyTol,powerTol);
    fprintf('====================================================\n');
end
end
function v=validate_model(p,cases,available)
v=struct('jacobian',0,'input',0,'beta',0,'coordinates',0);
C=p.E/(p.mass*p.Xeq);
for g=[.2,.8,1.4]
    x=[g;1.04];ps=1;h=1e-6;[A,B]=jacobians(g,p,x(2),ps);J=zeros(2);
    for k=1:2,e=zeros(2,1);e(k)=h;J(:,k)=(rhs(0,x+e,p,ps)-rhs(0,x-e,p,ps))/(2*h);end
    pp=p;pm=p;pp.U1=p.U1+h;pm.U1=p.U1-h;
    Bfd=(rhs(0,x,pp,ps)-rhs(0,x,pm,ps))/(2*h);
    v.jacobian=max(v.jacobian,norm(J-A,'fro')/norm(A,'fro'));v.input=max(v.input,norm(Bfd-B)/norm(B));
    flow=rhs(0,x,p,ps);
    [~,bd,bdd]=actual_beta(x,p,ps);
    [bp,dp]=actual_beta(x+h*flow,p,ps);[bm,dm]=actual_beta(x-h*flow,p,ps);
    v.beta=max(v.beta,max(abs([(bp-bm)/(2*h)-bd,(dp-dm)/(2*h)-bdd])./max(1,abs([bd,bdd]))));
    nu=x(2);q=nu-1;
    % Exact change between PU absolute-frequency and PU deviation states.
    absolute=[q;((ps-power_smib(g,p))/(1+q)-p.D*q)/p.mass];
    v.coordinates=max(v.coordinates,norm(absolute-flow,inf));
end
assert(v.jacobian<1e-7 && v.input<1e-7 && v.beta<1e-5 && v.coordinates<1e-9);
for j=find(available),assert(cases{j}.equilibriumResidual<1e-10);end
end
function plot_consistency(ax,g,s)
plot(ax,g,s.lambdaPaper,'b-','LineWidth',1.6,'DisplayName','Paper Eq15c / 2');hold(ax,'on');
plot(ax,g,s.lambdaDirect,'r--','LineWidth',1.3,'DisplayName','Direct: same reference metric');
plot(ax,g,s.lambdaCoordinate,'Color',[.1,.55,.2],'LineWidth',1,'DisplayName','Method B: actual-derivative metric');
plot(ax,g,s.lambdaMoving,'Color',[.6,.2,.65],'LineStyle',':','LineWidth',1,'DisplayName','Literal moving-angle audit');
yline(ax,0,'k--','HandleVisibility','off');grid(ax,'on');ylabel(ax,'lambda max (Sym F)');legend(ax,'Location','best');
end
function save_figure(fig,outDir,stem)
print(fig,fullfile(outDir,[stem,'.png']),'-dpng','-r180');savefig(fig,fullfile(outDir,[stem,'.fig']));
end
function v=legacy_sg_verification(outDir)
w=1;H=3.5;D=8;E=1.3;U=1;X=.3;p=w*E/(2*H*X);a2=-w*D/(2*H);
gd=[44.1,55.5,67,78.5,26.4,35.2,101.5,112.9];labels=[1,1,1,1,2,2,3,3];
ea=complex(zeros(2,8));em=zeros(2,8);er=zeros(2,8);ee=zeros(2,8);res=zeros(1,8);
for i=1:8
    g=deg2rad(gd(i));A=[0,1;-p*U*cos(g),a2];beta=-p*sin(g);bd=-p*w*cos(g);
    Theta=[-beta,0;-bd,-beta];M=Theta'*Theta;r=bd/beta;
    Sr=paper_sym(-p*U,a2,r);Se=paper_sym(A(2,1),a2,r);
    F=Theta*[0,1;-p*U,a2]/Theta;res(i)=norm(Sr-(F+F')/2,'fro');
    ea(:,i)=eig(A);em(:,i)=eig(M);er(:,i)=eig(Sr);ee(:,i)=eig(Se);
end
expected=[-.048537378,-.205570649,-.198325567,-.105355182,1.027017122,.274743914,.205298469,.432245894];
classes=all(max(er(:,labels==1),[],1)<0) && all(max(er(:,labels~=1),[],1)>0) && ...
    all(max(real(ea(:,labels<3)),[],1)<0) && all(max(real(ea(:,labels==3)),[],1)>0);
lo=atand(2*H/D);pass=classes && abs(lo-41.19)<.01 && max(res)<1e-12 && ...
    max(abs(max(er,[],1)-expected))<1e-8;
fprintf('SG LEGACY-CODE REPRODUCTION: gamma_rep=0 deg, gamma0 separate\n');
fprintf('Diagonal bounds %.9f < gamma0 < 90 deg; category reproduction PASS=%d\n',lo,pass);
fprintf('Scaling: legacy Eq15c eigenvalues = 2 * reported half-symmetric eigenvalues.\n');
tab=table(gd.',labels.',max(er,[],1).',max(ee,[],1).', ...
    'VariableNames',{'gamma0_deg','PaperCase','lambdaRepresentativeHalf','lambdaEquilibriumHalf'});disp(tab);
v=struct('pass',pass,'gammaRepresentative_deg',0,'diagonalBounds_deg',[lo,90], ...
    'table',tab,'maxDirectResidual',max(res));
writetable(tab,fullfile(outDir,'sg_validation_cases.csv'));save(fullfile(outDir,'sg_validation_results.mat'),'v');
fig=figure('Color','w','Position',[60,80,1450,550]);tl=tiledlayout(fig,1,3);
colors=[0,.35,.9;0,.55,.2;.85,.1,.1];markers={'o','^','x'};values={ea,em,er};
for n=1:3
    ax=nexttile(tl);hold(ax,'on');
    for j=1:3,z=values{n}(:,labels==j);plot(ax,real(z(:)),imag(z(:)),markers{j}, ...
        'Color',colors(j,:),'MarkerSize',8,'LineWidth',1.5,'DisplayName',sprintf('Case %d',j));end
    xline(ax,0,'k--','HandleVisibility','off');grid(ax,'on');axis(ax,'padded');
    names={'Eigenvalues of A(gamma0)','Eigenvalues of M(gamma0)','Eigenvalues of Sym(F), gamma=0'};
    title(ax,names{n});xlabel(ax,'Real part');ylabel(ax,'Imaginary part');legend(ax,'Location','best');
end
title(tl,'SG legacy verification PASS: representative angle differs from equilibrium angle');
save_figure(fig,outDir,'08_sg_validation_eigenvalues');
end
function write_report(outDir,r)
fid=fopen(fullfile(outDir,'SONUCLAR.md'),'w','n','UTF-8');assert(fid>0);cl=onCleanup(@() fclose(fid));
fprintf(fid,'# Pourazad-style per-unit VSM-GFM contraction model\n\n');
fprintf(fid,'Aynı MATLAB dosyası güncellendi. H=3.4, K=20, Kd=5, D=25; E=U1=1 ve Xeq=0.3. Ana nominal hız omega0_pu=1. Fiziksel taban 50 Hz ve omega_base=314.159265 rad/s yalnızca referans bilgisidir.\n\n');
fprintf(fid,'## Model ve normalizasyon\n\n');
fprintf(fid,'x=[gamma;nu]; gamma_dot=nu-1; nu_dot=[(Pstar-Pe)/nu-D(nu-1)]/(2H); Pe=(E U1/Xeq)sin(gamma). 1/nu korunur. Her gamma0 için nu0=1 ve Pstar=Pe(gamma0).\n\n');
fprintf(fid,'Bu, istenen makale ölçeğindeki analitik modeldir. Fiziksel EMT modelinde gamma_dot=omega_base(nu-1) olur. Yalnızca açı denkleminden omega_base kaldırılırken frekans denklemi aynı tutulduğu için bu değişiklik fiziksel modelin salt koordinat veya ortak zaman dönüşümü değildir. Grafik zamanı fiziksel saniye olarak yorumlanmaz.\n\n');
fprintf(fid,'A21=-Pe_gamma/(2H nu); A22=[-(Pstar-Pe)/nu^2-D]/(2H); B2=-E sin(gamma)/(2H Xeq nu). Dengede alpha2=-25/6.8=-3.676470588235. beta_dot_paper=-E cos(gamma0)/(2H Xeq); beta_dot_paper/beta=cot(gamma0).\n\n');
fprintf(fid,'Ana metrik bu makale katsayısıyla her denge için ayrı tasarlanır ve referans boyunca sabittir. Tam doğrudan F hesabıyla Eq15c/2 uyuşur. Gerçek durum türevi beta_gamma*(nu-1)+beta_nu*nu_dot ayrıca hesaplanır; dengede sıfırdır ve onun farklı metriği sertifika vermez.\n\n');
fprintf(fid,'Candidate contraction-oriented retuning; optimallik iddiası yoktur. D/(2H)=%.12g; D/H=%.12g.\n\n',r.parameters.D/r.parameters.mass,r.parameters.D/r.parameters.H);
fprintf(fid,'## Yeni sonuçlar\n\n');
for k=1:size(r.paperIntervals_deg,1)
    fprintf(fid,'- %.10f° < gamma0 < %.10f°; sınır güçleri %.8f, %.8f pu.\n',r.paperIntervals_deg(k,:),power_smib(r.paperIntervals_rad(k,:),r.parameters));
end
fprintf(fid,'\nHer sınırın sebebi boundary_causes.csv dosyasında; negatif köşegenler tek başına kullanılmadı, determinant ve tam özdeğer koşulu kontrol edildi.\n\n');
fprintf(fid,'| Case | gamma0 (deg) | Pstar (pu) | lambda max paper |\n|---|---:|---:|---:|\n');
for j=find(r.caseAvailable)
    z=r.cases{j};fprintf(fid,'| %s | %.9f | %.9f | %.9f |\n',z.name,rad2deg(z.gamma0),z.Pstar,max(z.eigPaper));
end
fprintf(fid,'\nNominal 1 pu güç açısı %.9f derecedir. Güçler ve 1 pu sınırı ayrıca grafiklenmiştir.\n\n',asind(.3));
fprintf(fid,'Nominal 1 pu sertifikası: %d; yerel kararlılık: %d; lambda_max(SymF)=%.12g.\n\n',r.nominalCertified,r.nominalStable,max(r.nominal.eigPaper));
fprintf(fid,'Önceden doğrulanmış eski PU ayarı: H=0.1, K=20, Kd=40, D=60; aralık [89.673050904,89.874923488] derece; nominal sertifikasız. Bu bilgi yalnızca karşılaştırmadır. Toplam aralık genişliği değişimi %+.9f derece; nominal lambda_max değişimi %+.9f.\n\n',r.tuningComparison.widthChange_deg,r.tuningComparison.nominalLambdaChange);
fprintf(fid,'## Doğrulama ve yörüngeler\n\n');
fprintf(fid,'Ana model Jacobianı, giriş Jacobianı, tam beta türevleri, PU mutlak/sapma koordinat eşdeğerliği ve tüm denge kalıntıları test edildi. Sertifika eşdeğerliği (s11<0, s22<0, det>0) ile lambda_max<0 arasında tüm taramada kontrol edildi.\n\n');
fprintf(fid,'12 farklı denge noktası matrislerden sınıflandırıldı; tüm matrisler MAT dosyalarında, sınıflandırma 12_point_gfm_classification.csv dosyasında saklanır. Beş temsilci gamma(0)=gamma0+1 derece ve nu(0)=1 ile başlatılır. Bu ortak test bozuntusu bize ait görselleştirme tercihidir; referans makaleden alınan bir büyüklük değildir. Figure-5 tarzı bozuntu şiddeti analizi daha sonra yapılacaktır. Ana figürde 120 derece için özgün dengeden ilk 20 derece ayrılmaya kadar aynı zaman indeksleriyle faz düzlemi ve L(t), ek 10_above_R2_full_time_diagnostic figüründe tam simülasyon gösterilir. L(t) özgün dengeye göredir.\n\n');
for j=find(r.caseAvailable),fprintf(fid,'- %s: %s; final gamma=%.9f derece; güç kalıntısı=%+.6g pu.\n',r.cases{j}.name,r.trajectories{j}.finalDiagnostics.interpretation,r.trajectories{j}.finalDiagnostics.gammaUnwrapped_deg,r.trajectories{j}.finalDiagnostics.powerResidual_pu);end
fprintf(fid,'Eski SI sonuçları outputs/diagnostics altında OLD SI-coordinate / mixed-normalization diagnostic; NOT used as thesis baseline result etiketiyle korunmuştur. Ana çıktılar yalnızca yeni PU hesabına aittir.\n');
end


function sweep=voltage_sag_sweep(p,outDir)
% Separate study. No angle jump: each run starts at its PRE-FAULT equilibrium.
% Fault at t=0; clearing at t=1 in NORMALIZED analytical time (not seconds).
% Nominal angle uses asin(0.3) to avoid the rounded-angle/Pstar=1 mismatch.
% Original 12-point inputs and every baseline matrix/figure remain unchanged.
gammas=[asin(p.Pstar_nominal/(p.E*p.U1/p.Xeq)),deg2rad(65)];
Usag_values=[.9,.8,.7,.6,.5,.4,.3,.2,.1];Tfault=1;
nSag=numel(Usag_values);
tol=struct('angle_deg',.1,'frequency_pu',1e-4,'power_pu',1e-4,'matrix',1e-9);
fprintf('\n========================================================\nSEPARATE VOLTAGE-SAG SWEEP: CERTIFIED REFERENCES\n========================================================\n');
fprintf('No angle jump; initial x=[gamma0;1]; fixed pre-fault Pstar. Tfault=1 normalized analytical time.\n');
fprintf('Nominal uses exact asin(0.3)=%.12f deg with Pstar=1; 17.457603 deg is its rounded label.\n',rad2deg(gammas(1)));
fprintf(['Instantaneous test uses each fixed pre-fault Theta_ref and M_ref=Theta_ref''*Theta_ref.\n' ...
    'F(t)=Theta_ref*A(gamma(t),nu(t),U1(t),Pstar)/Theta_ref, Theta_dot_ref=0.\n' ...
    'This is a consistent fixed-metric differential test away from equilibrium, NOT a re-designed equilibrium metric.\n' ...
    'WARNING: maintained/violated flags are numerical tests at stored samples, including both clearing sides.\n' ...
    'They are not a continuous-time or finite-region contraction proof; recovery is evaluated separately.\n']);
runs=cell(2,nSag);summary=table();
for k=1:2
    ref=point_model(gammas(k),p);
    assert(min(ref.eigM)>tol.matrix && max(ref.eigPaper)<-tol.matrix);
    ps=ref.Pstar;xd=[gammas(k);1];assert(norm(rhs(0,xd,p,ps),inf)==0);
    % Same horizon and sample grids for the nine severities at each point.
    tEnd=max(200,20/abs(max(real(ref.eigA))));
    faultTimes=linspace(0,Tfault,501).';postTimes=linspace(Tfault,tEnd,6001).';
    for z=1:nSag
        pf=p;pf.U1=Usag_values(z);
        opt=odeset('RelTol',1e-9,'AbsTol',1e-11,'MaxStep',.05);
        [tf,xf]=ode15s(@(t,x)rhs(t,x,pf,ps),faultTimes,xd,opt);
        [tp,xp]=ode15s(@(t,x)rhs(t,x,p,ps),postTimes,xf(end,:).',opt);
        % Duplicate t=Tfault records LEFT and RIGHT voltage/Jacobian limits.
        time=[tf;tp];x=[xf;xp];u=[repmat(pf.U1,numel(tf),1);ones(numel(tp),1)];
        fault=[true(numel(tf),1);false(numel(tp),1)];
        assert(max(abs(xf(end,:)-xp(1,:)))<1e-12 && all(x(:,2)>.5));
        gu=unwrap(x(:,1));assert(max(abs(gu-x(:,1)))<1e-10 && all(abs(diff(x(:,1)))<pi));
        gd=rad2deg(gu);errAngle=gd-rad2deg(ref.gamma0);nuerr=x(:,2)-1;
        L=vecnorm(ref.Theta*(x-xd.').',2,1).';
        pe=p.E*u/p.Xeq.*sin(gu);powerResidual=pe-ps;
        % FULL Jacobian at actual state and voltage, with Pstar held fixed.
        a=-p.E*u.*cos(gu)./(p.mass*p.Xeq.*x(:,2));
        d=(-(ps-pe)./x(:,2).^2-p.D)/p.mass;
        ratio=ref.betaDotPaper/ref.beta;
        s11=-ratio*ones(size(time));s22=d+ratio;
        s12=(1+a-d*ratio-ratio^2)/2;
        lambdaMax=(s11+s22)/2+hypot((s11-s22)/2,s12);
        lambdaMinM=repmat(min(ref.eigM),size(time));
        % Independently check selected samples by literal matrix multiplication.
        for ii=unique([1,numel(tf),numel(tf)+1,round(numel(time)/2),numel(time)])
            pInstant=p;pInstant.U1=u(ii);[A,~]=jacobians(gu(ii),pInstant,x(ii,2),ps);
            F=ref.Theta*A/ref.Theta;S=(F+F')/2;
            assert(abs(max(eig(S))-lambdaMax(ii))<1e-10);
        end
        pd=lambdaMinM>tol.matrix;nd=lambdaMax<-tol.matrix;cert=pd & nd;
        % Final 10% by elapsed time, not by sample count (fault grid is denser).
        tail=time>=.9*tEnd;
        recovered=all(abs(errAngle(tail))<tol.angle_deg & abs(nuerr(tail))<tol.frequency_pu);
        % Postfault same-power roots are gamma0+2*pi*n and pi-gamma0+2*pi*n.
        gm=mean(gu(tail));branchCandidates=[ref.gamma0+2*pi*round((gm-ref.gamma0)/(2*pi)), ...
            pi-ref.gamma0+2*pi*round((gm-(pi-ref.gamma0))/(2*pi))];
        [~,bi]=min(abs(branchCandidates-gm));branch=branchCandidates(bi);
        different=abs(rad2deg(branch-ref.gamma0))>=tol.angle_deg;
        another=different && all(abs(rad2deg(gu(tail)-branch))<tol.angle_deg) && ...
            all(abs(nuerr(tail))<tol.frequency_pu) && all(abs(powerResidual(tail))<tol.power_pu);
        if recovered,cls="Recovered original equilibrium";
        elseif another,cls="Transitioned to another equilibrium branch";
        else,cls="Did not recover original equilibrium";end
        postCert=cert(~fault);lastFailed=find(~postCert,1,'last');
        if isempty(lastFailed),recertTime=Tfault;
        elseif lastFailed<numel(tp),recertTime=tp(lastFailed+1);
        else,recertTime=NaN;end
        postSatisfied=all(cert(tail));faultViolated=any(~cert(fault));
        becameSatisfied=faultViolated && postSatisfied && isfinite(recertTime);
        row=table(rad2deg(ref.gamma0),ps,Usag_values(z),100*(1-Usag_values(z)),Tfault, ...
            max(gd),min(gd),max(abs(errAngle)),max(abs(nuerr)),max(L),gd(end),x(end,2), ...
            errAngle(end),nuerr(end),max(lambdaMax),min(lambdaMinM),all(pd),all(nd),all(cert), ...
            faultViolated,postSatisfied,becameSatisfied,recertTime,max(abs(errAngle(tail))), ...
            max(abs(nuerr(tail))),powerResidual(end),cls, ...
            "Fixed reference metric; sampled instantaneous test, not a regional proof", ...
            'VariableNames',{'gamma0_deg','Pstar_pu','Usag_pu','sagDepth_percent','faultDuration_normalized', ...
            'gammaMax_deg','gammaMin_deg','maxAngleExcursion_deg','maxAbsNuDeviation','maxL', ...
            'gammaFinal_deg','nuFinal','finalGammaError_deg','finalNuError','maxLambdaMaxSymF', ...
            'minLambdaMinM','metricRemainedPD','symFRemainedND','contractionMaintained', ...
            'contractionViolatedDuringFault','contractionSatisfiedPostfaultTail', ...
            'contractionBecameSatisfiedAfterClearing','postfaultSatisfiedFromTime', ...
            'tailMaxAngleError_deg','tailMaxAbsNuDeviation','finalPowerResidual_pu','recoveryClass','contractionEvaluationScope'});
        summary=[summary;row]; %#ok<AGROW>
        runs{k,z}=struct('time',time,'state',x,'gammaUnwrapped',gu,'voltage',u,'faultMask',fault, ...
            'xd',xd,'Pstar',ps,'reference',ref,'distance',L,'lambdaMaxSymF',lambdaMax, ...
            'lambdaMinM',lambdaMinM,'s11',s11,'s22',s22,'s12',s12, ...
            'powerResidual',powerResidual,'tailMask',tail,'summary',row);
    end
end
assert(height(summary)==2*nSag);
% Verify the severest case at each point against tighter solver tolerances.
for k=1:2
    r=runs{k,nSag};pf=p;pf.U1=Usag_values(end);tf=r.time(r.faultMask);tp=r.time(~r.faultMask);
    opt=odeset('RelTol',1e-10,'AbsTol',1e-12,'MaxStep',.025);
    [~,xf]=ode15s(@(t,x)rhs(t,x,pf,r.Pstar),tf,r.xd,opt);
    [~,xp]=ode15s(@(t,x)rhs(t,x,p,r.Pstar),tp,xf(end,:).',opt);
    xc=[xf;xp];e=max(abs(r.state-xc),[],'all');assert(e<1e-5);
    runs{k,nSag}.tighterSolverMaxAbsError=e;
    fprintf('Severest sag solver verification at %.9f deg: max state error %.3e.\n',rad2deg(r.reference.gamma0),e);
end
fprintf('\nVOLTAGE SAG SWEEP SUMMARY (all %d runs)\n',2*nSag);disp(summary);
writetable(summary,fullfile(outDir,'voltage_sag_sweep_certified_points.csv'));
sweep=struct('parameters',p,'Usag_values',Usag_values,'Tfault_normalized',Tfault,'tolerances',tol, ...
    'runs',{runs},'summary',summary,'method','Frozen pre-fault metric; full instantaneous Jacobian; sampled signs only');
save(fullfile(outDir,'voltage_sag_sweep_certified_points.mat'),'sweep','-v7.3');
colors=lines(nSag);
% Keep the original seven colors; give the two added severities distinct colors.
colors(8,:)=[.25,.25,.25];colors(9,:)=[.8,.2,.7];
for k=1:2
    fig=figure('Color','w','Position',[80,80,1400,800]);tl=tiledlayout(fig,2,1,'TileSpacing','compact');
    for panel=1:2
        ax=nexttile(tl);hold(ax,'on');
        for z=1:nSag
            r=runs{k,z};if panel==1,y=rad2deg(r.gammaUnwrapped);else,y=r.state(:,2)-1;end
            plot(ax,r.time,y,'Color',colors(z,:),'LineWidth',1.4,'DisplayName',sprintf('Usag = %.1f',Usag_values(z)));
        end
        xline(ax,Tfault,'k--','Clearing','HandleVisibility','off');grid(ax,'on');
        xlim(ax,[0,r.time(end)]);xlabel(ax,'Normalized analytical time');
        if panel==1,ylabel(ax,'gamma unwrapped (deg)');yline(ax,rad2deg(gammas(k)),'k:','HandleVisibility','off');
        else,ylabel(ax,'nu - 1 (pu)');end
        legend(ax,'Location','best','NumColumns',4);
    end
    title(tl,sprintf('Voltage-sag sweep: gamma0 = %.9f deg, fixed Pstar = %.9f pu, duration = 1 normalized time',rad2deg(gammas(k)),runs{k,1}.Pstar));
    save_figure(fig,outDir,sprintf('11_voltage_sag_point_%d',k));
end
fig=figure('Color','w','Position',[90,90,1300,550]);tl=tiledlayout(fig,1,2);
for panel=1:2
    ax=nexttile(tl);hold(ax,'on');
    for k=1:2
        rows=(k-1)*nSag+(1:nSag);if panel==1,y=summary.maxAngleExcursion_deg(rows);else,y=summary.maxAbsNuDeviation(rows);end
        plot(ax,summary.sagDepth_percent(rows),y,'-o','LineWidth',1.5,'DisplayName',sprintf('gamma0 = %.6f deg',rad2deg(gammas(k))));
    end
    xlabel(ax,'Sag depth (%)');grid(ax,'on');legend(ax,'Location','best');
    if panel==1,ylabel(ax,'Maximum angle excursion (deg)');else,ylabel(ax,'Maximum |nu - 1| (pu)');end
end
title(tl,'Discrete voltage-sag sweep; identical fault duration = 1 normalized analytical time');
save_figure(fig,outDir,'12_voltage_sag_severity_summary');
for k=1:2
    rows=(k-1)*nSag+(1:nSag);tab=summary(rows,:);
    fprintf('\nVOLTAGE SAG RESULTS: gamma0 = %.9f deg\n',rad2deg(gammas(k)));
    disp(tab(:,{'Usag_pu','sagDepth_percent','maxAngleExcursion_deg','maxAbsNuDeviation','contractionMaintained','recoveryClass'}));
    recovered=tab.recoveryClass=="Recovered original equilibrium";
    if any(recovered)
        lowest=min(tab.Usag_pu(recovered));
        fprintf('Strongest recovering sag among the tested discrete cases: Usag=%.1f pu, depth=%.0f%%.\n',lowest,100*(1-lowest));
    else,fprintf('Strongest recovering sag among the tested discrete cases: NONE.\n');end
end
fprintf('Recovery requires final 10%% in time within 0.1 deg and 1e-4 pu of ORIGINAL equilibrium.\n');
fprintf('No exact critical sag, voltage boundary or CCT is inferred; no fault-duration sweep was performed.\n');
end

function sweep=voltage_sag_duration_sweep(p,outDir)
% Separate study. No angle jump: each run starts at its PRE-FAULT equilibrium.
% Fault at t=0; clearing at t=Tfault in NORMALIZED analytical time (not seconds).
% Nominal angle uses asin(0.3) to avoid the rounded-angle/Pstar=1 mismatch.
% Original 12-point inputs and every baseline matrix/figure remain unchanged.
gammas=[asin(p.Pstar_nominal/(p.E*p.U1/p.Xeq)),deg2rad(65)];
Tfault_values=[1,2,3,4,5,7,10];Usag=.1;
nDur=numel(Tfault_values);
tol=struct('angle_deg',.1,'frequency_pu',1e-4,'power_pu',1e-4,'matrix',1e-9);
fprintf('\n========================================================\nSEPARATE FAULT-DURATION SWEEP: Usag=0.1 pu\n========================================================\n');
fprintf('No angle jump; initial x=[gamma0;1]; fixed pre-fault Pstar. Tfault=[1 2 3 4 5 7 10] normalized analytical time.\n');
fprintf('Nominal uses exact asin(0.3)=%.12f deg with Pstar=1; 17.457603 deg is its rounded label.\n',rad2deg(gammas(1)));
fprintf(['Instantaneous test uses each fixed pre-fault Theta_ref and M_ref=Theta_ref''*Theta_ref.\n' ...
    'F(t)=Theta_ref*A(gamma(t),nu(t),U1(t),Pstar)/Theta_ref, Theta_dot_ref=0.\n' ...
    'This is a consistent fixed-metric differential test away from equilibrium, NOT a re-designed equilibrium metric.\n' ...
    'WARNING: maintained/violated flags are numerical tests at stored samples, including both clearing sides.\n' ...
    'They are not a continuous-time or finite-region contraction proof; recovery is evaluated separately.\n']);
runs=cell(2,nDur);summary=table();
for k=1:2
    ref=point_model(gammas(k),p);
    assert(min(ref.eigM)>tol.matrix && max(ref.eigPaper)<-tol.matrix);
    ps=ref.Pstar;xd=[gammas(k);1];assert(norm(rhs(0,xd,p,ps),inf)==0);
    % A full postfault settling horizon is provided AFTER every clearing.
    postDuration=max(200,20/abs(max(real(ref.eigA))));
    for z=1:nDur
        Tfault=Tfault_values(z);tEnd=Tfault+postDuration;
        faultTimes=linspace(0,Tfault,ceil(Tfault/.002)+1).';
        postTimes=linspace(Tfault,tEnd,6001).';
        pf=p;pf.U1=Usag;
        opt=odeset('RelTol',1e-9,'AbsTol',1e-11,'MaxStep',.05);
        [tf,xf]=ode15s(@(t,x)rhs(t,x,pf,ps),faultTimes,xd,opt);
        [tp,xp]=ode15s(@(t,x)rhs(t,x,p,ps),postTimes,xf(end,:).',opt);
        % Duplicate t=Tfault records LEFT and RIGHT voltage/Jacobian limits.
        time=[tf;tp];x=[xf;xp];u=[repmat(pf.U1,numel(tf),1);ones(numel(tp),1)];
        fault=[true(numel(tf),1);false(numel(tp),1)];
        assert(max(abs(xf(end,:)-xp(1,:)))<1e-12 && all(x(:,2)>.5));
        gu=unwrap(x(:,1));assert(max(abs(gu-x(:,1)))<1e-10 && all(abs(diff(x(:,1)))<pi));
        gd=rad2deg(gu);errAngle=gd-rad2deg(ref.gamma0);nuerr=x(:,2)-1;
        L=vecnorm(ref.Theta*(x-xd.').',2,1).';
        pe=p.E*u/p.Xeq.*sin(gu);powerResidual=pe-ps;
        % FULL Jacobian at actual state and voltage, with Pstar held fixed.
        a=-p.E*u.*cos(gu)./(p.mass*p.Xeq.*x(:,2));
        d=(-(ps-pe)./x(:,2).^2-p.D)/p.mass;
        ratio=ref.betaDotPaper/ref.beta;
        s11=-ratio*ones(size(time));s22=d+ratio;
        s12=(1+a-d*ratio-ratio^2)/2;
        lambdaMax=(s11+s22)/2+hypot((s11-s22)/2,s12);
        lambdaMinM=repmat(min(ref.eigM),size(time));
        % Independently check selected samples by literal matrix multiplication.
        for ii=unique([1,numel(tf),numel(tf)+1,round(numel(time)/2),numel(time)])
            pInstant=p;pInstant.U1=u(ii);[A,~]=jacobians(gu(ii),pInstant,x(ii,2),ps);
            F=ref.Theta*A/ref.Theta;S=(F+F')/2;
            assert(abs(max(eig(S))-lambdaMax(ii))<1e-10);
        end
        pd=lambdaMinM>tol.matrix;nd=lambdaMax<-tol.matrix;cert=pd & nd;
        % Final 10% by elapsed time, not by sample count (fault grid is denser).
        tail=time>=.9*tEnd;
        recovered=all(abs(errAngle(tail))<tol.angle_deg & abs(nuerr(tail))<tol.frequency_pu);
        % Postfault same-power roots are gamma0+2*pi*n and pi-gamma0+2*pi*n.
        gm=mean(gu(tail));branchCandidates=[ref.gamma0+2*pi*round((gm-ref.gamma0)/(2*pi)), ...
            pi-ref.gamma0+2*pi*round((gm-(pi-ref.gamma0))/(2*pi))];
        [~,bi]=min(abs(branchCandidates-gm));branch=branchCandidates(bi);
        different=abs(rad2deg(branch-ref.gamma0))>=tol.angle_deg;
        another=different && all(abs(rad2deg(gu(tail)-branch))<tol.angle_deg) && ...
            all(abs(nuerr(tail))<tol.frequency_pu) && all(abs(powerResidual(tail))<tol.power_pu);
        if recovered,cls="Recovered original equilibrium";
        elseif another,cls="Transitioned to another equilibrium branch";
        else,cls="Did not recover original equilibrium";end
        postCert=cert(~fault);lastFailed=find(~postCert,1,'last');
        if isempty(lastFailed),recertTime=Tfault;
        elseif lastFailed<numel(tp),recertTime=tp(lastFailed+1);
        else,recertTime=NaN;end
        postSatisfied=all(cert(tail));faultViolated=any(~cert(fault));
        becameSatisfied=faultViolated && postSatisfied && isfinite(recertTime);
        row=table(rad2deg(ref.gamma0),ps,Usag,100*(1-Usag),Tfault, ...
            max(gd),min(gd),max(abs(errAngle)),max(abs(nuerr)),max(L),gd(end),x(end,2), ...
            errAngle(end),nuerr(end),max(lambdaMax),min(lambdaMinM),all(pd),all(nd),all(cert), ...
            faultViolated,postSatisfied,becameSatisfied,recertTime,max(abs(errAngle(tail))), ...
            max(abs(nuerr(tail))),powerResidual(end),cls, ...
            "Fixed reference metric; sampled instantaneous test, not a regional proof", ...
            'VariableNames',{'gamma0_deg','Pstar_pu','Usag_pu','sagDepth_percent','faultDuration_normalized', ...
            'gammaMax_deg','gammaMin_deg','maxAngleExcursion_deg','maxAbsNuDeviation','maxL', ...
            'gammaFinal_deg','nuFinal','finalGammaError_deg','finalNuError','maxLambdaMaxSymF', ...
            'minLambdaMinM','metricRemainedPD','symFRemainedND','contractionMaintained', ...
            'contractionViolatedDuringFault','contractionSatisfiedPostfaultTail', ...
            'contractionBecameSatisfiedAfterClearing','postfaultSatisfiedFromTime', ...
            'tailMaxAngleError_deg','tailMaxAbsNuDeviation','finalPowerResidual_pu','recoveryClass','contractionEvaluationScope'});
        ic=numel(tf);ir=ic+1;
        row.gammaAtClearing_deg=gd(ic);row.nuAtClearing=x(ic,2);
        row.deltaGammaClearing_deg=errAngle(ic);row.deltaNuClearing=nuerr(ic);
        firstBad=find(~cert,1,'first');
        if isempty(firstBad),row.firstContractionViolationTime=NaN;
        else,row.firstContractionViolationTime=time(firstBad);end
        % Integrate positive-time fraction of piecewise-linear sampled lambda.
        % Duplicate clearing samples have dt=0: do not smear the voltage jump.
        row.totalTimeWithPositiveLambdaMaxSymF=positive_lambda_duration(time,lambdaMax);
        row.lambdaMaxSymFAtClearing=lambdaMax(ir); % RIGHT limit, restored U1=1.
        row.lambdaMaxSymFBeforeClearing=lambdaMax(ic); % LEFT limit, U1=0.1.
        row.totalSimulationTime_normalized=tEnd;
        row.postfaultDuration_normalized=postDuration;
        summary=[summary;row]; %#ok<AGROW>
        runs{k,z}=struct('time',time,'state',x,'gammaUnwrapped',gu,'voltage',u,'faultMask',fault, ...
            'xd',xd,'Pstar',ps,'reference',ref,'distance',L,'lambdaMaxSymF',lambdaMax, ...
            'lambdaMinM',lambdaMinM,'s11',s11,'s22',s22,'s12',s12, ...
            'powerResidual',powerResidual,'tailMask',tail,'summary',row);
    end
end
assert(height(summary)==2*nDur);
% Verify each duration with tighter tolerances; classifications remain separate.
for k=1:2
    for z=1:nDur
    r=runs{k,z};pf=p;pf.U1=Usag;tf=r.time(r.faultMask);tp=r.time(~r.faultMask);
    opt=odeset('RelTol',1e-10,'AbsTol',1e-12,'MaxStep',.025);
    [~,xf]=ode15s(@(t,x)rhs(t,x,pf,r.Pstar),tf,r.xd,opt);
    [~,xp]=ode15s(@(t,x)rhs(t,x,p,r.Pstar),tp,xf(end,:).',opt);
    xc=[xf;xp];e=max(abs(r.state-xc),[],'all');assert(e<1e-5);
    runs{k,z}.tighterSolverMaxAbsError=e;
    fprintf('Duration solver verification: gamma0=%.9f deg, Tfault=%.0f, max state error=%.3e.\n',rad2deg(r.reference.gamma0),Tfault_values(z),e);
    end
end
fprintf('\nVOLTAGE SAG DURATION SWEEP SUMMARY (all %d runs)\n',2*nDur);disp(summary);
writetable(summary,fullfile(outDir,'voltage_sag_duration_sweep_Usag_0p1.csv'));
sweep=struct('parameters',p,'Usag_pu',Usag,'Tfault_values_normalized',Tfault_values,'tolerances',tol, ...
    'runs',{runs},'summary',summary,'method','Frozen pre-fault metric; full instantaneous Jacobian; sampled signs only');
save(fullfile(outDir,'voltage_sag_duration_sweep_Usag_0p1.mat'),'sweep','-v7.3');

colors=lines(nDur);
for k=1:2
    fig=figure('Color','w','Position',[80,80,1400,800]);tl=tiledlayout(fig,2,1,'TileSpacing','compact');
    for panel=1:2
        ax=nexttile(tl);hold(ax,'on');
        for z=1:nDur
            r=runs{k,z};ic=find(r.faultMask,1,'last');
            if panel==1,y=rad2deg(r.gammaUnwrapped);else,y=r.state(:,2)-1;end
            plot(ax,r.time,y,'Color',colors(z,:),'LineWidth',1.4,'DisplayName',sprintf('Tfault = %g',Tfault_values(z)));
            plot(ax,r.time(ic),y(ic),'o','Color',colors(z,:),'MarkerSize',5,'HandleVisibility','off');
        end
        xlim(ax,[0,runs{k,end}.time(end)]);grid(ax,'on');xlabel(ax,'Normalized analytical time');
        if panel==1,ylabel(ax,'gamma unwrapped (deg)');yline(ax,rad2deg(gammas(k)),'k:','HandleVisibility','off');
        else,ylabel(ax,'nu - 1 (pu)');end
        legend(ax,'Location','best','NumColumns',4);
    end
    title(tl,sprintf('Duration sweep: gamma0 = %.9f deg, Usag = 0.1 pu; circles mark clearing',rad2deg(gammas(k))));
    save_figure(fig,outDir,sprintf('13_voltage_sag_duration_point_%d',k));
end
fig=figure('Color','w','Position',[80,80,1600,550]);tl=tiledlayout(fig,1,3);
fields={'maxAngleExcursion_deg','maxAbsNuDeviation','deltaGammaClearing_deg'};
ylabels={'Maximum angle excursion (deg)','Maximum |nu - 1| (pu)','Angle displacement at clearing (deg)'};
for panel=1:3
    ax=nexttile(tl);hold(ax,'on');
    for k=1:2
        rows=(k-1)*nDur+(1:nDur);
        plot(ax,Tfault_values,summary.(fields{panel})(rows),'-o','LineWidth',1.5,'DisplayName',sprintf('gamma0 = %.6f deg',rad2deg(gammas(k))));
    end
    grid(ax,'on');xlabel(ax,'Fault duration (normalized analytical time)');ylabel(ax,ylabels{panel});legend(ax,'Location','best');
end
title(tl,'Coarse fault-duration sweep, Usag = 0.1 pu');
save_figure(fig,outDir,'14_voltage_sag_duration_summary');
fig=figure('Color','w','Position',[80,80,1150,600]);ax=axes(fig);hold(ax,'on');
labels=["Recovered original equilibrium","Transitioned to another equilibrium branch","Did not recover original equilibrium"];
for k=1:2
    rows=(k-1)*nDur+(1:nDur);codes=nan(1,nDur);
    for z=1:nDur,codes(z)=find(labels==summary.recoveryClass(rows(z)))-1;end
    markers={'o','s'};
    plot(ax,Tfault_values,codes,['-',markers{k}],'LineWidth',1.5,'MarkerSize',9-k*2, ...
        'DisplayName',sprintf('gamma0 = %.6f deg',rad2deg(gammas(k))));
end
yticks(ax,0:2);yticklabels(ax,labels);ylim(ax,[-.2,2.2]);xticks(ax,Tfault_values);grid(ax,'on');
xlabel(ax,'Fault duration (normalized analytical time)');legend(ax,'Location','best');
title(ax,'Recovery relative to ORIGINAL equilibrium; lines connect discrete tests only');
save_figure(fig,outDir,'15_voltage_sag_duration_recovery');
for k=1:2
    tab=summary((k-1)*nDur+(1:nDur),:);ok=tab.recoveryClass==labels(1);
    fprintf('\nFAULT-DURATION RESULTS: gamma0=%.9f deg, Usag=0.1 pu\n',rad2deg(gammas(k)));
    disp(tab(:,{'faultDuration_normalized','deltaGammaClearing_deg','maxAngleExcursion_deg', ...
        'maxAbsNuDeviation','contractionMaintained','firstContractionViolationTime','recoveryClass'}));
    if any(ok),fprintf('Longest tested duration recovering ORIGINAL equilibrium: %g normalized time units.\n',max(Tfault_values(ok)));
    else,fprintf('No tested duration recovered ORIGINAL equilibrium.\n');end
    bad=find(~ok,1,'first');
    if isempty(bad),fprintf('No recovery loss observed within tested duration range.\n');
    else
        fprintf('First tested duration not recovering original equilibrium: %g normalized time units.\n',Tfault_values(bad));
        for z=2:nDur
            if ok(z-1) && ~ok(z)
                fprintf('Recovery transition is bracketed between Tfault = %g and Tfault = %g normalized time units.\n',Tfault_values(z-1),Tfault_values(z));
            end
        end
    end
end
fprintf(['Contraction flags remain sampled fixed-reference-metric diagnostics, independent of recovery.\n' ...
    'lambdaMaxSymFAtClearing is the post-clearing RIGHT limit; the LEFT limit is also saved.\n' ...
    'Positive-lambda duration is estimated from linear interpolation between samples.\n' ...
    'Only the seven requested durations were tested. No automatic refinement or exact CCT claim.\n']);
end

function duration=positive_lambda_duration(t,lambda)
% Exact positive-duration measure of the piecewise-linear sampled trace.
a=lambda(1:end-1);b=lambda(2:end);fraction=zeros(size(a));
both=a>0 & b>0;fraction(both)=1;
falling=a>0 & b<=0;fraction(falling)=a(falling)./(a(falling)-b(falling));
rising=a<=0 & b>0;fraction(rising)=b(rising)./(b(rising)-a(rising));
duration=sum(diff(t).*fraction);
end


function figure5_gfm_postprocess(base,sweep,outDir)
% New post-processing only: existing trajectories, classifications and intervals
% are read unchanged. Red HORIZONTAL boundaries require an explicit nu rule.
% For each trial gamma, use the largest lambda over ALL stored FAULT-ON nu
% samples at the chosen Usag, with the original fixed Theta_ref and Pstar.
% This intersection of sampled frequency slices is not a regional proof,
% nor is it the pointwise gamma boundary along the time-varying trajectory.
% It is only a sampled instantaneous contraction boundary (fixed metric).
p=sweep.parameters;spec=[1,.5;1,.3;2,.5;2,.3];caseIDs=[2,2,4,4];
I=base.paperIntervals_deg;examples=cell(1,4);companion=table();
gridDeg=linspace(-180,180,36001);g=deg2rad(gridDeg);
fprintf('\nFIGURE-5-STYLE POST-PROCESSING\n');
fprintf(['Blue: previously computed steady-state contraction boundaries (family of reference metrics).\n' ...
    'Red: fixed pre-fault metric; Usag fixed; intersection over all stored fault-on nu samples.\n' ...
    'Gamma scan [-180,180] deg, step 0.01 deg; crossings use linear interpolation only.\n' ...
    'These sampled instantaneous contraction boundaries are not stability limits or a regional proof.\n']);
for q=1:4
    k=spec(q,1);z=find(abs(sweep.Usag_values-spec(q,2))<1e-12,1);r=sweep.runs{k,z};
    ref=r.reference;nuSamples=unique(r.state(r.faultMask,2));U=spec(q,2);
    ratio=ref.betaDotPaper/ref.beta;lambdaEnvelope=-inf(size(g));
    for nu=nuSamples.'
        pe=p.E*U/p.Xeq*sin(g);
        a=-p.E*U*cos(g)/(p.mass*p.Xeq*nu);
        d=(-(r.Pstar-pe)/nu^2-p.D)/p.mass;
        off=(1+a-d*ratio-ratio^2)/2;
        lam=d/2+hypot(-ratio-d/2,off);
        lambdaEnvelope=max(lambdaEnvelope,lam);
    end
    % Connected sign-negative components; no equilibrium-region roots are changed.
    neg=lambdaEnvelope<0 & min(ref.eigM)>sweep.tolerances.matrix;
    starts=find(diff([false,neg,false])==1);stops=find(diff([false,neg,false])==-1)-1;
    pairs=zeros(0,2);
    for h=1:numel(starts)
        lo=starts(h);hi=stops(h);
        if lo==1 || hi==numel(g),continue;end % scan-edge clipping is not a boundary
        x1=gridDeg(lo-1:lo);v1=lambdaEnvelope(lo-1:lo);
        x2=gridDeg(hi:hi+1);v2=lambdaEnvelope(hi:hi+1);
        pairs(end+1,:)=[x1(1)-v1(1)*diff(x1)/diff(v1),x2(1)-v2(1)*diff(x2)/diff(v2)]; %#ok<AGROW>
    end
    bd=[NaN,NaN];status="No clean boundary pair";
    if ~isempty(pairs)
        % Choose the component nearest the actual fault-on angular excursion.
        faultG=rad2deg(r.gammaUnwrapped(r.faultMask));target=mean([min(faultG),max(faultG)]);
        distance=max(pairs(:,1)-target,0)+max(target-pairs(:,2),0);
        [~,h]=min(distance);bd=pairs(h,:);status="Nearest complete sampled component";
    else
        fprintf('Case %d, Usag=%.1f: No clean sampled instantaneous gamma boundary found for this case with the current fixed-reference-metric formulation.\n',caseIDs(q),U);
    end
    region=find(rad2deg(ref.gamma0)>I(:,1) & rad2deg(ref.gamma0)<I(:,2),1);
    assert(~isempty(region));
    row=r.summary(:,{'gamma0_deg','Usag_pu','faultDuration_normalized','gammaMax_deg', ...
        'maxAngleExcursion_deg','maxAbsNuDeviation','contractionMaintained','recoveryClass'});
    row.CaseLabel="Case "+string(caseIDs(q));row.Example=string(char('A'+q-1));
    row.sampledGammaMin_deg=bd(1);row.sampledGammaMax_deg=bd(2);
    row.faultNuMin=min(nuSamples);row.faultNuMax=max(nuSamples);row.boundaryStatus=status;
    row.caseLabel=row.CaseLabel;row.regionLabel="Region "+string(region);
    row.gamma0_region_min_deg=I(region,1);row.gamma0_region_max_deg=I(region,2);
    row.gammaInstantaneousMin_deg=bd(1);row.gammaInstantaneousMax_deg=bd(2);
    row.gammaTrajectoryMin_deg=min(rad2deg(r.gammaUnwrapped));
    row.gammaTrajectoryMax_deg=max(rad2deg(r.gammaUnwrapped));
    row.boundaryEvaluationScope="Fixed reference metric; intersection over fault-on nu samples; not a regional proof";
    companion=[companion;row]; %#ok<AGROW>
    examples{q}=struct('run',r,'caseID',caseIDs(q),'region',region,'steadyBoundaries_deg',I(region,:), ...
        'sampledBoundaries_deg',bd,'allSampledComponents_deg',pairs,'gammaScan_deg',gridDeg, ...
        'worstSampledLambda',lambdaEnvelope,'faultNuSamples',nuSamples,'boundaryStatus',status);
end
boundarySummary=companion(:,{'caseLabel','regionLabel','gamma0_deg','Usag_pu', ...
    'gamma0_region_min_deg','gamma0_region_max_deg','gammaInstantaneousMin_deg','gammaInstantaneousMax_deg', ...
    'gammaTrajectoryMin_deg','gammaTrajectoryMax_deg','contractionMaintained','recoveryClass','boundaryEvaluationScope'});
disp(boundarySummary);
fig=figure('Color','w','Position',[60,100,1600,550]);ax=axes(fig);hold(ax,'on');ylim(ax,[-.7,1.5]);xlim(ax,[0,125]);
shade_intervals(ax,I);
for k=1:size(I,1)
    for h=1:2
        boundary=I(k,h);xline(ax,boundary,'b--','HandleVisibility','off');
        boundNames={'min','max'};
        text(ax,boundary,1.05+.25*(h-1),sprintf('gamma0_%s_R%d = %.6f deg',boundNames{h},k,boundary), ...
            'HorizontalAlignment','center','FontSize',8,'BackgroundColor','w','Interpreter','none');
    end
    text(ax,mean(I(k,:)),.75,sprintf('Steady-state contraction-certified\nRegion %d',k),'HorizontalAlignment','center','FontWeight','bold');
end
yline(ax,0,'k-','HandleVisibility','off');
for j=1:5
    r=base.cases{j};ma=max(real(r.eigA));cert=min(r.eigM)>base.classificationTolerance && max(r.eigPaper)<-base.classificationTolerance;
    if cert,color=[0,.55,.2];marker='o';elseif ma>base.classificationTolerance,color=[.85,.1,.1];marker='x';else,color=[.9,.5,.05];marker='s';end
    plot(ax,rad2deg(r.gamma0),0,marker,'Color',color,'MarkerSize',10,'LineWidth',2,'HandleVisibility','off');
    % Stagger nearby Case 1/2 labels for clarity.
    if mod(j,2),height=.3;else,height=-.35;end
    text(ax,rad2deg(r.gamma0),height,sprintf('Case %d\n%.6g deg',j,rad2deg(r.gamma0)),'HorizontalAlignment','center');
end
plot(ax,NaN,NaN,'s','Color',[.9,.5,.05],'DisplayName','Locally stable, not certified');
plot(ax,NaN,NaN,'o','Color',[0,.55,.2],'DisplayName','Contraction-certified');
plot(ax,NaN,NaN,'x','Color',[.85,.1,.1],'DisplayName','Locally unstable');
yticks(ax,[]);xlabel(ax,'gamma0 (deg)');legend(ax,'Location','southoutside','Orientation','horizontal');
title(ax,'Figure 5A - Steady-State Contraction-Certified Operating Regions');
save_figure(fig,outDir,'figure5A_gamma0_contraction_regions');
exportgraphics(fig,fullfile(outDir,'figure5A_gamma0_contraction_regions.pdf'),'ContentType','vector');
% Common y limits within each reference pair include every displayed boundary.
limits=zeros(2,2);
for k=1:2
    vals=[];
    for q=2*k-1:2*k
        e=examples{q};vals=[vals,e.sampledBoundaries_deg,rad2deg(e.run.gammaUnwrapped).']; %#ok<AGROW>
    end
    vals=vals(isfinite(vals));pad=.08*(max(vals)-min(vals));limits(k,:)=[min(vals)-pad,max(vals)+pad];
end
fig=figure('Color','w','Position',[50,40,1650,1450]);tl=tiledlayout(fig,3,2,'TileSpacing','compact');
% TOP: extend the existing +1-degree TEST display to show Case-5 escape.
% Only display data are selected; full simulations and diagnostics are intact.
displayEnd=75;
axtop=nexttile(tl,[1,2]);hold(axtop,'on');cmap=lines(5);
positions={'below Region 1','inside Region 1','between regions','inside Region 2','above Region 2'};
displayedGamma=[];
for j=1:5
    tr=base.trajectories{j};ref=base.cases{j};
    displayMask=tr.time<=displayEnd;
    gammaDisplay=rad2deg(tr.gammaUnwrapped(displayMask));
    displayedGamma=[displayedGamma;gammaDisplay(:)]; %#ok<AGROW>
    plot(axtop,tr.time(displayMask),gammaDisplay,'Color',cmap(j,:),'LineWidth',1.6, ...
        'DisplayName',sprintf('Case %d: %.6g deg (%s)',j,rad2deg(ref.gamma0),positions{j}));
    plot(axtop,tr.time(1),rad2deg(tr.gammaUnwrapped(1)),'o','Color',cmap(j,:), ...
        'MarkerSize',4,'HandleVisibility','off');
end
for k=1:2
    for h=1:2,yline(axtop,I(k,h),'b--','LineWidth',1,'HandleVisibility','off');end
end
xlabel(axtop,'Normalized analytical time');ylabel(axtop,'gamma (deg)');grid(axtop,'on');
title(axtop,'Five representative operating points: gamma(t) and Case 5 departure');
xlim(axtop,[0,displayEnd]);
% Reserve a clear band above the curves for the five-case legend.
ylim(axtop,[0,1.35*max([displayedGamma;I(:)])]);
legend(axtop,'Location','northwest','NumColumns',3,'FontSize',9);
mainAxes=gobjects(1,4);
for q=1:4
    e=examples{q};r=e.run;ax=nexttile(tl);mainAxes(q)=ax;hold(ax,'on');
    plot(ax,r.time,rad2deg(r.gammaUnwrapped),'k','LineWidth',1.6,'DisplayName','Voltage-sag gamma(t)');
    trBase=base.trajectories{e.caseID};
    plot(ax,trBase.time,rad2deg(trBase.gammaUnwrapped),'--','Color',[.45,.45,.45], ...
        'LineWidth',.9,'DisplayName','Baseline +1 deg test');
    yline(ax,rad2deg(r.reference.gamma0),'Color',[.4,.4,.4],'LineWidth',.7,'DisplayName','Original gamma0');
    for h=1:2
        if h==1,visibility='on';else,visibility='off';end
        if isfinite(e.sampledBoundaries_deg(h))
            yline(ax,e.sampledBoundaries_deg(h),'r--','LineWidth',1.3, ...
                'DisplayName','Sampled instantaneous bounds','HandleVisibility',visibility);
        end
    end
    ylim(ax,limits(spec(q,1),:));
    % Focus on disturbance and recovery; full data remain unchanged.
    timeWindows=[40,80];xlim(ax,[0,timeWindows(spec(q,1))]);yl=ylim(ax);tf=r.summary.faultDuration_normalized;
    patch(ax,[0,tf,tf,0],[yl(1),yl(1),yl(2),yl(2)],[1,.75,.3], ...
        'FaceAlpha',.3,'EdgeColor','none','DisplayName','Fault on');
    xline(ax,tf,':','Color',[.6,.4,0],'HandleVisibility','off');
    if r.summary.contractionMaintained,condition='Sampled contraction maintained';
    elseif r.summary.contractionViolatedDuringFault,condition='Sampled contraction violated during fault';
    else,condition='Sampled contraction not maintained';end
    title(ax,sprintf('%c: Case %d, gamma0=%.6f deg, Usag=%.1f, Tfault=1',char('A'+q-1),e.caseID,rad2deg(r.reference.gamma0),r.summary.Usag_pu));
    subtitle(ax,{condition,char(r.summary.recoveryClass)});
    if any(~isfinite(e.sampledBoundaries_deg)),text(ax,.02,.55,'No clean sampled boundary pair','Units','normalized','Color','r');end
    xlabel(ax,'Normalized analytical time');ylabel(ax,'gamma (deg)');grid(ax,'on');
    legend(ax,'Location','northwest','FontSize',8);
end
title(tl,{'Figure 5B - Temporal gamma(t) for five operating points, with voltage-sag examples for certified cases', ...
    'Red dashed: sampled instantaneous contraction boundaries, fixed reference metric; not a regional proof'});
% Early-time insets expose the sag response despite wide gamma-boundary ranges.
drawnow;
for q=1:4
    pos=mainAxes(q).Position;e=examples{q};r=e.run;k=spec(q,1);
    insetY=.35;if q>2,insetY=.12;end
    axz=axes(fig,'Position',[pos(1)+.57*pos(3),pos(2)+insetY*pos(4),.39*pos(3),.30*pos(4)]);
    hold(axz,'on');sel=r.time<=5;
    plot(axz,r.time(sel),rad2deg(r.gammaUnwrapped(sel)),'k','LineWidth',1.4);
    trBase=base.trajectories{e.caseID};ss=trBase.time<=5;
    plot(axz,trBase.time(ss),rad2deg(trBase.gammaUnwrapped(ss)),'--','Color',[.45,.45,.45],'LineWidth',.9);
    pair=companion(2*k-1:2*k,:);lo=min(pair.gamma0_deg)-.15;hi=max(pair.gammaMax_deg)+.15;
    ylim(axz,[lo,hi]);xlim(axz,[0,5]);
    patch(axz,[0,1,1,0],[lo,lo,hi,hi],[1,.75,.3],'FaceAlpha',.3,'EdgeColor','none');
    for y=e.sampledBoundaries_deg
        if y>lo && y<hi,yline(axz,y,'r--');end
    end
    yline(axz,rad2deg(r.reference.gamma0),'Color',[.4,.4,.4]);
    title(axz,'Early time: 0 to 5','FontSize',8);set(axz,'FontSize',8);grid(axz,'on');box(axz,'on');
end
% Top-panel inset: Cases 1-4 retain their colors and steady-state guides.
% Position it in the open area above the early Case-5 trace, below the legend.
posTop=axtop.Position;
axTopZoom=axes(fig,'Position',[posTop(1)+.13*posTop(3), ...
    posTop(2)+.35*posTop(4),.33*posTop(3),.32*posTop(4)]);
hold(axTopZoom,'on');
for j=1:4
    tr=base.trajectories{j};sel=tr.time<=15;
    plot(axTopZoom,tr.time(sel),rad2deg(tr.gammaUnwrapped(sel)), ...
        'Color',cmap(j,:),'LineWidth',1.3);
end
for k=1:2
    for h=1:2,yline(axTopZoom,I(k,h),'b--','LineWidth',.7);end
end
xlim(axTopZoom,[0,15]);ylim(axTopZoom,[8,86]);
title(axTopZoom,'Cases 1-4: early time 0-15','FontSize',8);
set(axTopZoom,'FontSize',8);grid(axTopZoom,'on');box(axTopZoom,'on');
save_figure(fig,outDir,'figure5B_temporal_gamma_with_sag_examples');
exportgraphics(fig,fullfile(outDir,'figure5B_temporal_gamma_with_sag_examples.pdf'),'ContentType','vector');
writetable(boundarySummary,fullfile(outDir,'figure5B_boundary_summary.csv'));
end

function sweep=hd_parameter_sweep(outDir)
% Additive GFM adaptation of Fig. 6/7: no existing data or tuning is changed.
% D_eff is swept directly; K=20 and Kd=D_eff-K is only its interpretation.
% Instantaneous slices use nu=U1=1, fixed midpoint metric and fixed Pstar.
Hc=[.1:.1:1,1.25:.25:6];Dc=20:2:70;
Hf=2.5:.05:4.5;Df=20:.25:35;
fprintf('\nH-D_eff GFM contraction-boundary sensitivity (two-stage sweep)\n');
coarse=hd_grid(Hc,Dc,.25);fine=hd_grid(Hf,Df,.1);
highlight=[hd_pair(.1,60,.1);hd_pair(3.4,25,.1)];
names={'H','D_eff','numberOfCertifiedRegions','gamma0_min_R1','gamma0_max_R1', ...
    'gamma0_min_R2','gamma0_max_R2','gamma_min_R1','gamma_max_R1','gamma_min_R2','gamma_max_R2'};
fineTable=array2table(fine.rows,'VariableNames',names);
highlightTable=array2table(highlight,'VariableNames',names);
% Independent denser detection and direct matrix checks at both highlighted tunings.
for j=1:2
    dense=hd_pair(highlight(j,1),highlight(j,2),.005);
    assert(isequal(isnan(dense),isnan(highlight(j,:))));
    ok=isfinite(dense);assert(max(abs(dense(ok)-highlight(j,ok)))<1e-6);
end
b=load(fullfile(outDir,'baseline_gfm_results.mat'),'results');
assert(max(abs(highlight(2,4:7)-reshape(b.results.paperIntervals_deg.',1,[])))<1e-6, ...
    'Selected tuning must reproduce the verified contraction intervals.');
fprintf('Highlighted tunings: angle boundaries in degrees\n');disp(highlightTable);
fprintf(['Instantaneous domain: [-180,180] deg; connected component containing the regional midpoint.\n' ...
    'Fixed reference metric, U1=nu=1, fixed midpoint Pstar; sampled instantaneous test, not a regional proof.\n' ...
    'Regions sorted by increasing angle; absent regions and unclosed instantaneous pairs are NaN.\n']);
sweep=struct('coarse',coarse,'fine',fine,'fineTable',fineTable,'highlightTable',highlightTable, ...
    'K',20,'E',1,'Xeq',.3,'U1',1,'nu',1,'instantaneousDomain_deg',[-180,180], ...
    'scope','Sampled instantaneous contraction boundaries (fixed reference metric); not a regional proof', ...
    'regionRule','Connected intervals in increasing gamma0 order; missing regions NaN', ...
    'validation','Highlighted tunings cross-checked with 0.005 deg detection and direct Theta*A/Theta');
writetable(fineTable,fullfile(outDir,'HD_eff_sweep_contraction_boundaries_fine.csv'));
writetable(highlightTable,fullfile(outDir,'HD_eff_sweep_highlighted_tunings.csv'));
save(fullfile(outDir,'HD_eff_sweep_contraction_boundaries.mat'),'sweep','-v7.3');
fprintf('Selected tuning: H=3.4, D_eff=25; certified regions=%d.\n',highlight(2,3));
fprintf('Original tuning: H=0.1, D_eff=60.\n');
fprintf('H-D_eff sensitivity sweep completed. Numerical results and exports retained; H-D_eff visualization disabled.\n');
end

function grid=hd_grid(Hvalues,Dvalues,step)
rows=nan(numel(Hvalues)*numel(Dvalues),11);n=0;
for h=Hvalues
    for d=Dvalues,n=n+1;rows(n,:)=hd_pair(h,d,step);end
end
grid=struct('H',Hvalues,'D_eff',Dvalues,'gamma0DetectionStep_deg',step,'rows',rows);
fprintf('Completed %d parameter pairs (gamma0 detection step %.3g deg).\n',n,step);
end

function row=hd_pair(H,D,step)
C=1/(2*H*.3);b=D/(2*H);
f=@(g)hd_steady_lambda(g,C,b);
% Endpoints: gamma0=0 has singular metric; gamma0=90 cannot be strictly ND.
% Small numerical endpoint offset is used only to bracket roots, never as a certificate.
I=hd_intervals(f,[1e-7,step:step:90],true);
row=[H,D,size(I,1),nan(1,8)];
for region=1:min(2,size(I,1))
    g0=mean(I(region,:));beta=-C*sind(g0);c=-C*cosd(g0);
    Theta=[-beta,0;-c,-beta];M=Theta'*Theta;
    [~,pd]=chol(M);assert(pd==0,'Midpoint metric must be positive definite.');
    row(2+2*region:3+2*region)=I(region,:);
    r=c/beta;ps=sin(deg2rad(g0))/.3;
    instant=@(g)hd_instant_lambda(g,C,b,r,ps,H);
    J=hd_intervals(instant,-180:step:180,false);
    ii=find(J(:,1)<g0 & J(:,2)>g0,1);
    if ~isempty(ii) && J(ii,1)>-180+1e-7 && J(ii,2)<180-1e-7
        row(6+2*region:7+2*region)=J(ii,:);
    end
    % Audit the scalar lambda expression against the actual fixed-metric matrix.
    for g=[g0,g0+1]
        a=-C*cosd(g);d=-b-(ps-sind(g)/.3)/(2*H);
        F=Theta*[0,1;a,d]/Theta;direct=max(eig((F+F')/2));
        assert(abs(direct-instant(g))<1e-8*max(1,abs(direct)));
    end
end
end

function z=hd_steady_lambda(g,C,b)
r=cotd(g);off=(1-C*cosd(g)+b*r-r.^2)/2;
z=-b/2+hypot(b/2-r,off);
z(g<=0 | g>=90)=max(z(g<=0 | g>=90),eps);
end

function z=hd_instant_lambda(g,C,b,r,ps,H)
a=-C*cosd(g);d=-b-(ps-sind(g)/.3)/(2*H);
z=d/2+hypot(-r-d/2,(1+a-d*r-r^2)/2);
end

function I=hd_intervals(f,g,searchNarrow)
y=f(g);extra=[];
% Refine local minima as well as sign changes: a narrow certified island can
% fall wholly between grid nodes (especially for the original low-H tuning).
valleys=find(y(2:end-1)<=y(1:end-2) & y(2:end-1)<=y(3:end))+1;
if searchNarrow,valleys=unique([2,valleys,numel(g)-1]);end
for j=valleys
    [x,v]=fminbnd(f,g(j-1),g(j+1),optimset('TolX',1e-11,'Display','off'));
    if v<0,extra(end+1)=x;end %#ok<AGROW>
end
% Also resolve narrow uncertified gaps so nearby intervals are not merged
% merely because the initial angular grid skips a positive local maximum.
peaks=find(y(2:end-1)>=y(1:end-2) & y(2:end-1)>=y(3:end))+1;
for j=peaks
    [x,v]=fminbnd(@(x)-f(x),g(j-1),g(j+1),optimset('TolX',1e-11,'Display','off'));
    if v<0,extra(end+1)=x;end %#ok<AGROW>
end
g=unique([g,extra]);y=f(g);roots=g(y==0);
idx=find((y(1:end-1)<0)~=(y(2:end)<0));
for j=idx
    roots(end+1)=fzero(f,g(j:j+1),optimset('TolX',1e-11,'Display','off')); %#ok<AGROW>
end
edges=unique([g(1),roots,g(end)]);I=zeros(0,2);
for j=1:numel(edges)-1
    if f(mean(edges(j:j+1)))<0,I(end+1,:)=edges(j:j+1);end %#ok<AGROW>
end
end

function hd_boundary_sanity(sweep,outDir)
% Additive audit only: never overwrite numerical sweeps or final 3D surfaces.
% The saved instantaneous domain ALREADY equals [-180,180]; this is a denser
% independent re-search, not a fictitious expansion of a narrower old domain.
oldDomain=sweep.instantaneousDomain_deg;assert(isequal(oldDomain,[-180,180]));
logID=fopen(fullfile(outDir,'HD_boundary_sanity_checks_log.txt'),'w');
cleanup=onCleanup(@()fclose(logID)); %#ok<NASGU>
hd_sanity_log(logID,'Current instantaneous gamma scan: gamma_scan_min=%g deg; gamma_scan_max=%g deg.\n',oldDomain);
hd_sanity_log(logID,'Original steps: coarse=%.3g, fine=%.3g deg. Diagnostic step=0.05 deg with root refinement.\n', ...
    sweep.coarse.gamma0DetectionStep_deg,sweep.fine.gamma0DetectionStep_deg);
hd_sanity_log(logID,'No wider range is requested: the existing range is already +/-180 deg.\n');
rows=[sweep.coarse.rows;sweep.fine.rows;sweep.highlightTable{:,:}];
stages=[repmat("coarse",size(sweep.coarse.rows,1),1);repmat("fine",size(sweep.fine.rows,1),1);repmat("highlight",2,1)];
steps=[repmat(sweep.coarse.gamma0DetectionStep_deg,size(sweep.coarse.rows,1),1); ...
    repmat(sweep.fine.gamma0DetectionStep_deg,size(sweep.fine.rows,1),1);.1;.1];
N=size(rows,1);gap=rows(:,6)-rows(:,5);fresh=nan(N,4);limited=false(N,4);
status=strings(N,4);rootResidual=nan(N,4);crossing=false(N,4);intervals=cell(N,2);
for n=1:N
    H=rows(n,1);D=rows(n,2);C=1/(2*H*.3);b=D/(2*H);tol=max(.25,2*steps(n));
    for region=1:2
        gc=4+2*(region-1);ic=8+2*(region-1);dest=2*region-1:2*region;
        old=rows(n,ic:ic+1);limited(n,dest)=abs(old-oldDomain)<=tol;
        if any(~isfinite(rows(n,gc:gc+1)))
            status(n,dest)="Unavailable: reference region absent";continue
        end
        g0=mean(rows(n,gc:gc+1));r=cotd(g0);ps=sind(g0)/.3;
        beta=-C*sind(g0);Theta=[-beta,0;C*cosd(g0),-beta];
        [~,pd]=chol(Theta'*Theta);assert(pd==0);
        fun=@(g)hd_instant_lambda(g,C,b,r,ps,H);
        J=hd_intervals(fun,-180:.05:180,false);intervals{n,region}=J;
        % Select ONLY the connected interval containing this exact midpoint.
        % Never union disconnected components or choose another periodic branch.
        k=find(J(:,1)<g0 & J(:,2)>g0,1);
        if ~isempty(k)
            for edge=1:2
                q=dest(edge);z=J(k,edge);
                if z>-180+1e-7 && z<180-1e-7
                    delta=1e-5;
                    if edge==1,cross=fun(z-delta)>0 && fun(z+delta)<0;
                    else,cross=fun(z-delta)<0 && fun(z+delta)>0;end
                    fresh(n,q)=z;rootResidual(n,q)=abs(fun(z));crossing(n,q)=cross;
                end
            end
        end
        for edge=1:2
            q=dest(edge);
            if ~isfinite(fresh(n,q))
                status(n,q)="No clean finite boundary found in extended range";
            elseif crossing(n,q) && isfinite(old(edge)) && abs(fresh(n,q)-old(edge))<1e-6
                status(n,q)="True internal boundary";
            elseif crossing(n,q) && limited(n,q) && abs(fresh(n,q)-old(edge))>tol
                status(n,q)="Previous scan-limit truncation";
            else
                status(n,q)="Boundary discrepancy: further investigation required";
            end
        end
    end
end
names={'H','D_eff','numberOfCertifiedRegions','gamma0_min_R1','gamma0_max_R1','gamma0_min_R2','gamma0_max_R2'};
checks=array2table(rows(:,1:7),'VariableNames',names);checks.stage=stages;checks.regionGap_deg=gap;
checks.scanLimitTolerance_deg=max(.25,2*steps);
for region=1:2
    q=2*region-1:2*region;oldcol=8+2*(region-1);
    for edge=1:2
        if edge==1,word='min';capital='Min';else,word='max';capital='Max';end
        checks.(sprintf('old_gamma_%s_R%d',word,region))=rows(:,oldcol+edge-1);
        checks.(sprintf('new_gamma_%s_R%d',word,region))=fresh(:,q(edge));
        checks.(sprintf('gamma%sR%d_scanLimited',capital,region))=limited(:,q(edge));
        checks.(sprintf('gamma%sR%d_status',capital,region))=status(:,q(edge));
        checks.(sprintf('gamma%sR%d_rootResidual',capital,region))=rootResidual(:,q(edge));
        checks.(sprintf('gamma%sR%d_signCrossing',capital,region))=crossing(:,q(edge));
    end
end
% Explicit comparisons on every horizontal/vertical neighbor, not the plot.
edges=[hd_neighbor_audit(sweep.coarse,"coarse");hd_neighbor_audit(sweep.fine,"fine")];
% A 10-degree neighbor jump is declared the diagnostic 'large cliff'.
% Export every edge too, so this threshold is auditable rather than hidden.
large=edges.upperBoundaryJump_deg>=10;merge=edges.isTwoToOne;
fine=edges.stage=="fine";
hd_sanity_log(logID,'Neighbor comparison: large cliff threshold = 10 deg.\n');
hd_sanity_log(logID,'Fine grid: %d 2-to-1 edges; %d large cliffs; %d large cliffs away from topology change.\n', ...
    nnz(fine&merge),nnz(fine&large),nnz(fine&large&~merge));
hd_sanity_log(logID,'Coarse grid: %d 2-to-1 edges; %d large cliffs; %d large cliffs away from topology change.\n', ...
    nnz(~fine&merge),nnz(~fine&large),nnz(~fine&large&~merge));
% Refine three fine-grid transition brackets selected across the boundary.
trans=find(fine&merge);[~,order]=sort(edges.H_two(trans));trans=trans(order);
chosen=trans(unique(round(linspace(1,numel(trans),min(3,numel(trans))))));
examples=table();evidence=true;
for j=1:numel(chosen)
    e=edges(chosen(j),:);p2=[e.H_two,e.D_two];p1=[e.H_one,e.D_one];lo=0;hi=1;
    for it=1:28
        t=(lo+hi)/2;p=p2+t*(p1-p2);a=hd_pair(p(1),p(2),.025);
        if a(3)==2,lo=t;elseif a(3)==1,hi=t;else,error('Unexpected transition topology.');end
    end
    tc=(lo+hi)/2;
    % 'At' is the numerically bracketed transition; strict ND excludes the
    % touching point itself, so its region count is tolerance-sensitive.
    ts=[max(0,tc-.01),tc,min(1,tc+.01)];phase=["before","at (numeric bracket)","after"];
    for z=1:3
        p=p2+ts(z)*(p1-p2);a=hd_pair(p(1),p(2),.01);
        tab=array2table(a(1:7),'VariableNames',names);tab.regionGap_deg=a(6)-a(5);
        tab.transitionExample=j;tab.position=phase(z);tab.bracketWidthParameterFraction=hi-lo;
        examples=[examples;tab]; %#ok<AGROW>
        hd_sanity_log(logID,'Example %d %s: H=%.10g D_eff=%.10g count=%d R1=[%.9g, %.9g] R2=[%.9g, %.9g] gap=%.9g deg.\n', ...
            j,char(phase(z)),a(1:7),a(6)-a(5));
    end
    p=p2+lo*(p1-p2);before=hd_pair(p(1),p(2),.01);
    p=p2+hi*(p1-p2);after=hd_pair(p(1),p(2),.01);
    evidence=evidence && before(3)==2 && after(3)==1 && before(6)-before(5)<.01 ...
        && abs(before(4)-after(4))<1e-4 && abs(before(7)-after(5))<1e-4;
end
coincides=any(fine&large&merge) && ~any(fine&large&~merge) && evidence;
if coincides,conclusion='Region-1 upper-boundary cliff coincides with merging of the two certified intervals.';
else,conclusion='Region-1 upper-boundary cliff does NOT coincide with region merging and requires further boundary-tracking investigation.';end
hd_sanity_log(logID,'\n1) %s\n',conclusion);
truncated=status=="Previous scan-limit truncation";
hd_sanity_log(logID,'2) Scan-limit-near flags=%d; definitively truncated boundaries=%d.\n',nnz(limited),nnz(truncated));
if ~any(truncated(:)),hd_sanity_log(logID,'3) No scan-truncated surfaces / H-D_eff areas identified.\n');
else,hd_sanity_log(logID,'3) Affected areas are explicitly flagged in HD_boundary_sanity_checks.csv.\n');end
hd_sanity_log(logID,'Internal crossing checks passed: %d / %d finite re-found boundaries.\n',nnz(crossing),nnz(isfinite(fresh)));
hd_sanity_log(logID,'Maximum change old vs re-found finite boundaries: %.4g deg.\n',max(abs(fresh(:)-reshape(rows(:,8:11),[],1)),[],'omitnan'));
sel=N;hd_sanity_log(logID,'4) Selected tuning H=3.4, D_eff=25:\n');
for region=1:2
    gc=4+2*(region-1);q=2*region-1:2*region;
    hd_sanity_log(logID,'Region %d: gamma0=[%.12f, %.12f] deg; instantaneous=[%.12f, %.12f] deg.\n', ...
        region,rows(sel,gc:gc+1),fresh(sel,q));
    hd_sanity_log(logID,'Lower: %s; upper: %s; scan-limited flags=[%d %d].\n', ...
        char(status(sel,q(1))),char(status(sel,q(2))),limited(sel,q));
end
hd_sanity_log(logID,'No finite endpoint is replaced by a search edge. Unavailable endpoint means: No finite sampled instantaneous boundary identified within +/-180 deg.\n');
hd_sanity_log(logID,'Final 3D figures and original sweep files remain unchanged; diagnostics are separate.\n');
writetable(checks,fullfile(outDir,'HD_boundary_sanity_checks.csv'));
writetable(edges,fullfile(outDir,'HD_region_merge_neighbor_diagnostics.csv'));
writetable(examples,fullfile(outDir,'HD_region_merge_examples.csv'));
save(fullfile(outDir,'HD_boundary_sanity_checks.mat'),'checks','edges','examples','intervals','conclusion','oldDomain','evidence');
hd_sanity_log(logID,'Selected tuning numberOfCertifiedRegions=%d; original tuning H=0.1, D_eff=60.\n',rows(sel,3));
hd_sanity_log(logID,'H-D_eff sensitivity sweep completed. Numerical results and exports retained; H-D_eff visualization disabled.\n');
end

function t=hd_neighbor_audit(grid,stage)
r=grid.rows;nd=numel(grid.D_eff);nh=numel(grid.H);pairs=[];
for h=1:nh
    for d=1:nd
        k=(h-1)*nd+d;
        if d<nd,pairs(end+1,:)=[k,k+1];end %#ok<AGROW>
        if h<nh,pairs(end+1,:)=[k,k+nd];end %#ok<AGROW>
    end
end
a=r(pairs(:,1),:);b=r(pairs(:,2),:);
m=(a(:,3)==2 & b(:,3)==1)|(a(:,3)==1 & b(:,3)==2);
t=table(repmat(stage,size(a,1),1),a(:,1),a(:,2),b(:,1),b(:,2),a(:,3),b(:,3), ...
    abs(a(:,5)-b(:,5)),m,'VariableNames',{'stage','H_A','D_A','H_B','D_B','countA','countB','upperBoundaryJump_deg','isTwoToOne'});
t.H_two=nan(height(t),1);t.D_two=t.H_two;t.H_one=t.H_two;t.D_one=t.H_two;
for k=find(m).'
    if a(k,3)==2,two=a(k,:);one=b(k,:);else,two=b(k,:);one=a(k,:);end
    t.H_two(k)=two(1);t.D_two(k)=two(2);t.H_one(k)=one(1);t.D_one(k)=one(2);
end
end

function hd_sanity_log(fid,varargin)
fprintf(1,varargin{:});fprintf(fid,varargin{:});
end
