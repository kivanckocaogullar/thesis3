from pathlib import Path
p=Path('outputs/baseline_gfm_contraction.m');s=p.read_text(encoding='utf-8')
a=s.index('selected=nan(1,3);');b=s.index('% Previously verified OLD',a)
s=s[:a]+'''% These are 12 different steady-state operating points, NOT disturbance magnitudes.
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
fprintf('\\n========================================================\\n12-POINT VSM-GFM OPERATING-POINT CLASSIFICATION\\n========================================================\\n');
fprintf('H=%.1f s; K=%.0f; Kd=%.0f; D_eff=%.0f; classification tol=%.3g\\n',H,K,Kd,D,tol);
print_intervals('Verified contraction intervals',intervals);disp(pointResults(:,1:11));
writetable(pointResults(:,1:11),fullfile(outDir,'12_point_gfm_classification.csv'));
save(fullfile(outDir,'12_point_gfm_classification.mat'),'pointResults','pointModels','tol');
representativeIndices=[2,4,8,10,12];selected=deg2rad(testAngles(representativeIndices));
caseNames={'Below R1 - 12 deg','Inside R1 - nominal 1 pu','Between R1/R2 - 35 deg','Inside R2 - 65 deg','Above R2 - 120 deg'};
caseLabels=cellstr(pointResults.stabilityClass(representativeIndices));
fprintf('Five representative trajectory angles (deg):\\n');disp(testAngles(representativeIndices));
fprintf('COMMON perturbation for all five: gamma(0)=gamma0+1 deg; nu(0)=1. No severity sweep.\\n');
''' + s[b:]
s=s.replace('cases=cell(1,3);trajectories=cell(1,3);','cases=cell(1,5);trajectories=cell(1,5);')
s=s.replace('s=point_model(selected(j),par);s.name=', 's=pointModels{representativeIndices(j)};s.name=')
for line in ["    if j==1,assert(min(s.eigM)>0 && max(s.eigPaper)<0);end\n","    if j==2,assert(max(real(s.eigA))<0 && max(s.eigPaper)>0);end\n","    if j==3,assert(max(real(s.eigA))>0);end\n"]:s=s.replace(line,'')
s=s.replace("colors=[0,.35,.9;0,.55,.2;.85,.1,.1];markers={'o','^','x'};\nfig=figure", "colors=[.5,.25,.75;0,.35,.9;.9,.5,.05;0,.55,.2;.85,.1,.1];markers={'s','o','^','d','x'};\nfig=figure",1)
a=s.index('    if n>=2');b=s.index("title(tl,'Pourazad-style",a)
s=s[:a]+'end\n'+s[b:]
a=s.index('if available(1)');b=s.index("save_figure(fig,outDir,'06_",a)
s=s[:a]+'''for j=[2,4]
    q=cases{j};gi=linspace(q.gamma0-deg2rad(2),q.gamma0+deg2rad(2),1001);
    plot(ax,rad2deg(gi),frozen_lambda(gi,q,par),'Color',colors(j,:),'DisplayName',caseNames{j});
end
yline(ax,0,'k--');grid(ax,'on');legend(ax,'Location','best');
xlabel(ax,'Instantaneous gamma (deg)');ylabel(ax,'lambda max (Sym F)');
title(ax,'Certified references: nu=1 slice, not a trajectory certificate');
''' + s[b:]
s=s.replace("[60,60,1500,850]","[20,60,2400,900]")
s=s.replace('if j==3,phaseAngle=tr.gammaUnwrapped;end','phaseAngle=tr.gammaUnwrapped;')
s=s.replace("    if j==3\n        xlabel(ax,'gamma unwrapped (deg)');", "    plot(ax,rad2deg(phaseAngle(end)),tr.state(end,2)-1,'gd','MarkerSize',7,'LineWidth',1.5);\n    if true\n        xlabel(ax,'gamma unwrapped (deg)');")
s=s.replace("'Nonlinear PU trajectories: normalized analytical time / frequency convention'", "'Five operating points: common +1 deg perturbation; red circle = equilibrium, black square = initial, green diamond = final'")
s=s.replace("'boundaryCauses',boundaryTable,'nominal'", "'pointResults',pointResults,'pointModels',{pointModels},'representativeIndices',representativeIndices,'classificationTolerance',tol, ...\n    'boundaryCauses',boundaryTable,'nominal'")
s=s.replace('for j=1:3\n    if available(j)', 'for j=1:5\n    if available(j)')
s=s.replace("fprintf('Same MATLAB file updated;", "fprintf('\\n========================================================\\nFIVE REPRESENTATIVE NONLINEAR TRAJECTORIES\\n========================================================\\n');\nfor j=1:5,fprintf('%s: gamma0=%.9f deg; classification=%s\\n',caseNames{j},rad2deg(selected(j)),trajectories{j}.finalDiagnostics.interpretation);end\nfprintf('Same MATLAB file updated;")
a=s.index('% Case I: an invariant metric ball');b=s.index('xinit=[g0+delta;1];',a)
s=s[:a]+'''% Identical initial displacement for every representative. Certification
% may differ, but the disturbance is never reduced for a selected point.
g0=s.gamma0;delta=deg2rad(1);cert=[NaN,NaN];inv=cert;rate=NaN;certified=false;R=NaN;dw=NaN;
if s.contractionCertified
    % One common candidate rectangle; a bound is displayed only when its
    % invariant metric ball contains the unchanged +1-degree initial state.
    dg=min(deg2rad(15),.9*min(g0,pi/2-g0));dw=.1;
    [rate,R]=rectangle_certificate(s,p,dg,dw);
    certified=rate<0 && norm(s.Theta*[delta;0])<R;
    cert=g0+[-dg,dg];inv=cert;
end
if ~certified,fprintf('%s: no validated envelope for common +1 deg displacement.\\n',s.name);end
''' + s[b:]
s=s.replace("assert(err<1e-5);", "")
s=s.replace('if j==3\n    % The integrated angle','if true\n    % The integrated angle')
s=s.replace("'Case III angle sampling too coarse for unwrap.'","'Angle sampling too coarse for unwrap.'")
s=s.replace('    if abs(dnu)<=frequencyTol && abs(rp)<=powerTol','    if ~all(isfinite([finalDeg,finalNu,rp,err])) || err>=1e-5\n        interpretation=\'Numerically inconclusive\';\n    elseif abs(dnu)<=frequencyTol && abs(rp)<=powerTol')
s=s.replace("'returns to original equilibrium'","'Returns to original equilibrium'").replace("'approaches another equilibrium branch'","'Approaches another equilibrium branch'").replace("'still drifting and not converged'","'Still drifting / divergent'")
s=s.replace('CASE III FINAL-STATE DIAGNOSTICS','REPRESENTATIVE FINAL-STATE DIAGNOSTICS')
s=s.replace("    fprintf('Original equilibrium gamma0", "    fprintf('%s\\nInitial perturbed gamma = %.12f deg\\nL_final = %.12g\\n',s.name,rad2deg(xinit(1)),L(end));\n    fprintf('Original equilibrium gamma0")
s=s.replace('Pstar Case III','Pstar reference')
s=s.replace("'deltaGamma_deg',dg,'Pe_pu'", "'initialGamma_deg',rad2deg(xinit(1)),'L_final',L(end),'deltaGamma_deg',dg,'Pe_pu'")
s=s.replace("fprintf(fid,'Yörüngeler tam nonlinear PU denklemleriyle hesaplanır. Case I üstel sınırı, bir nonlinear açı-frekans dikdörtgeninde doğrulanan değişmez metrik topuna aittir. Başlangıç bozuntusu bu yerel sertifikaya uygun seçilir; geniş bozuntu testi olarak sunulmaz. Case II/III için azalan teorik sınır çizilmez.\\n\\n');", "fprintf(fid,'12 farklı denge noktası matrislerden sınıflandırıldı; tüm matrisler MAT dosyalarında, sınıflandırma 12_point_gfm_classification.csv dosyasında saklanır. Beş temsilci için ortak +1 derece açı bozuntusu ve nu(0)=1 kullanılır. Bunlar farklı bozuntu şiddetleri değildir. Üstel sınır yalnızca ortak başlangıcı içeren değişmez metrik topu doğrulanırsa gösterilir. L(t) özgün dengeye göredir.\\n\\n');\nfor j=find(r.caseAvailable),fprintf(fid,'- %s: %s; final gamma=%.9f derece; güç kalıntısı=%+.6g pu.\\n',r.cases{j}.name,r.trajectories{j}.finalDiagnostics.interpretation,r.trajectories{j}.finalDiagnostics.gammaUnwrapped_deg,r.trajectories{j}.finalDiagnostics.powerResidual_pu);end")
p.write_text(s,encoding='utf-8')
