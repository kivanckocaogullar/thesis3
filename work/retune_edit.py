from pathlib import Path
p=Path('outputs/baseline_gfm_contraction.m')
s=p.read_text(encoding='utf-8')
s=s.replace('%% CURRENT-LIMITER-FREE VSM-GFM: PURELY INDUCTIVE SMIB BASELINE','%% CURRENT-LIMITER-FREE RETUNED VSM-GFM BASELINE\n% Candidate contraction-oriented retuning; no optimization is claimed.')
s=s.replace('H=0.1;K=20;Kd=40;D=K+Kd;','H=3.4;K=20;Kd=5;D=K+Kd;')
s=s.replace('assert(omega0_pu==1 && H==.1 && K==20 && Kd==40 && D==60 && b==300);','assert(omega0_pu==1 && H==3.4 && K==20 && Kd==5 && D==25);')
s=s.replace("fprintf('SI equivalents", "fprintf('Candidate contraction-oriented retuning: D/(2H)=%.12g; D/H=%.12g; Pmax=%.12g pu\\n',b,D/H,Pmax);\nfprintf('SI equivalents")
s=s.replace('alpha2=-D/(2*H)=-300','alpha2=-D/(2*H) (printed numerically below)')
s=s.replace("fprintf('Nominal Pstar=1 pu angle", "fprintf('alpha2=%.12g\\n',-b);\nfprintf('Nominal Pstar=1 pu angle")
a=s.index('selected=nan(1,3);');e=s.index('available=isfinite(selected)',a)
s=s[:a]+'''nominal=asin(Pstar_nominal/Pmax);snom=point_model(nominal,par);
nominalStable=max(real(snom.eigA))<0;
nominalCertified=min(snom.eigM)>0 && max(snom.eigPaper)<0;
fprintf('\\nNominal 1-pu operating point: gamma0=%.12f deg; residual=%.3e\\n',rad2deg(nominal),snom.equilibriumResidual);
for fld={'A','eigA','M','eigM','SymF_paper','eigPaper'}
    fprintf('%s=\\n',fld{1});disp(snom.(fld{1}));
end
fprintf('lambda_max(SymF)=%.12g; det(SymF)=%.12g\\n',max(snom.eigPaper),det(snom.SymF_paper));
yesno={'NO','YES'};
fprintf('LOCALLY STABLE = %s\\nCONTRACTION CERTIFIED = %s\\n',yesno{1+nominalStable},yesno{1+nominalCertified});
selected=nan(1,3);
if nominalCertified
    selected(1)=nominal;caseNames{1}='Case I - nominal 1 pu';
elseif ~isempty(intervals)
    [~,j]=max(diff(intervals,1,2));selected(1)=mean(intervals(j,:));
    fprintf('Retuning objective NOT achieved: nominal is not certified.\\n');
end
% Select the middle of the widest stable noncertified component, based on this scan.
mask=scan.maxRealEigA<0 & scan.lambdaPaper>1e-8;
starts=find(diff([false,mask,false])==1);stops=find(diff([false,mask,false])==-1)-1;
if ~isempty(starts)
    [~,j]=max(gamma0(stops)-gamma0(starts));selected(2)=gamma0(round((starts(j)+stops(j))/2));
end
unstable=find(scan.maxRealEigA>1e-8);
if ~isempty(unstable),selected(3)=gamma0(unstable(round((numel(unstable)+1)/2)));end
% Previously verified OLD PU baseline summary, never used in the new scan.
old=struct('H',.1,'K',20,'Kd',40,'D',60,'intervals_deg',[89.6730509036269,89.8749234883390], ...
    'nominalLambda',337.117533861893,'nominalCertified',false,'source','Previously verified OLD PU baseline; not recalculated');
comparison=struct('old',old,'newIntervals_deg',rad2deg(intervals),'newNominalLambda',max(snom.eigPaper), ...
    'newNominalCertified',nominalCertified,'widthChange_deg',sum(rad2deg(diff(intervals,1,2)))-diff(old.intervals_deg), ...
    'nominalLambdaChange',max(snom.eigPaper)-old.nominalLambda);
fprintf('\\n====================================================\\nCONTROLLER TUNING COMPARISON\\n====================================================\\n');
fprintf('OLD: H=0.1 K=20 Kd=40 D=60; nominal certified=NO\\n%s\\n',old.source);
print_intervals('OLD interval',deg2rad(old.intervals_deg));
fprintf('NEW: H=%.1f K=%.0f Kd=%.0f D=%.0f; nominal certified=%s\\n',H,K,Kd,D,yesno{1+nominalCertified});
print_intervals('NEW intervals',intervals);
fprintf('Total width: OLD %.12g deg; NEW %.12g deg; change %+.12g deg\\n',diff(old.intervals_deg),sum(rad2deg(diff(intervals,1,2))),comparison.widthChange_deg);
fprintf('Nominal lambda_max: OLD %.12g; NEW %.12g; change %+.12g\\n',old.nominalLambda,max(snom.eigPaper),comparison.nominalLambdaChange);
intervalTable=table();
for k=1:size(intervals,1)
    gd=rad2deg(intervals(k,:));pw=power_smib(intervals(k,:),par);
    row=table(k,gd(1),gd(2),pw(1),pw(2),diff(gd),boundaryTable.Cause(2*k-1),boundaryTable.Cause(2*k), ...
        'VariableNames',{'Interval','gamma_min_deg','gamma_max_deg','Pstar_min_edge','Pstar_max_edge','width_deg','lowerCause','upperCause'});
    intervalTable=[intervalTable;row]; %#ok<AGROW>
end
fprintf('Full interval summaries:\\n');disp(intervalTable);
writetable(intervalTable,fullfile(outDir,'contraction_intervals.csv'));
''' + s[e:]
s=s.replace("title(ax,'Positive reference metric; each gamma0 has a separate M');", "title(ax,'Positive reference metric; each gamma0 has a separate M');\nxline(ax,rad2deg(nominal),'m-.','Nominal 1 pu');")
s=s.replace("plot_consistency(ax,gamma0_deg,overview);xlim(ax,[0,120]);", "plot_consistency(ax,gamma0_deg,overview);xlim(ax,[0,120]);\nshade_intervals(ax,rad2deg(intervals));xline(ax,rad2deg(nominal),'m-.','Nominal 1 pu','HandleVisibility','off');")
s=s.replace("shade_intervals(ax,rad2deg(intervals));\n    for boundary", "shade_intervals(ax,rad2deg(intervals));\n    if nominal>=intervals(n,1) && nominal<=intervals(n,2),xline(ax,rad2deg(nominal),'m-.','Nominal 1 pu');end\n    for boundary")
s=s.replace("title(tl,'Contraction construction and derivative consistency');", "title(tl,'Pourazad-style per-unit VSM-GFM contraction model: retuned candidate');")
s=s.replace("title(ax,[names{n},': full scan, sign-preserving scale']);", "title(ax,[names{n},': full scan, sign-preserving scale']);\n    for g=rad2deg(intervals(:)).',xline(ax,g,'k:');end\n    shade_intervals(ax,rad2deg(intervals));")
s=s.replace("title(ax,'Purely inductive SMIB: mathematical operating powers are retained');legend(ax,'Location','best');", "title(ax,'Retuned VSM candidate: power and certified intervals');legend(ax,'Location','best');\nshade_intervals(ax,rad2deg(intervals));xline(ax,rad2deg(nominal),'m-.','Nominal 1 pu','HandleVisibility','off');")
s=s.replace("'status','POURAZAD_STYLE_PU_VSM_GFM_COMPLETE'", "'status','RETUNED_CURRENT_LIMITER_FREE_VSM_GFM_COMPLETE'")
s=s.replace("'boundaryCauses',boundaryTable);", "'boundaryCauses',boundaryTable,'nominal',snom,'nominalCertified',nominalCertified, ...\n    'nominalStable',nominalStable,'tuningComparison',comparison,'intervalTable',intervalTable);")
s=s.replace("fprintf('Same MATLAB file updated; outputs: %s\\n',outDir);diary off;", """fprintf('====================================================\\nRETUNED CURRENT-LIMITER-FREE VSM-GFM SUMMARY\\n====================================================\\n');
fprintf('Candidate contraction-oriented retuning\\nH=%.1f s; K=%.0f; Kd=%.0f; D_eff=%.0f\\n',H,K,Kd,D);
fprintf('Nominal Pstar=1 pu; gamma0=%.12f deg\\n',rad2deg(nominal));disp(intervalTable);
fprintf('Nominal contraction certified: %s\\nlambda_max(SymF): %.12g\\nLocal stability at nominal: %s\\n',yesno{1+nominalCertified},max(snom.eigPaper),yesno{1+nominalStable});
for j=1:3
    if available(j),fprintf('%s = %.12f deg\\n',caseNames{j},rad2deg(selected(j)));
    else,fprintf('%s: category unavailable\\n',caseNames{j});end
end
fprintf('Same MATLAB file updated; outputs: %s\\n',outDir);diary off;""")
s=s.replace('H=0.1, K=20, Kd=40, D=60;', 'H=3.4, K=20, Kd=5, D=25;')
s=s.replace('Dengede alpha2=-300.', 'Dengede alpha2=-25/6.8=-3.676470588235.')
s=s.replace("fprintf(fid,'## Yeni sonuçlar", "fprintf(fid,'Candidate contraction-oriented retuning; optimallik iddiası yoktur. D/(2H)=%.12g; D/H=%.12g.\\n\\n',r.parameters.D/r.parameters.mass,r.parameters.D/r.parameters.H);\nfprintf(fid,'## Yeni sonuçlar")
s=s.replace("fprintf(fid,'## Doğrulama ve yörüngeler", """fprintf(fid,'Nominal 1 pu sertifikası: %d; yerel kararlılık: %d; lambda_max(SymF)=%.12g.\\n\\n',r.nominalCertified,r.nominalStable,max(r.nominal.eigPaper));
fprintf(fid,'Önceden doğrulanmış eski PU ayarı: H=0.1, K=20, Kd=40, D=60; aralık [89.673050904,89.874923488] derece; nominal sertifikasız. Bu bilgi yalnızca karşılaştırmadır. Toplam aralık genişliği değişimi %+.9f derece; nominal lambda_max değişimi %+.9f.\\n\\n',r.tuningComparison.widthChange_deg,r.tuningComparison.nominalLambdaChange);
fprintf(fid,'## Doğrulama ve yörüngeler""")
p.write_text(s,encoding='utf-8')
