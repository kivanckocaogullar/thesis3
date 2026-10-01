repoRoot=fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(genpath(fullfile(repoRoot,'model','tdmsquared-master','tdmsquared-master','Library')));
load_system(fullfile(repoRoot,'model','GFMVSG1.slx'));
g='GFMVSG1/Grid_Forming_Converter';
fprintf('MAIN MASK\n'); names=get_param(g,'MaskNames');vals=get_param(g,'MaskValues');for k=1:numel(names),fprintf('%s = %s\n',names{k},vals{k});end
b=find_system(g,'LookUnderMasks','all','FollowLinks','on','Type','Block');
for k=1:numel(b)
 if contains(b{k},{'Cascaded_PI_controller','Operation_Mode_Selector','Active_Damping'})
 fprintf('\nBLOCK %s TYPE %s REF %s\n',b{k},get_param(b{k},'BlockType'),get_param(b{k},'ReferenceBlock'));
 dp=get_param(b{k},'DialogParameters');if isstruct(dp),fn=fieldnames(dp);for j=1:numel(fn),v=get_param(b{k},fn{j});if ischar(v),fprintf(' %s = %s\n',fn{j},v);end;end;end
 end
end
r=load('../baseline_SCR5_XR10/SCR5_XR10_baseline_results.mat','results');p=r.results.parameters;n=r.results.network;
gd=[rad2deg(p.gammaNominal);mean(r.results.regions_deg,2);17.457603;65];V=p.E*exp(1i*(deg2rad(gd)+p.phi));Ic=n.icE*V+n.icU*p.U;Ig=n.igE*V+n.igU*p.U;Vg=n.vgE*V+n.vgU*p.U;If=Ic-n.Yc*V;Vb=V+n.Zconverter*Ic;W=Vb+.2*Ic;
Ibase=p.S/(sqrt(3)*p.Vlv);tab=table(gd,real(Vg.*conj(Ig)),abs(Ic),1.2-abs(Ic),abs(Ic)*Ibase,abs(Ig),abs(n.Yc*V),abs(If),abs(Vb),abs(W),'VariableNames',{'gamma_deg','Pstar','Iconv_pu','currentMargin_pu','Iconv_RMS_A','Igrid_pu_HVbase','Ishunt_pu','IfilterGrid_pu_LVbase','Vbridge_required_pu','preActiveDamping_required_pu'});disp(tab);writetable(tab,'nominal_feasibility.csv');
fprintf('BASE %.12g RMS_A; PEAK %.12g A; Imax RMS %.12g A; PEAK %.12g A\n',Ibase,sqrt(2)*Ibase,1.2*Ibase,1.2*sqrt(2)*Ibase);

