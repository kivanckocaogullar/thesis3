function validate_retention_results
base=fullfile(getenv('USERPROFILE'),'Documents','MATLAB2');
addpath(fullfile(base,'Fault_Retention_0p60','scripts'));
report=table();
for code={'60','80'}
 root=fullfile(base,['Fault_Retention_0p' code{1}]);out=fullfile(root,'results');a=load(fullfile(out,'fct_analysis.mat'));
 assert(numel(a.runs)==5&&all(a.results.horizon_s==8)&&~any(a.results.tighter));
 for k=1:5
  r=a.runs{k};tc=a.FCT(k);q=a.results(k,:);assert(r.fault.t(1)==0&&abs(r.fault.t(end)-tc)<1e-12&&abs(r.postfault.t(1)-tc)<1e-12&&r.postfault.t(end)==8);
  assert(all(isfinite([r.fault.x;r.postfault.x]),'all'));assert(max(abs(r.fault.x(end,:)-r.postfault.x(1,:)))<1e-12);
  dp=fct_retention_engine('electrical',r.postfault.x(:,1),a.p,a.n);tail=r.postfault.t>=7.2;
  recovered=max(abs(rad2deg(r.postfault.x(tail,1)-a.cfg.gammaEq)))<.1&&max(abs(r.postfault.x(tail,2)-1))<1e-4&&abs(dp.P(end)-1)<1e-4;
  assert(recovered==(q.classification=="RECOVERED"));
  ints=a.limiterIntervals(a.limiterIntervals.FCT_s==tc,:);
  assert(abs(sum(ints.Duration_s(ints.Stage=="FAULT"))-r.fault.activeDuration)<1e-10);
  assert(abs(sum(ints.Duration_s(ints.Stage=="POSTFAULT"))-r.postfault.activeDuration)<1e-10);
  assert(a.details.MaxActualCurrent_pu(k)<=a.p.Imax+1e-9&&a.details.MaxVoltageCommand_pu(k)<=a.p.Vmax+1e-9);
 end
 selected=struct('mass',a.p.mass,'D',a.p.D,'Imax',a.p.Imax,'Vmax',a.p.Vmax,'Pstar',a.p.Pstar,'omega',a.p.omega);
 config=struct('cfg',a.cfg,'real_model_parameters',selected,'initial_state',a.x0,'nominal_restored_voltage',a.p.U,'fault_voltage',a.pf.U,'switch_angle_deg',rad2deg(a.switchGamma));
 % Full complex network p/n is preserved losslessly in the MAT file.
 fid=fopen(fullfile(out,'numerical_configuration.json'),'w','n','UTF-8');fprintf(fid,'%s',jsonencode(config,'PrettyPrint',true));fclose(fid);
 report=[report;table(a.Uret,5,true,true,true,true,'VariableNames',{'Uret','Cases','Full8sHorizon','StateContinuity','RecoveryCriteriaUnchanged','CurrentVoltageLimitsPass'})];
end
ref=load(fullfile(base,'Fault_Retention_0p60','scripts','inputs','reference_040.mat'));
t=readtable(fullfile(base,'Fault_Retention_0p60','results','combined_040_060_080.csv'),'TextType','string');t=t(t.Uret==.4,:);
assert(max(abs(t.GammaClear_deg-ref.results.gamma_clear_deg))<1e-11&&max(abs(t.NuClear_pu-ref.results.nu_clear))<1e-13&&all(t.Recovery==ref.results.classification));
for code={'60','80'},writetable(report,fullfile(base,['Fault_Retention_0p' code{1}],'results','validation.csv'));end
disp(report);disp('PASS: 10 runs; original 0.40 results preserved; no CCT or metric calculation.');

end
