%% Update only contraction membership in existing saved results; no ODE run.
root=fileparts(mfilename('fullpath'));addpath(root);file=fullfile(root,'fixed_Pstar_sag_sweep_results.mat');a=load(file);
regionRoot=fullfile(getenv('USERPROFILE'),'Documents','MATLAB2','GFM Widened Contraction','contraction region');w=load(fullfile(regionRoot,'inputs','model_inputs.mat'),'p','n');w.p.Pstar=1;w.p.omega=100*pi;assert(isequaln(a.p,w.p)&&isequaln(a.n,w.n));
[Cpost,Cinactive]=load_widened_regions();summary=a.summary;S=readtable(fullfile(root,'fixed_Pstar_sag_sweep.csv'),'TextType','string');assert(isequal(S.Uret,summary.Uret)&&max(abs(S.gamma_clear_deg-summary.gamma_clear_deg))<1e-10&&all(S.Tclear_s==a.cfg.tClear));
for k=1:height(summary)
 d=readtable(fullfile(root,sprintf('trajectory_Uret_%.2f.csv',summary.Uret(k))),'TextType','string');ix=find(d.stage=="postfault",1);assert(abs(d.time_s(ix)-a.cfg.tClear)<1e-12&&abs(d.gamma_deg(ix)-summary.gamma_clear_deg(k))<1e-9&&abs(d.nu(ix)-summary.nu_clear(k))<1e-12);
 active=d.raw_current_pu(ix)>a.p.Imax;assert(active==summary.limiter_active_postclear(k));g=summary.gamma_clear_deg(k);v=summary.nu_clear(k);
 summary.inside_Cpost(k)=active&&g>=Cpost(1,1)&&g<=Cpost(1,2)&&v>=Cpost(2,1)&&v<=Cpost(2,2);
 summary.inside_inactive_Cpost(k)=~active&&g>=Cinactive(1,1)&&g<=Cinactive(1,2)&&v>=Cinactive(2,1)&&v<=Cinactive(2,2);
 summary.inside_any_Cpost(k)=summary.inside_Cpost(k)||summary.inside_inactive_Cpost(k);
end
membership=summary(:,{'Uret','Tclear_s','gamma_clear_deg','nu_clear','limiter_active_postclear','inside_Cpost','inside_inactive_Cpost','inside_any_Cpost','classification'});
writetable(membership,fullfile(root,'widened_region_membership.csv'));writetable(summary,fullfile(root,'fixed_Pstar_sag_sweep.csv'));save(file,'summary','Cpost','Cinactive','-append');
check=load(file);assert(isequaln(check.cfg,a.cfg)&&isequaln(check.verification,a.verification));
for k=1:numel(a.runs),assert(isequaln(check.runs{k}.data,a.runs{k}.data)&&isequaln(check.runs{k}.clearState,a.runs{k}.clearState)&&isequaln(check.runs{k}.currentEvents,a.runs{k}.currentEvents));end
fid=fopen(fullfile(root,'fixed_Pstar_sag_sweep_summary.txt'),'w');guard=onCleanup(@()fclose(fid));
fprintf(fid,'Saved nominal sag sweep; Tclear=%.6g s.\n',a.cfg.tClear);
fprintf(fid,'Region membership updated from saved trajectories. See RUN_STATUS.txt for simulation provenance.\n');
fprintf(fid,'Widened active rectangle: gamma [%.10g %.10g] deg, nu [%.10g %.10g].\n',Cpost(1,:),Cpost(2,:));
fprintf(fid,'Widened inactive rectangle: gamma [%.10g %.10g] deg, nu [%.10g %.10g].\n',Cinactive(1,:),Cinactive(2,:));
fprintf(fid,'inside_Cpost means active-branch membership; inside_any_Cpost includes either correct-mode rectangle. Switching gap excluded. Membership is not a recovery guarantee.\n');
fprintf(fid,'Uret gamma_clear_deg nu_clear recovery inside_active inside_inactive inside_either\n');
for k=1:height(summary),s=summary(k,:);fprintf(fid,'%.2f %.12f %.12f %s %d %d %d\n',s.Uret,s.gamma_clear_deg,s.nu_clear,s.classification,s.inside_Cpost,s.inside_inactive_Cpost,s.inside_any_Cpost);end
clear guard;disp(membership);
plot_fixed_Pstar_nominal_sag_sweep_SCR5_XR10();
