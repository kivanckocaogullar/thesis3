root=fileparts(mfilename('fullpath'));out=fullfile(root,'radial_detailed_controller_audit');assert(~isfolder(out));mkdir(out);diary(fullfile(out,'model_inspection.txt'));
repoRoot=fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(repoRoot,'model','tdmsquared-master','tdmsquared-master','Library')));
load_system(fullfile(repoRoot,'model','GFMVSG1.slx'));m='GFMVSG1';
fprintf('Solver: %s, type %s, fixedstep %s\n',get_param(m,'Solver'),get_param(m,'SolverType'),get_param(m,'FixedStep'));
g=[m '/Grid_Forming_Converter'];names=get_param(g,'MaskNames');vals=get_param(g,'MaskValues');for k=1:numel(names),fprintf('MASK %s = %s\n',names{k},vals{k});end
b=find_system(m,'LookUnderMasks','all','FollowLinks','on','Type','Block');
for k=1:numel(b)
 path=b{k};
 if contains(path,{'Virtual_Synchronous_Machine','PLL','measurement_transformations','Three-Phase Series RLC','PI_freeze','Integrator_Freeze_Logic'})||strcmp(get_param(path,'MaskType'),'Three-Phase Series RLC Branch')
  fprintf('\nBLOCK %s TYPE %s REF %s\n',path,get_param(path,'BlockType'),get_param(path,'ReferenceBlock'));
  dp=get_param(path,'DialogParameters');if isstruct(dp),fn=fieldnames(dp);for j=1:numel(fn),v=get_param(path,fn{j});if ischar(v),fprintf(' %s = %s\n',fn{j},v);end,end,end
 end
end
fprintf('\nTOOLS\n');disp(which('findop'));disp(which('linearize'));disp(which('power_analyze'));disp(which('power_initstates'));
close_system(m,0);diary off;
