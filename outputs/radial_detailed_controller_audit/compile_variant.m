root=fileparts(mfilename('fullpath'));repoRoot=fileparts(fileparts(root));diary(fullfile(root,'variant_compile.txt'));addpath(root);addpath(genpath(fullfile(repoRoot,'model','tdmsquared-master','tdmsquared-master','Library')));
m='GFMVSG1_Xeq0350_audit';load_system(fullfile(root,[m '.slx']));
fprintf('Loaded %s at %s\n',m,char(datetime('now')));
watch=timer('ExecutionMode','fixedSpacing','Period',20,'TimerFcn',@watchdog);start(watch);
try
 set_param(m,'SimulationCommand','update');fprintf('COMPILE SUCCESS at %s\n',char(datetime('now')));
 op=operpoint(m);fprintf('Operating point contains %d state blocks.\n',numel(op.States));save(fullfile(root,'untrimmed_state_inventory.mat'),'op');
catch e
 fprintf('COMPILE ERROR\n%s\n',getReport(e,'extended','hyperlinks','off'));
end
stop(watch);delete(watch);close_system(m,0);diary off;
function watchdog(~,~)
fprintf('WATCH %s\n',char(datetime('now')));ff=findall(groot,'Type','figure');for f=ff',fprintf('FIGURE %s\n',get(f,'Name'));end
end

