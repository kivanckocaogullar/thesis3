function run_fault_retention_0p80(stage,outputRoot)
% stage: 'all' (default), 'compute', 'plots', or 'postprocess'.
if nargin<1,stage='all';end
root=fileparts(fileparts(mfilename('fullpath')));
shared=fullfile(fileparts(root),'Fault_Retention_0p60','scripts');addpath(shared);
if nargin>1,root=outputRoot;for sub={'results','figures'},if ~isfolder(fullfile(root,sub{1})),mkdir(fullfile(root,sub{1}));end;end;end
if any(strcmp(stage,{'all','compute'}))&&~isfile(fullfile(root,'results','fct_analysis.mat')),fct_retention_engine('run',root,.80,[.10 .25 .50 .75 1.00]);end
if strcmp(stage,'postprocess'),fct_retention_engine('postprocess',root);end
if any(strcmp(stage,{'all','plots'}))
 plot_retention_FCT(fullfile(root,'results'),fullfile(root,'figures'));
 plot_retention_containment(fullfile(root,'results'),fullfile(root,'figures'));
 fct_retention_engine('currentplots',root);
end
end
