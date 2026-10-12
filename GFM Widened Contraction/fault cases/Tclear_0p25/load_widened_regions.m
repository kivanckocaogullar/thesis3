function [Cpost,Cinactive]=load_widened_regions()
regionRoot=fullfile(getenv('USERPROFILE'),'Documents','MATLAB2','GFM Widened Contraction','contraction region');
a=load(fullfile(regionRoot,'metrics','active','metric.mat'),'bounds');i=load(fullfile(regionRoot,'metrics','inactive','metric.mat'),'bounds');
for mode={'active','inactive'},c=jsondecode(fileread(fullfile(regionRoot,'metrics',mode{1},'interval_verification.json')));assert(strcmp(c.status,'VERIFIED')&&c.unresolved_boxes==0&&c.counterexample_boxes==0);end
Cpost=[rad2deg(a.bounds(:,1))';a.bounds(:,2)'];Cinactive=[rad2deg(i.bounds(:,1))';i.bounds(:,2)'];
end
