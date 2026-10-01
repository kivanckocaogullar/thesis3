root=fileparts(mfilename('fullpath'));repoRoot=fileparts(fileparts(root));diary(fullfile(root,'variant_build.txt'));
addpath(genpath(fullfile(repoRoot,'model','tdmsquared-master','tdmsquared-master','Library')));addpath(root);
old=fullfile(repoRoot,'model','GFMVSG1.slx');
efile=which('Subsys_Conv_Electrical_Part_SPS.slx');ep='GFM_Electrical_Xeq0350';mp='GFMVSG1_Xeq0350_audit';
assert(~isfile(fullfile(root,[ep '.slx']))&&~isfile(fullfile(root,[mp '.slx'])));
copyfile(efile,fullfile(root,[ep '.slx']));copyfile(old,fullfile(root,[mp '.slx']));
load_system(fullfile(root,[ep '.slx']));
filt=[ep '/LC_harmonic_filter'];trafo=[ep '/Transformer_Model'];pf=get_param(filt,'PortHandles');pt=get_param(trafo,'PortHandles');
react=[ep '/Additional_Interfacing_Reactor'];
add_block([filt '/converter_side_inductor'],react,'BranchType','L','Inductance','8.58989549982e-5','Position',[780 310 850 390]);
pr=get_param(react,'PortHandles');
for k=1:3
 lh=get_param(pf.RConn(k),'Line');assert(lh~=-1);delete_line(lh);
 add_line(ep,pf.RConn(k),pr.LConn(k),'autorouting','on');add_line(ep,pr.RConn(k),pt.LConn(k),'autorouting','on');
end
save_system(ep,fullfile(root,[ep '.slx']));
load_system(fullfile(root,[mp '.slx']));g=[mp '/Grid_Forming_Converter'];set_param(g,'LinkStatus','inactive');set_param([g '/Electrical_Part'],'ReferencedSubsystem',ep);
assert(strcmp(get_param([g '/Electrical_Part/Additional_Interfacing_Reactor'],'Inductance'),'8.58989549982e-5'));
save_system(mp,fullfile(root,[mp '.slx']));
fprintf('Created isolated variant %s and isolated electrical reference %s. Original files not saved.\n',mp,ep);
close_system(mp,0);close_system(ep,0);diary off;
