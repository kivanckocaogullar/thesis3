repoRoot=fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(genpath(fullfile(repoRoot,'model','tdmsquared-master','tdmsquared-master','Library')));
load_system(fullfile(repoRoot,'model','GFMVSG1.slx'));
g='GFMVSG1/Grid_Forming_Converter';
for s={g,[g '/Electrical_Part'],[g '/Electrical_Part/LC_harmonic_filter']}
 b=s{1};fprintf('BLOCK %s\n',b);for key={'LinkStatus','ReferenceBlock','ReferencedSubsystem'},try,disp([key{1} ': ' get_param(b,key{1})]);catch,end,end
 bs=find_system(b,'SearchDepth',1,'Type','Block');disp(bs);
end
fprintf('Initial model compile attempt (no simulation):\n');try,set_param('GFMVSG1','SimulationCommand','update');fprintf('COMPILE PASS\n');catch e,fprintf('COMPILE FAILURE: %s\n',getReport(e,'extended','hyperlinks','off'));end
close_system('GFMVSG1',0);
