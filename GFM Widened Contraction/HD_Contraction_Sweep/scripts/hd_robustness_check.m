function T=hd_robustness_check(root)
%HD_ROBUSTNESS_CHECK Solver-sensitivity test at the non-monotone active nu_max points.
% Same rectangle (active gamma interval, nu in [0.9845, 1.05]) for H = 8 and
% D = 40..44, fitted with initial LMI grids 7, 9, 11 and 13. Everything else unchanged.
addpath(fullfile(root,'scripts'));cfg=hd_default_config(root);a=load(cfg.modelInputs,'p','n');p0=a.p;n=a.n;
g=cfg.gammaB{2};T=table();
for D=40:44
 p=hd_set_HD(p0,8,D-p0.K);
 for gN=[7 9 11 13]
  c=cfg;c.initGrid=gN;R=check_contraction_HD([g(1) .9845;g(2) 1.05],2,p,n,c);
  T=[T;table(8,D,gN,R.status,R.tmin,R.minM,R.maxQ,R.runtime_s,'VariableNames',{'H_s','D','init_grid','status','feasp_tmin','dense_minM','dense_maxQ','runtime_s'})]; %#ok<AGROW>
  fprintf('H=8 D=%d grid=%d: %s tmin=%.2e\n',D,gN,R.status,R.tmin);
 end
end
writetable(T,fullfile(root,'results','diagnostics','robustness_active_numax_H8.csv'));
end
