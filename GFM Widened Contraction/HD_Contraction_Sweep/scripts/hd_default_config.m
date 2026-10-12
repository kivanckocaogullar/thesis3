function cfg=hd_default_config(root)
%HD_DEFAULT_CONFIG Settings of the H-Kd contraction sweep (single source of truth).
% root: HD_Contraction_Sweep folder. Baseline data are read, never written.
cfg.root=root;
cfg.baselineDir=fullfile(root,'..','GFM Contraction Analsys','limiter region','widened_nu_20261012_130655');
cfg.modelInputs=fullfile(root,'..','GFM Widened Contraction','contraction region','inputs','model_inputs.mat');
% Parameter grid (K fixed at the stored value 20; D = K + Kd).
cfg.H=[1 3.4 5.5 8];cfg.Kd=[0 5 15 30];cfg.tag='grid4x4_v2';
% Contraction method - unchanged from the baseline.
cfg.lambda=[1 0.1];            % [inactive active] s^-1, fixed
cfg.epsilon=1e-5;cfg.upper=10; % M > eps*I, M < 10*I, Q < -eps*I
cfg.fitMode=[1 0];             % mode flag used by the original fits (inactive forced / active auto)
cfg.initGrid=9;cfg.testGrid=101;cfg.denseGrid=301;cfg.maxRounds=8;cfg.addPerRound=30;
cfg.maxDenseRounds=4;          % documented extension: 301x301 violations fed back to the LMI grid
cfg.infeasTol=1e-7;            % feasp tmin in [0, infeasTol] -> inconclusive; > infeasTol -> infeasible
cfg.modeGridN=2001;            % 1-D gamma grid for limiter-mode / Vmax checks
% Baseline rectangles (exact edges from the saved, interval-verified metrics).
bi=load(fullfile(cfg.baselineDir,'inactive_metric.mat'),'bounds');ba=load(fullfile(cfg.baselineDir,'active_metric.mat'),'bounds');
cfg.baseBounds={bi.bounds,ba.bounds};
cfg.gammaFixedInactiveHi=bi.bounds(2,1);cfg.gammaFixedActiveLo=ba.bounds(1,1);
% Study A: nu band fixed, gamma edges searched from the baseline edges.
cfg.nuBandA=[.98 1.02];
cfg.gammaBaseInactive=bi.bounds(1,1);cfg.gammaBaseActive=ba.bounds(2,1);       % 10 deg, 43 deg
cfg.gammaAnchorInactive=deg2rad(20);cfg.gammaAnchorActive=deg2rad(26);           % smallest rectangles
cfg.tolGamma=deg2rad(.1);
% Study B: gamma intervals fixed at the baseline, nu edges searched.
cfg.gammaB={bi.bounds(:,1)',ba.bounds(:,1)'};
cfg.nuBaseHalfWidth=.02;cfg.nuCoreHalfWidth=1e-3;cfg.nuDomain=[.95 1.05];cfg.tolNu=5e-4;
% Gamma search-domain scan.
cfg.gammaScanDeg=-60:.01:100;cfg.gammaFloorCap=-60;cfg.gammaCeilCap=100;
end
