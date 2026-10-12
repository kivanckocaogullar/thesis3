function A=search_gamma_bounds_HD(p,n,cfg,dom)
%SEARCH_GAMMA_BOUNDS_HD Study A: angle bounds with the nu band fixed (cfg.nuBandA).
% Inactive branch: upper edge fixed at the baseline inactive edge (switch - 0.05 deg);
%   lower edge searched from the baseline edge (10 deg): outward to dom.inactiveFloor if
%   the baseline rectangle is feasible, otherwise inward to cfg.gammaAnchorInactive.
% Active branch: lower edge fixed at the baseline active edge (switch + 0.05 deg);
%   upper edge searched from the baseline edge (43 deg): outward to dom.activeCeiling,
%   otherwise inward to cfg.gammaAnchorActive.
% The fixed edges come from the limiter transition, which is independent of H and D.
nu=cfg.nuBandA;gi=cfg.gammaFixedInactiveHi;ga=cfg.gammaFixedActiveLo;
A.inactive=hd_bound_search(@(x)[x nu(1);gi nu(2)],cfg.gammaBaseInactive,cfg.gammaAnchorInactive,dom.inactiveFloor,cfg.tolGamma,1,p,n,cfg);
A.active=hd_bound_search(@(x)[ga nu(1);x nu(2)],cfg.gammaBaseActive,cfg.gammaAnchorActive,dom.activeCeiling,cfg.tolGamma,2,p,n,cfg);
A.inactive.fixed_edge=gi;A.active.fixed_edge=ga;
end
