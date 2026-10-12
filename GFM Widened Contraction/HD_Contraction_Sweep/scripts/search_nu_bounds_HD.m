function B=search_nu_bounds_HD(p,n,cfg,dom) %#ok<INUSD>
%SEARCH_NU_BOUNDS_HD Study B: frequency bounds with the gamma intervals fixed (cfg.gammaB).
% Reproducible asymmetric rule, per branch:
%  1) symmetric half-width delta*: largest [1-delta,1+delta] feasible, searched from the
%     baseline half-width cfg.nuBaseHalfWidth (0.02) outward to the domain limit, or
%     inward to cfg.nuCoreHalfWidth if the baseline band is not feasible;
%  2) lower extension: nu_max held at 1+delta*, nu_min searched down to cfg.nuDomain(1);
%  3) upper extension: nu_min held at 1-delta*, nu_max searched up to cfg.nuDomain(2);
%  4) combined rectangle [nu_min, nu_max] tested once and its status reported.
% Steps 2 and 3 both start from the symmetric core, so the result does not depend on
% the order of the one-sided searches.
B=struct();lab={'inactive','active'};
for b=1:2
 g=cfg.gammaB{b};S=struct();
 S.sym=hd_bound_search(@(dl)[g(1) 1-dl;g(2) 1+dl],cfg.nuBaseHalfWidth,cfg.nuCoreHalfWidth,min(1-cfg.nuDomain(1),cfg.nuDomain(2)-1),cfg.tolNu,b,p,n,cfg);
 if S.sym.status~="feasible"
  S.status=S.sym.status;S.nu_min=NaN;S.nu_max=NaN;S.delta=NaN;S.lower=[];S.upper=[];S.combined=[];
 else
  dl=S.sym.bound;S.delta=dl;core=S.sym.best;
  S.lower=hd_bisect(@(x)[g(1) x;g(2) 1+dl],1-dl,cfg.nuDomain(1),cfg.tolNu,b,p,n,cfg,core);
  S.upper=hd_bisect(@(x)[g(1) 1-dl;g(2) x],1+dl,cfg.nuDomain(2),cfg.tolNu,b,p,n,cfg,core);
  S.nu_min=S.lower.bound;S.nu_max=S.upper.bound;
  if S.nu_min==1-dl&&S.nu_max==1+dl,S.combined=core;
  else,S.combined=check_contraction_HD([g(1) S.nu_min;g(2) S.nu_max],b,p,n,cfg);end
  S.status="feasible";
 end
 B.(lab{b})=S;
end
end
