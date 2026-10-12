function dom=hd_domain_limits(p,n,cfg)
%HD_DOMAIN_LIMITS Search-domain limits for the gamma searches (independent of H, D).
% inactiveFloor: lowest gamma below the equilibrium for which the inactive branch is
%   mode-consistent (|q| < Imax) and |u| < Vmax, limited to cfg.gammaFloorCap.
% activeCeiling: smallest of (a) the active-branch saddle equilibrium gamma_u where
%   Pe = Pstar beyond the current-limited peak (no contraction rectangle can contain a
%   non-attracting equilibrium), (b) loss of active-mode consistency or of |u| < Vmax,
%   (c) cfg.gammaCeilCap. A 0.05-deg margin is kept from (a) and (b).
gd=cfg.gammaScanDeg;g=deg2rad(gd);d=hd_coupled(g,1,p,n,0);P=d.P;
valid=~d.active&abs(d.u)<p.Vmax;eq=fzero(@(x)hd_coupled(x,1,p,n,0).P-p.Pstar,deg2rad([15 24]));
k0=find(gd<=rad2deg(eq),1,'last');k=k0;while k>1&&valid(k-1),k=k-1;end
dom.inactiveFloorReason="gamma scan cap";if k>1,dom.inactiveFloorReason="inactive mode / voltage validity";end
dom.inactiveFloor=max(deg2rad(gd(k)+.05),deg2rad(cfg.gammaFloorCap));if deg2rad(cfg.gammaFloorCap)>deg2rad(gd(k)+.05),dom.inactiveFloorReason="gamma floor cap";end
sw=fzero(@(x)hd_coupled(x,1,p,n,1).rawI-p.Imax,deg2rad([24 25]));
ka=find(gd>rad2deg(sw)+.06,1);va=d.active&abs(d.u)<p.Vmax;kk=ka;while kk<numel(gd)&&va(kk+1),kk=kk+1;end
ku=find(gd>rad2deg(sw)&P<p.Pstar,1);gu=fzero(@(x)hd_coupled(x,1,p,n,0).P-p.Pstar,g([ku-1 ku]));
cands=[rad2deg(gu)-.05,gd(kk)-.05,cfg.gammaCeilCap];[c,ic]=min(cands);
reasons=["active-branch saddle equilibrium gamma_u","active mode / voltage validity","gamma ceiling cap"];
dom.activeCeiling=deg2rad(c);dom.activeCeilingReason=reasons(ic);
dom.gamma_eq=eq;dom.gamma_switch=sw;dom.gamma_u=gu;dom.gamma_valid_inactive_low=deg2rad(gd(k));dom.gamma_valid_active_high=deg2rad(gd(kk));
end
