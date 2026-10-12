function S=hd_bound_search(makeBounds,xBase,xMin,xLimit,tol,branch,p,n,cfg)
%HD_BOUND_SEARCH Two-directional edge search starting from the baseline edge.
% xBase: baseline edge value (the baseline rectangle). xMin: edge of the smallest
% rectangle of the family. xLimit: search-domain limit.
%  - baseline rectangle feasible  -> bisection outward between xBase and xLimit;
%  - otherwise                    -> bisection inward between xMin and xBase.
% Starting from the baseline avoids very small anchor rectangles, for which the
% coordinate normalization makes Mdot large and the finite-grid fit less reliable.
R0=check_contraction_HD(makeBounds(xBase),branch,p,n,cfg);R0.edge_value=xBase;
if R0.status=="feasible"
 S=hd_bisect(makeBounds,xBase,xLimit,tol,branch,p,n,cfg,R0);S.direction="outward from baseline";
else
 S=hd_bisect(makeBounds,xMin,xBase,tol,branch,p,n,cfg,[],R0);S.direction="inward from baseline";
 if S.limited_by=="search_domain_limit",S.limited_by="contraction_feasibility";end % cannot happen (R0 infeasible); kept explicit
end
S.base_status=R0.status;
if isempty(S.probes),S.probes=R0;else,S.probes=[R0 S.probes];end
S.n_probes=numel(S.probes);
end
