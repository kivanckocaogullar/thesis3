function S=hd_bisect(makeBounds,xAnchor,xLimit,tol,branch,p,n,cfg,anchorR,limitR)
%HD_BISECT Bracketed boundary search for one moving rectangle edge.
% makeBounds(x) returns the rectangle for edge value x. xAnchor: edge value of a
% rectangle expected feasible; xLimit: the other end of the search interval.
% Feasibility is monotone under domain inclusion for the exact problem (a metric valid
% on a rectangle is valid on any sub-rectangle, and the quadratic metric family is
% invariant under the affine coordinate normalization), so bisection is justified.
% The numerical procedure (finite adaptive grids, feasp accuracy) may violate this;
% all probes are stored, and inconclusive probes inside the bracket are flagged.
% Optional anchorR / limitR: already evaluated results for xAnchor / xLimit.
probes=struct([]);
    function R=probe(x)
        R=check_contraction_HD(makeBounds(x),branch,p,n,cfg);R.edge_value=x;
        if isempty(probes),probes=R;else,probes(end+1)=R;end
    end
    function tf=isinc(R),tf=startsWith(R.status,"inconclusive")||R.status=="solver_failed";end
if nargin>=9&&~isempty(anchorR),Ra=anchorR;Ra.edge_value=xAnchor;else,Ra=probe(xAnchor);end
S=struct('status',"",'bound',NaN,'limited_by',"",'best',[],'anchor',xAnchor,'limit',xLimit,'tol',tol, ...
 'bracket',[NaN NaN],'inconclusive_in_bracket',false,'n_probes',0,'probes',[]);
if Ra.status~="feasible"
 if isinc(Ra),S.status="inconclusive";else,S.status="no_feasible_rectangle";end
 S.limited_by="anchor rectangle "+Ra.status;S.best=Ra;S.probes=probes;S.n_probes=numel(probes);return
end
if nargin>=10&&~isempty(limitR),Rl=limitR;else,Rl=probe(xLimit);end
if Rl.status=="feasible"
 S.status="feasible";S.bound=xLimit;S.limited_by="search_domain_limit";S.best=Rl;S.bracket=[xLimit xLimit];
else
 good=xAnchor;bad=xLimit;best=Ra;incon=isinc(Rl);
 while abs(bad-good)>tol
  mid=(good+bad)/2;Rm=probe(mid);
  if Rm.status=="feasible",good=mid;best=Rm;else,bad=mid;incon=incon||isinc(Rm);end
 end
 S.status="feasible";S.bound=good;S.limited_by="contraction_feasibility";S.best=best;S.bracket=[good bad];S.inconclusive_in_bracket=incon;
end
S.probes=probes;S.n_probes=numel(probes);
end
