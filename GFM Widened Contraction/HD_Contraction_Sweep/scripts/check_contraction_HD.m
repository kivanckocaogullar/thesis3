function R=check_contraction_HD(bounds,branch,p,n,cfg)
%CHECK_CONTRACTION_HD Fit and grid-check one candidate rectangle with the original method.
% bounds = [gamma_lo nu_lo; gamma_hi nu_hi] (gamma in rad). branch = 1 inactive, 2 active.
% Procedure (as in widen_nu_contraction.m / radial_*_metric scripts):
%   9x9 initial LMI grid -> feasp (accepted if tmin < 0, as in the original) ->
%   101x101 test -> add up to 30 worst nodes, max 8 rounds -> 301x301 dense check of
%   M > eps*I, M < upper*I, Q < -eps*I.
% Documented extension: if the 301x301 check fails, its worst violating nodes are added
% to the LMI grid and the procedure repeats, at most cfg.maxDenseRounds times. The
% acceptance criterion (301x301 with the eps margins) is unchanged.
% Added (cheap, gamma-only) checks: limiter-mode consistency and |u| < Vmax on a 1-D
% gamma grid; they do not depend on nu, H or D.
% Status: feasible | infeasible_lmi | mode_or_voltage_invalid | inconclusive_lmi |
%         inconclusive_refinement | inconclusive_dense | solver_failed
% NUMERICALLY GRID-CHECKED - NOT INTERVAL-CERTIFIED.
t0=tic;mode=cfg.fitMode(branch);lam=cfg.lambda(branch);ep=cfg.epsilon;up=cfg.upper;
center=mean(bounds,1);scale=diff(bounds,1)/2;
R=struct('status',"",'detail',"",'branch',branch,'bounds',bounds,'gamma_deg',rad2deg(bounds(:,1))','nu',bounds(:,2)', ...
 'lambda',lam,'epsilon',ep,'upper',up,'H',p.H,'D',p.D,'mass',p.mass,'C',[],'center',center,'scale',scale, ...
 'tmin',NaN,'rounds',0,'dense_rounds',0,'lmi_points',0,'minM',NaN,'maxM',NaN,'maxQ',NaN,'maxCond',NaN, ...
 'current_margin',NaN,'voltage_margin',NaN,'runtime_s',NaN,'error',"");
try
 % Limiter-mode and voltage-command validity (gamma only).
 g=linspace(bounds(1,1),bounds(2,1),cfg.modeGridN);d=hd_coupled(g,1,p,n,0);
 if branch==1,R.current_margin=min(p.Imax-d.rawI);modeOK=~any(d.active);
 else,R.current_margin=min(d.rawI-p.Imax);modeOK=all(d.active);end
 R.voltage_margin=min(p.Vmax-abs(d.u));
 if ~modeOK||R.current_margin<=0||R.voltage_margin<=0
  R.status="mode_or_voltage_invalid";R.detail=sprintf('current margin %.3g, voltage margin %.3g',R.current_margin,R.voltage_margin);R.runtime_s=toc(t0);return
 end
 [gg,vv]=ndgrid(linspace(-1,1,cfg.initGrid));grid=[gg(:),vv(:)];
 [gg,vv]=ndgrid(linspace(-1,1,cfg.testGrid));test=[gg(:),vv(:)];
 [gg,vv]=ndgrid(linspace(-1,1,cfg.denseGrid));dgrid=[gg(:),vv(:)];
 for dr=0:cfg.maxDenseRounds
  R.dense_rounds=dr;bad=1;
  for it=1:cfg.maxRounds
   R.rounds=R.rounds+1;R.lmi_points=size(grid,1);
   [C,tm]=hd_fit_metric(grid,lam,ep,up,center,scale,p,n,mode);R.tmin=tm;
   if ~(tm<0)
    if tm>cfg.infeasTol,R.status="infeasible_lmi";R.detail=sprintf('feasp tmin=%.3g > %.1g on %d LMI nodes',tm,cfg.infeasTol,size(grid,1));
    else,R.status="inconclusive_lmi";R.detail=sprintf('feasp tmin=%.3g in [0, %.1g]',tm,cfg.infeasTol);end
    R.runtime_s=toc(t0);return
   end
   z=hd_verify_metric(C,test,lam,center,scale,p,n,mode);
   scores=max([ep-z.minM,z.maxM-up,z.maxQ+ep],[],2);bad=find(scores>0);if isempty(bad),break;end
   [~,ix]=sort(scores,'descend');grid=unique([grid;test(ix(1:min(cfg.addPerRound,numel(bad))),:)],'rows');
  end
  R.C=C;
  if ~isempty(bad),R.status="inconclusive_refinement";R.detail=sprintf('%d of %d test nodes still violate after %d rounds',numel(bad),size(test,1),cfg.maxRounds);R.runtime_s=toc(t0);return,end
  dense=hd_verify_metric(C,dgrid,lam,center,scale,p,n,mode);
  R.minM=min(dense.minM);R.maxM=max(dense.maxM);R.maxQ=max(dense.maxQ);R.maxCond=max(dense.maxM./dense.minM);
  if R.minM>ep&&R.maxM<up&&R.maxQ<-ep,R.status="feasible";R.detail=sprintf('passes 301x301 dense check (%d dense feedback rounds)',dr);R.runtime_s=toc(t0);return,end
  sd=max([ep-dense.minM,dense.maxM-up,dense.maxQ+ep],[],2);bd=find(sd>0);[~,ix]=sort(sd(bd),'descend');
  grid=unique([grid;dgrid(bd(ix(1:min(cfg.addPerRound,numel(bd)))),:)],'rows');
 end
 R.status="inconclusive_dense";R.detail=sprintf('dense minM %.3g maxM %.3g maxQ %.3g after %d dense feedback rounds',R.minM,R.maxM,R.maxQ,cfg.maxDenseRounds);
catch e
 R.status="solver_failed";R.error=string(e.message);R.detail="exception";
end
R.runtime_s=toc(t0);
end
