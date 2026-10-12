function out=run_HD_contraction_sweep(varargin)
%RUN_HD_CONTRACTION_SWEEP Reproduce the H-Kd contraction-bound sweep, tables and figures.
%   run_HD_contraction_sweep                               % 4x4 grid (tag 'grid4x4_v2')
%   run_HD_contraction_sweep('Kd',0:30,'tag','dense_Kd')   % user-specified grid
% Options: 'H', 'Kd', 'tag', 'workers', 'plots' (true/false), 'sensitivity' (true/false).
% K is fixed at the stored value (20); D = K + Kd. lambda, M < 10I, the metric
% structure and the limiter model are unchanged. Completed (H,Kd) points are
% checkpointed in results/diagnostics/checkpoints_<tag>/ and skipped on restart.
% All results are NUMERICALLY GRID-CHECKED - NOT INTERVAL-CERTIFIED.
root=fileparts(fileparts(mfilename('fullpath')));addpath(fullfile(root,'scripts'));
cfg=hd_default_config(root);
ip=inputParser;ip.addParameter('H',cfg.H);ip.addParameter('Kd',cfg.Kd);ip.addParameter('tag',cfg.tag);
ip.addParameter('workers',4);ip.addParameter('plots',true);ip.addParameter('sensitivity',true);ip.parse(varargin{:});o=ip.Results;
cfg.H=o.H(:)';cfg.Kd=o.Kd(:)';cfg.tag=o.tag;
a=load(cfg.modelInputs,'p','n');p0=a.p;n=a.n;cfg.K=p0.K;cfg.D=cfg.K+cfg.Kd;
dres=fullfile(root,'results','diagnostics');ck=fullfile(dres,['checkpoints_' cfg.tag]);if ~isfolder(ck),mkdir(ck);end
logf=fullfile(dres,['run_log_' cfg.tag '.txt']);diary(logf);cleanup=onCleanup(@()diary('off'));
fprintf('=== H-Kd contraction sweep [%s] %s\n',cfg.tag,char(datetime('now')));
%% Phase 1: baseline consistency (saved, interval-verified metrics; read only)
dom=hd_domain_limits(p0,n,cfg);base=baseline_check(cfg,p0,n);
writetable(struct2table(base),fullfile(dres,'baseline_check.csv'));
dl=struct('inactive_floor_deg',rad2deg(dom.inactiveFloor),'inactive_floor_reason',dom.inactiveFloorReason, ...
 'active_ceiling_deg',rad2deg(dom.activeCeiling),'active_ceiling_reason',dom.activeCeilingReason,'gamma_eq_deg',rad2deg(dom.gamma_eq), ...
 'gamma_switch_deg',rad2deg(dom.gamma_switch),'gamma_u_deg',rad2deg(dom.gamma_u),'note','Independent of H and Kd (electrical closure does not use m or D).');
fid=fopen(fullfile(dres,'domain_limits.json'),'w');fwrite(fid,jsonencode(dl,PrettyPrint=true),'char');fclose(fid);
%% Phases 2-3: sweep with checkpointing
[HH,KK]=ndgrid(cfg.H,cfg.Kd);N=numel(HH);sig=config_signature(cfg);
files=arrayfun(@(k)fullfile(ck,sprintf('H%s_Kd%s.mat',num2tag(HH(k)),num2tag(KK(k)))),1:N,'UniformOutput',false);
todo=find(~cellfun(@(f)valid_checkpoint(f,sig),files));fprintf('%d of %d (H,Kd) points to compute\n',numel(todo),N);
if ~isempty(todo)
 nw=o.workers;if nw>1&&isempty(gcp('nocreate')),parpool('local',nw);end
 Hs=HH(todo);Ks=KK(todo);fs=files(todo);
 parfor (j=1:numel(todo),nw*(nw>1))
  res=compute_point(Hs(j),Ks(j),p0,n,cfg,dom);res.signature=sig;
  tmp=[fs{j}(1:end-4) '_part.mat'];parsave(tmp,res);movefile(tmp,fs{j},'f');
  fprintf('done H=%g Kd=%g (%.0f s)\n',Hs(j),Ks(j),res.runtime_s);
 end
end
for k=N:-1:1,t=load(files{k},'res');res(k)=t.res;end;res=reshape(res,size(HH));
%% Initialization sensitivity at the baseline (Study A only)
sens=table();if o.sensitivity,sens=init_sensitivity(cfg,p0,n,dom,dres);end
%% Tables, JSON exports, MAT
out=assemble(res,cfg,dom,base,sens,p0,n,root);
if o.plots
 plot_HD_gamma_surfaces(root,cfg.tag);plot_HD_nu_surfaces(root,cfg.tag);plot_HD_feasibility_maps(root,cfg.tag);plot_HD_region_comparison(root,cfg.tag);plot_HD_line_profiles(root,cfg.tag);
end
fprintf('=== finished %s\n',char(datetime('now')));
end

function res=compute_point(H,Kd,p0,n,cfg,dom)
t0=tic;p=hd_set_HD(p0,H,Kd);res.H=H;res.Kd=Kd;res.K=p.K;res.D=p.D;res.mass=p.mass;
res.A=search_gamma_bounds_HD(p,n,cfg,dom);res.B=search_nu_bounds_HD(p,n,cfg,dom);
res.runtime_s=toc(t0);res.finished=char(datetime('now'));res.signature='';
end
function parsave(f,res),save(f,'res'),end
function s=num2tag(x),s=strrep(sprintf('%g',x),'.','p');end
function s=config_signature(cfg)
s=sprintf('v2|lam%g_%g|eps%g|up%g|grids%d_%d_%d_%d_%d_%d|tol%g|nuA%g_%g|tolG%g|tolN%g|nuDom%g_%g|core%g|nubase%g|anch%g_%g|base%g_%g|K%g', ...
 cfg.lambda,cfg.epsilon,cfg.upper,cfg.initGrid,cfg.testGrid,cfg.denseGrid,cfg.maxRounds,cfg.addPerRound,cfg.maxDenseRounds,cfg.infeasTol,cfg.nuBandA, ...
 cfg.tolGamma,cfg.tolNu,cfg.nuDomain,cfg.nuCoreHalfWidth,cfg.nuBaseHalfWidth,cfg.gammaAnchorInactive,cfg.gammaAnchorActive,cfg.gammaBaseInactive,cfg.gammaBaseActive,cfg.K);
end
function ok=valid_checkpoint(f,sig)
ok=false;if ~isfile(f),return,end
try,r=load(f);ok=isfield(r,'res')&&isfield(r.res,'signature')&&strcmp(r.res.signature,sig);catch,ok=false;end
end

function base=baseline_check(cfg,p0,n)
% Saved widened metrics evaluated with the copied verify function at the stored H, D.
lab={'inactive','active'};base=struct([]);
for b=1:2
 m=load(fullfile(cfg.baselineDir,[lab{b} '_metric.mat']));
 [g,v]=ndgrid(linspace(-1,1,cfg.denseGrid));z=hd_verify_metric(m.C,[g(:),v(:)],m.lambda,m.center,m.scale,p0,n,cfg.fitMode(b));
 cert=jsondecode(fileread(fullfile(cfg.baselineDir,'interval_verification_20261012',lab{b},'interval_verification.json')));
 r=struct('branch',string(lab{b}),'H',p0.H,'K',p0.K,'Kd',p0.Kd,'D',p0.D,'lambda',m.lambda,'gamma_lo_deg',rad2deg(m.bounds(1,1)),'gamma_hi_deg',rad2deg(m.bounds(2,1)), ...
  'nu_lo',m.bounds(1,2),'nu_hi',m.bounds(2,2),'dense_minM',min(z.minM),'dense_maxM',max(z.maxM),'dense_maxQ',max(z.maxQ), ...
  'saved_dense_minM',min(m.dense.minM),'saved_dense_maxQ',max(m.dense.maxQ), ...
  'reproduced',abs(min(z.minM)-min(m.dense.minM))<1e-12&&abs(max(z.maxQ)-max(m.dense.maxQ))<1e-12,'interval_status',string(cert.status));
 if isempty(base),base=r;else,base(end+1)=r;end %#ok<AGROW>
 fprintf('baseline %s: dense minM %.6g maxQ %.6g reproduced=%d interval=%s\n',lab{b},r.dense_minM,r.dense_maxQ,r.reproduced,r.interval_status);
end
end

function T=init_sensitivity(cfg,p0,n,dom,dres)
% Study A at the baseline (H, Kd) with different initial LMI grids.
f=fullfile(dres,['init_sensitivity_' cfg.tag '.csv']);if isfile(f),T=readtable(f,'TextType','string');return,end
p=hd_set_HD(p0,p0.H,p0.Kd);grids=[7 9 11];T=table();
for gN=grids
 c=cfg;c.initGrid=gN;A=search_gamma_bounds_HD(p,n,c,dom);
 T=[T;table(gN,rad2deg(A.inactive.bound),A.inactive.status,rad2deg(A.active.bound),A.active.status, ...
  'VariableNames',{'init_grid','gamma_min_inactive_deg','status_inactive','gamma_max_active_deg','status_active'})]; %#ok<AGROW>
end
writetable(T,f);disp(T);
end

function out=assemble(res,cfg,dom,base,sens,p0,n,root)
[nH,nK]=size(res);lab={'inactive','active'};
G=struct();U=struct();SA=strings(nH,nK,2);SB=strings(nH,nK,2);probes=table();rows=table();
for b=1:2
 G.(['min_' lab{b}])=nan(nH,nK);G.(['max_' lab{b}])=nan(nH,nK);G.(['limited_' lab{b}])=strings(nH,nK);
 U.(['min_' lab{b}])=nan(nH,nK);U.(['max_' lab{b}])=nan(nH,nK);U.(['delta_' lab{b}])=nan(nH,nK);
 U.(['combined_' lab{b}])=strings(nH,nK);U.(['limited_min_' lab{b}])=strings(nH,nK);U.(['limited_max_' lab{b}])=strings(nH,nK);
end
for i=1:nH
 for k=1:nK
  r=res(i,k);p=hd_set_HD(p0,r.H,r.Kd);ctag=sprintf('H%s_Kd%s',strrep(sprintf('%g',r.H),'.','p'),strrep(sprintf('%g',r.Kd),'.','p'));
  for b=1:2
   A=r.A.(lab{b});SA(i,k,b)=A.status;G.(['limited_' lab{b}])(i,k)=A.limited_by;
   if A.status=="feasible"
    if b==1,G.min_inactive(i,k)=rad2deg(A.bound);G.max_inactive(i,k)=rad2deg(A.fixed_edge);
    else,G.min_active(i,k)=rad2deg(A.fixed_edge);G.max_active(i,k)=rad2deg(A.bound);end
    hd_export_metric_json(A.best,p,n,fullfile(root,'results','angle_bounds',['metrics_json_' cfg.tag],ctag,lab{b},'metric.json'));
   end
   Bb=r.B.(lab{b});SB(i,k,b)=Bb.status;
   if Bb.status=="feasible"
    U.(['min_' lab{b}])(i,k)=Bb.nu_min;U.(['max_' lab{b}])(i,k)=Bb.nu_max;U.(['delta_' lab{b}])(i,k)=Bb.delta;
    U.(['combined_' lab{b}])(i,k)=Bb.combined.status;U.(['limited_min_' lab{b}])(i,k)=Bb.lower.limited_by;U.(['limited_max_' lab{b}])(i,k)=Bb.upper.limited_by;
    Rb=Bb.combined;if Rb.status~="feasible",Rb=Bb.sym.best;end
    hd_export_metric_json(Rb,p,n,fullfile(root,'results','frequency_bounds',['metrics_json_' cfg.tag],ctag,lab{b},'metric.json'));
   end
   probes=[probes;probe_rows(r,b,'A',A.probes);probe_rows(r,b,'B_sym',Bb.sym.probes)]; %#ok<AGROW>
   if Bb.status=="feasible"
    probes=[probes;probe_rows(r,b,'B_lower',Bb.lower.probes);probe_rows(r,b,'B_upper',Bb.upper.probes)]; %#ok<AGROW>
    if ~(Bb.nu_min==1-Bb.delta&&Bb.nu_max==1+Bb.delta),probes=[probes;probe_rows(r,b,'B_combined',Bb.combined)];end %#ok<AGROW>
   end
  end
  Ai=r.A.inactive.best;Aa=r.A.active.best;
  rows=[rows;table(r.H,r.K,r.Kd,r.D,SA(i,k,1),SA(i,k,2),SB(i,k,1),SB(i,k,2),G.min_inactive(i,k),G.max_inactive(i,k),G.min_active(i,k),G.max_active(i,k), ...
   U.min_inactive(i,k),U.max_inactive(i,k),U.combined_inactive(i,k),U.min_active(i,k),U.max_active(i,k),U.combined_active(i,k), ...
   Ai.minM,Ai.maxQ,Aa.minM,Aa.maxQ,r.runtime_s,warn_text(r), ...
   'VariableNames',{'H_s','K','Kd','D','A_status_inactive','A_status_active','B_status_inactive','B_status_active', ...
   'gamma_min_inactive_deg','gamma_max_inactive_deg_fixed','gamma_min_active_deg_fixed','gamma_max_active_deg', ...
   'nu_min_inactive','nu_max_inactive','B_combined_inactive','nu_min_active','nu_max_active','B_combined_active', ...
   'A_inactive_minM','A_inactive_maxQ','A_active_minM','A_active_maxQ','runtime_s','warnings'})]; %#ok<AGROW>
 end
end
tag=cfg.tag;dres=fullfile(root,'results','diagnostics');
writetable(rows,fullfile(dres,['summary_table_' tag '.csv']));writetable(probes,fullfile(dres,['all_probes_' tag '.csv']));
angle=rows(:,{'H_s','K','Kd','D','A_status_inactive','A_status_active','gamma_min_inactive_deg','gamma_max_inactive_deg_fixed', ...
 'gamma_min_active_deg_fixed','gamma_max_active_deg','A_inactive_minM','A_inactive_maxQ','A_active_minM','A_active_maxQ'});
freq=rows(:,{'H_s','K','Kd','D','B_status_inactive','B_status_active','nu_min_inactive','nu_max_inactive','B_combined_inactive', ...
 'nu_min_active','nu_max_active','B_combined_active'});
writetable(angle,fullfile(root,'results','angle_bounds',['angle_bounds_' tag '.csv']));
writetable(freq,fullfile(root,'results','frequency_bounds',['frequency_bounds_' tag '.csv']));
Hgrid=cfg.H;Kdgrid=cfg.Kd;Dgrid=cfg.D;
save(fullfile(root,'results','angle_bounds',['angle_bounds_' tag '.mat']),'Hgrid','Kdgrid','Dgrid','G','SA','cfg','dom');
save(fullfile(root,'results','frequency_bounds',['frequency_bounds_' tag '.mat']),'Hgrid','Kdgrid','Dgrid','U','SB','cfg','dom');
save(fullfile(dres,['sweep_results_' tag '.mat']),'res','cfg','dom','base','sens','G','U','SA','SB','Hgrid','Kdgrid','Dgrid','-v7.3');
out=struct('rows',rows,'G',G,'U',U,'SA',SA,'SB',SB);disp(rows(:,1:8));
end
function T=probe_rows(r,b,study,P)
T=table();if isempty(P),return,end
for q=1:numel(P)
 x=P(q);ev=NaN;if isfield(x,'edge_value'),ev=x.edge_value;end
 T=[T;table(r.H,r.Kd,r.D,b,string(study),x.status,string(x.detail),x.gamma_deg(1),x.gamma_deg(2),x.nu(1),x.nu(2),ev,x.lambda,x.tmin, ...
  x.rounds,x.lmi_points,x.minM,x.maxM,x.maxQ,x.maxCond,x.current_margin,x.voltage_margin,x.runtime_s,string(x.error), ...
  'VariableNames',{'H_s','Kd','D','branch','search','status','detail','gamma_lo_deg','gamma_hi_deg','nu_lo','nu_hi','edge_value','lambda', ...
  'feasp_tmin','refine_rounds','lmi_points','minM','maxM','maxQ','maxCondM','current_margin_pu','voltage_margin_pu','runtime_s','error'})]; %#ok<AGROW>
end
end
function w=warn_text(r)
w="";lab={'inactive','active'};
for b=1:2
 A=r.A.(lab{b});if A.inconclusive_in_bracket,w=w+"A-"+lab{b}+": inconclusive probe in bracket; ";end
 Bb=r.B.(lab{b});
 if Bb.status=="feasible"
  if Bb.lower.inconclusive_in_bracket||Bb.upper.inconclusive_in_bracket||Bb.sym.inconclusive_in_bracket,w=w+"B-"+lab{b}+": inconclusive probe in bracket; ";end
  if Bb.combined.status~="feasible",w=w+"B-"+lab{b}+": one-sided nu extents not jointly feasible ("+Bb.combined.status+"); ";end
 end
end
end
