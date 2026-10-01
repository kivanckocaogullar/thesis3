out=fileparts(mfilename('fullpath'));load(fullfile(out,'results.mat'),'p','n','cfg','runs','summary');audit=table();
for s=1:2,for k=1:5
 tr=runs{k,s};T=tr.data;pf=p;pf.U=p.U*cfg.retained(s);nf=rebuild(0,pf,n);
 powerError=0;timeError=0;currentEventResidual=0;continuity=NaN;
 for stage=1:2
  sg=tr.fault;pp=pf;nn=nf;if stage==2,sg=tr.post;pp=p;nn=n;end
  if isempty(sg)||~isfield(sg.sol,'x'),continue;end
  tt=linspace(sg.t(1)+1e-6,sg.t(end)-1e-6,101);[xx,dx]=deval(sg.sol,tt);
  timeError=max(timeError,max(abs(dx(1,:)-100*pi*(xx(2,:)-1))));
  d=coupled(xx(1,:),1,pp,nn,0);
  for j=1:numel(tt)
   node=nn.nodalY\[d.vc(j)/nn.Zt1;pp.U/(nn.Zgrid+nn.Zbreaker)];ig=(node(2)-pp.U)/(nn.Zgrid+nn.Zbreaker);vg=node(2)-nn.Zbreaker*ig;
   powerError=max(powerError,abs(real(vg*conj(ig))-d.P(j)));
  end
  for j=find(sg.sol.ie==1)
   d=coupled(sg.sol.ye(1,j),1,pp,nn,0);currentEventResidual=max(currentEventResidual,abs(d.rawI-1.2));
  end
 end
 if ~isempty(tr.post),continuity=max(abs(tr.fault.x(end,:)-tr.post.x(1,:)));end
 assert(powerError<1e-10&&timeError<2e-5&&currentEventResidual<1e-8&&continuity<1e-13);
 post=T.stage=="postfault";in=T.certifiedRegion>0&post;maxQ=max(T.lambdaMaxQc(in));minM=min(T.lambdaMinM(in));
 assert(maxQ<-1e-5&&minM>1e-5&&min(T.voltageMargin_pu)>-1e-8);
 en=T.energyNumerical(post);energyIncrease=max(diff(en));assert(energyIncrease<2e-7);
 field=100*pi*(tr.x(:,2)-1);assert(max(abs(T.frequencyDeviation_Hz-50*(tr.x(:,2)-1)))<1e-12);
 activeEvents=sum(tr.currentEvents.type=="activation");releaseEvents=sum(tr.currentEvents.type=="release");
 audit=[audit;table(cfg.caseNames(k),100*(1-cfg.retained(s)),powerError,timeError,currentEventResidual,continuity,maxQ,minM,energyIncrease,activeEvents,releaseEvents,'VariableNames',{'Case','sagDepth_percent','directNodalPowerError','physicalTimeDerivativeError','limiterEventResidual','clearingStateDiscontinuity','worstQcInside','minMetricInside','maxNumericalEnergyIncrease_postfault','activationEvents','releaseEvents'})]; %#ok<AGROW>
end,end
writetable(audit,fullfile(out,'model_consistency_audit.csv'));disp(audit);
% Summary counts describe sampled evidence, not rigorous flow containment.
summary.formalFaultRecoveryCertificate=repmat("Not established: numerical flow only",height(summary),1);
summary.firstDomainEntry_sampled_s=nan(height(summary),1);
for row=1:height(summary)
 k=find(cfg.caseNames==summary.Case(row));s=find(abs(100*(1-cfg.retained)-summary.sagDepth_percent(row))<1e-8);T=runs{k,s}.data;ix=find(T.stage=="postfault"&T.certifiedRegion>0,1);if ~isempty(ix),summary.firstDomainEntry_sampled_s(row)=T.time_s(ix);end
end
writetable(summary,fullfile(out,'comparison_summary.csv'));
save(fullfile(out,'audit_results.mat'),'audit','summary');
settings=rmfield(cfg,{'mi','ma','terminal','energyGamma','energyW'});settings.Pstar=1;settings.H=3.4;settings.K=20;settings.Kd=5;settings.Imax=1.2;settings.umax=1.3;settings.Xeq=imag(n.Zeq);settings.externalSCR=5;settings.externalXR=10;settings.LV_added_reactor_mH=.0858989549982;settings.sourceNominal_U=p.U;settings.physicalTime=true;settings.QV_approximation='instantaneous Eref=1-0.1*Qgrid';settings.energyMembership='numerical quadrature/interpolation, not a validated flow certificate';
fid=fopen(fullfile(out,'settings.json'),'w');fwrite(fid,jsonencode(settings),'char');fclose(fid);
function d=coupled(g,vs,p,n,mode)
 E=vs+zeros(size(g));
 for it=1:60
  d=raw(g,E,p,n,mode);h=E+.1*d.Q-vs;he=1+.1*d.QE;
  if max(abs(h),[],'all')<2e-14,break;end
  E=E-max(-.1,min(.1,h./he));
 end
 assert(all(E>0,'all')&&max(abs(h),[],'all')<1e-11&&all(abs(he)>1e-6,'all'));
 d.E=E;d.hE=he;d.Eg=-.1*d.Qg./he;d.totalPg=d.Pg+d.PE.*d.Eg;
end
function d=raw(g,E,p,n,mode)
 rot=exp(1i*(g+angle(n.Vth)));v=E.*rot;q=(v-n.Vth)/n.Zeq;qg=1i*v/n.Zeq;qE=rot/n.Zeq;r=abs(q);
 active=r>p.Imax;if mode==1,active=false(size(q));elseif mode==2,active=true(size(q));end
 ic=q;icg=qg;icE=qE;
 if any(active(:))
  ii=active;qa=q(ii);ra=r(ii);
  ic(ii)=p.Imax*qa./ra;
  icg(ii)=p.Imax*(qg(ii)./ra-qa.*real(conj(qa).*qg(ii))./ra.^3);
  icE(ii)=p.Imax*(qE(ii)./ra-qa.*real(conj(qa).*qE(ii))./ra.^3);
 end
 vc=n.Vth+n.Zeq*ic;vg=n.vgE*vc+n.vgU*p.U;ig=n.igE*vc+n.igU*p.U;
 vcg=n.Zeq*icg;vcE=n.Zeq*icE;
 S=vg.*conj(ig);Sg=n.vgE*vcg.*conj(ig)+vg.*conj(n.igE*vcg);SE=n.vgE*vcE.*conj(ig)+vg.*conj(n.igE*vcE);
 d=struct('P',real(S),'Q',imag(S),'Pg',real(Sg),'Qg',imag(Sg),'PE',real(SE),'QE',imag(SE),'vc',vc,'vg',vg,'ig',ig,'ic',ic,'u',vc+(n.Zconverter+p.AD)*ic,'rawI',r,'active',active);
end
function n=rebuild(x,p,n0)
n=n0;z1=n0.Zt1+1i*x;z2=n0.Zt2;zg=n0.Zgrid;zb=n0.Zbreaker;
y=[1/z1+1/z2+n0.Ym,-1/z2;-1/z2,1/z2+1/(zg+zb)];co=y\diag([1/z1,1/(zg+zb)]);
iv=co(2,1)/(zg+zb);iu=(co(2,2)-1)/(zg+zb);av=co(2,1)-zb*iv;au=co(2,2)-zb*iu;
icv=(1-co(1,1))/z1+n0.Yc;icu=-co(1,2)/z1;
n.Zt1=z1;n.nodalY=y;n.nodeCoefficients=co;n.vgE=av;n.vgU=au;n.igE=iv;n.igU=iu;n.icE=icv;n.icU=icu;n.Zeq=1/icv;n.Vth=-icu/icv*p.U;
end

function r=find_roots(g,y,f)
ix=find(y(1:end-1).*y(2:end)<0);r=g(abs(y)<1e-12);for j=ix',r(end+1,1)=fzero(f,g(j+[0 1]),optimset('TolX',1e-12));end;r=uniquetol(r,1e-9,'DataScale',1);
end





