%% Three real limiter strategies; ideal raw-reference proxy, not full PI EMT.
% Inspected source: Library/Shared Functions/dq_Limiter.slx
% axis_prioritization: system_102 radial, system_121 d, system_152 q.
% The norm comparator is > Imax; priority-axis comparator is abs(axis)>Imax.
% Zero-based switch index = normExceeded + priorityAxisExceeded:
% 0: raw r; 1: preserve priority axis and sign(other)*sqrt(Imax^2-axis^2);
% 2: sign(priority)*Imax on priority axis and zero on the other.
% Thus d priority: sd=clip(rd,+/-I), sq=clip(rq,+/-sqrt(I^2-sd^2));
% q priority: sq=clip(rq,+/-I), sd=clip(rd,+/-sqrt(I^2-sq^2)).
% For radial: s=min(1,I/abs(r))*r. Equality has inactive flag.
%
% APPROXIMATION retained: q=(Vtarget-Vth)/Zeq proxies the raw PI reference.
% r_d+j*r_q=q*exp(-j*theta), theta=gamma+arg(Vth), is in the VSM frame.
% Apply limiter in that frame, then Iconv=exp(j*theta)*(sd+j*sq).
% Ideal current tracking; PI/freeze states omitted. Reconstruct Vc=Vth+Zeq*I
% and grid measurements for every branch. Voltage-command limit is checked,
% NOT enforced by a second ideal limiter. No metric, figures, or output files.
root=fileparts(mfilename('fullpath'));
a=load(fullfile(root,'current_limiter_SCR5_XR10','current_limiter_SCR5_XR10_results.mat'),'results');p=a.results.parameters;n0=a.results.network;
assert(p.Imax==1.2&&p.Vmax==1.3&&p.AD==.2&&p.H==3.4&&p.K==20&&p.Kd==5&&p.D==25);
x=fzero(@(z)imag(rebuild(z,p,n0).Zeq)-.350,[0 .6],optimset('TolX',1e-13));n=rebuild(x,p,n0);
assert(abs(imag(n.Zeq)-.350)<1e-12&&n.Zgrid==n0.Zgrid);
sw=fzero(@(g)abs(map(deg2rad(g),p,n,1,0).q)-p.Imax,[0 180],optimset('TolX',1e-12));
names=["Constant ratio d/q","d priority","q priority"];results=cell(3,1);rows=cell(3,1);
for strategy=1:3
 left=map(deg2rad(sw),p,n,strategy,1);right=map(deg2rad(sw),p,n,strategy,2);
 assert(abs(left.Pe-right.Pe)<1e-11&&abs(left.vc-right.vc)<1e-11);
 % Split at all priority component thresholds BEFORE root refinement.
 gg=linspace(sw,180,180001)';raw=map(deg2rad(gg),p,n,1,0).raw;
 cuts=[sw;180];
 if strategy>1
  axis=real(raw);if strategy==3,axis=imag(raw);end
  for threshold=[-p.Imax,p.Imax]
   cuts=[cuts;roots_grid(gg,axis-threshold,@(g)priority_axis(g,p,n,strategy)-threshold)]; %#ok<AGROW>
  end
 end
 cuts=uniquetol(sort(cuts),1e-8,'DataScale',1);stable=[];saddles=[];voltage=[];usable=[];exportUsable=[];slopeRoots=[];
 maxU=0;minU=Inf;maxI=0;validation=struct('PeDerivative',0,'SecondDerivative',0,'KCL',0);
 for j=1:numel(cuts)-1
  lo=cuts(j);hi=cuts(j+1);epsg=min(1e-7,(hi-lo)/1000);g=linspace(lo+epsg,hi-epsg,max(1001,ceil((hi-lo)/.001)))';d=map(deg2rad(g),p,n,strategy,2);
  % Every smooth slope extremum is a potential bracket endpoint.
  ex=roots_grid(g,d.Pgg,@(z)map(deg2rad(z),p,n,strategy,2).Pgg);stops=sort([g(1);ex;g(end)]);sv=map(deg2rad(stops),p,n,strategy,2).Pg;
  zr=roots_grid(stops,sv,@(z)map(deg2rad(z),p,n,strategy,2).Pg);slopeRoots=[slopeRoots;zr];
  uv=roots_grid(g,abs(d.u)-p.Vmax,@(z)abs(map(deg2rad(z),p,n,strategy,2).u)-p.Vmax);
  pr=roots_grid(g,d.Pe,@(z)map(deg2rad(z),p,n,strategy,2).Pe);
  edges=uniquetol(sort([lo;hi;zr;uv;pr]),1e-8,'DataScale',1);
  for q=1:numel(edges)-1
   v=map(deg2rad(mean(edges(q:q+1))),p,n,strategy,2);piece=edges(q:q+1)';
   if v.Pg>0,stable=[stable;piece];else,saddles=[saddles;piece];end %#ok<AGROW>
   if abs(v.u)<=p.Vmax,voltage=[voltage;piece];end %#ok<AGROW>
   if v.Pg>0&&abs(v.u)<=p.Vmax,usable=[usable;piece];if v.Pe>0,exportUsable=[exportUsable;piece];end;end %#ok<AGROW>
  end
  % Refine extrema of required voltage magnitude as well as sampling it.
  ud=real(conj(d.u).*d.ug);ur=roots_grid(g,ud,@(z)voltage_slope(z,p,n,strategy));ux=map(deg2rad([lo;hi;ur]),p,n,strategy,2);
  maxU=max([maxU;abs(ux.u);abs(d.u)]);minU=min([minU;abs(ux.u);abs(d.u)]);maxI=max([maxI;abs(d.ic)]);
  % Derivative validation excludes activation and priority-axis corners.
  for z=linspace(lo+.02*(hi-lo),hi-.02*(hi-lo),15)
   % Square-root priority corners have unbounded one-sided derivatives.
   % Use a five-point stencil inside the branch and scaled error criteria.
   gr=deg2rad(z);h=min(1e-5,.002*deg2rad(min(z-lo,hi-z)));v=map(gr,p,n,strategy,2);vp=map(gr+h,p,n,strategy,2);vm=map(gr-h,p,n,strategy,2);vp2=map(gr+2*h,p,n,strategy,2);vm2=map(gr-2*h,p,n,strategy,2);
   validation.PeDerivative=max(validation.PeDerivative,abs((-vp2.Pe+8*vp.Pe-8*vm.Pe+vm2.Pe)/(12*h)-v.Pg)/(1+abs(v.Pg)));
   validation.SecondDerivative=max(validation.SecondDerivative,abs((-vp2.Pg+8*vp.Pg-8*vm.Pg+vm2.Pg)/(12*h)-v.Pgg)/(1+abs(v.Pgg)));
   nodes=n.nodeCoefficients*[v.vc;p.U];kcl=(v.vc-nodes(1))/n.Zt1+n.Yc*v.vc-v.ic;validation.KCL=max(validation.KCL,abs(kcl));
   J=[0 1;-v.Pg/p.mass -p.D/p.mass];assert((v.Pg>0&&max(real(eig(J)))<0)||(v.Pg<0&&det(J)<0));
  end
 end
 assert(validation.PeDerivative<1e-6&&validation.SecondDerivative<1e-4&&validation.KCL<1e-11&&maxI<=p.Imax+1e-12, ...
 'Strategy %d: derivative %.12g, second %.12g, KCL %.12g, I %.12g',strategy,validation.PeDerivative,validation.SecondDerivative,validation.KCL,maxI);
 % Exact activation-side derivatives checked using one-sided differences.
 for side=1:2
  h=1e-7;sgn=2*side-3;g=deg2rad(sw);v=map(g,p,n,strategy,side);v1=map(g+sgn*h,p,n,strategy,side);v2=map(g+sgn*2*h,p,n,strategy,side);
  assert(abs(sgn*(-3*v.Pe+4*v1.Pe-v2.Pe)/(2*h)-v.Pg)<1e-6);
 end
 % Do not join intervals across priority corners (ordinary derivative absent).
 stable=join_within_branches(stable,cuts);saddles=join_within_branches(saddles,cuts);voltage=join_within_branches(voltage,cuts);usable=join_within_branches(usable,cuts);exportUsable=join_within_branches(exportUsable,cuts);
 rows{strategy}={names(strategy),sw,left.Pe,left.Pg,right.Pg,interval_text(stable),maxI,minU,maxU,interval_text(usable),interval_text(exportUsable)};
 results{strategy}=struct('network',n,'strategy',names(strategy),'switch_deg',sw,'left',left,'right',right,'priorityCorners_deg',cuts,'slopeRoots_deg',slopeRoots, ...
  'stableActive_deg',stable,'saddleActive_deg',saddles,'voltageFeasibleActive_deg',voltage,'feasibleStableActive_deg',usable,'feasibleExportStable_deg',exportUsable,'validation',validation);
end
fprintf(' Strategy             sw(deg)    Pe_sw     Pe_g^-     Pe_g^+    Stable active (deg)*                         max|I|   |u| range      Feasible stable (deg)*                      Feasible stable/export (deg)*\n');
fprintf('                                  pu      pu/rad     pu/rad     *sw<gamma<=180; open smooth intervals         pu        pu\n');
for k=1:3,r=rows{k};fprintf(' %-20s %9.6f %9.6f %+10.6f %+10.6f %-43s %.4f  %.4f..%.4f  %-43s %s\n',r{:});end

function n=rebuild(x,p,n0)
n=n0;z1=n0.Zt1+1i*x;z2=n0.Zt2;zg=n0.Zgrid;zb=n0.Zbreaker;
y=[1/z1+1/z2+n0.Ym,-1/z2;-1/z2,1/z2+1/(zg+zb)];co=y\diag([1/z1,1/(zg+zb)]);
iv=co(2,1)/(zg+zb);iu=(co(2,2)-1)/(zg+zb);av=co(2,1)-zb*iv;au=co(2,2)-zb*iu;icv=(1-co(1,1))/z1+n0.Yc;icu=-co(1,2)/z1;
n.Zt1=z1;n.nodalY=y;n.nodeCoefficients=co;n.vgE=av;n.vgU=au;n.igE=iv;n.igU=iu;n.icE=icv;n.icU=icu;n.Zeq=1/icv;n.Vth=-icu/icv*p.U;
end
function d=map(g,p,n,strategy,mode)
theta=g+angle(n.Vth);rot=exp(1i*theta);v=p.E*rot;q=(v-n.Vth)/n.Zeq;
% r(gamma)=(E-|Vth| exp(-j gamma))/Zeq is the raw-reference proxy
% expressed in the VSM frame, not in the global network phasor frame.
r=(p.E-abs(n.Vth)*exp(-1i*g))/n.Zeq;rg=1i*abs(n.Vth)*exp(-1i*g)/n.Zeq;rgg=abs(n.Vth)*exp(-1i*g)/n.Zeq;
active=abs(r)>p.Imax;if mode==1,active=false(size(g));elseif mode==2,active=true(size(g));end
s=r;sg=rg;sgg=rgg;
if strategy==1
 a=active;if any(a(:)),ra=abs(r(a));qa=r(a);ga=rg(a);gga=rgg(a);ag=real(conj(qa).*ga)./ra;agg=(abs(ga).^2+real(conj(qa).*gga))./ra-ag.^2./ra;
 s(a)=p.Imax*qa./ra;sg(a)=p.Imax*(ga./ra-qa.*ag./ra.^2);sgg(a)=p.Imax*(gga./ra-2*ga.*ag./ra.^2-qa.*agg./ra.^2+2*qa.*ag.^2./ra.^3);end
else
 a=real(r);ag=real(rg);agg=real(rgg);b=imag(r);if strategy==3,a=imag(r);ag=imag(rg);agg=imag(rgg);b=real(r);end
 full=active & abs(a)>=p.Imax;partial=active & ~full;
 sa=a;sag=ag;sagg=agg;sb=b;sbg=zeros(size(g));sbgg=sbg;
 sa(full)=sign(a(full))*p.Imax;sag(full)=0;sagg(full)=0;sb(full)=0;
 f=sqrt(max(0,p.Imax^2-a(partial).^2));sig=sign(b(partial));
 sb(partial)=sig.*f;sbg(partial)=-sig.*a(partial).*ag(partial)./f;
 sbgg(partial)=-sig.*((ag(partial).^2+a(partial).*agg(partial))./f+a(partial).^2.*ag(partial).^2./f.^3);
 if strategy==2,s(active)=sa(active)+1i*sb(active);sg(active)=sag(active)+1i*sbg(active);sgg(active)=sagg(active)+1i*sbgg(active);
 else,s(active)=sb(active)+1i*sa(active);sg(active)=sbg(active)+1i*sag(active);sgg(active)=sbgg(active)+1i*sagg(active);end
end
ic=rot.*s;icg=rot.*(sg+1i*s);icgg=rot.*(sgg+2i*sg-s);vc=n.Vth+n.Zeq*ic;vcg=n.Zeq*icg;vcgg=n.Zeq*icgg;
vg=n.vgE*vc+n.vgU*p.U;ig=n.igE*vc+n.igU*p.U;vg1=n.vgE*vcg;ig1=n.igE*vcg;vg2=n.vgE*vcgg;ig2=n.igE*vcgg;
d=struct('raw',r,'q',q,'ic',ic,'vc',vc,'Vgrid',vg,'Igrid',ig,'active',active,'Pe',real(vg.*conj(ig)), ...
 'Pg',real(vg1.*conj(ig)+vg.*conj(ig1)),'Pgg',real(vg2.*conj(ig)+2*vg1.*conj(ig1)+vg.*conj(ig2)), ...
 'u',vc+(n.Zconverter+p.AD)*ic,'ug',vcg+(n.Zconverter+p.AD)*icg);
end
function v=priority_axis(g,p,n,strategy)
r=map(deg2rad(g),p,n,1,0).raw;v=real(r);if strategy==3,v=imag(r);end
end
function v=voltage_slope(g,p,n,strategy)
d=map(deg2rad(g),p,n,strategy,2);v=real(conj(d.u).*d.ug);
end
function r=roots_grid(g,y,f)
ix=find(y(1:end-1).*y(2:end)<0);r=g(abs(y)<1e-12);for k=ix',r(end+1,1)=fzero(f,g(k+[0 1]),optimset('TolX',1e-12));end;r=uniquetol(r,1e-9,'DataScale',1);
end
function out=join_within_branches(in,corners)
out=zeros(0,2);for j=1:size(in,1),if ~isempty(out)&&abs(out(end,2)-in(j,1))<1e-8&&all(abs(corners-in(j,1))>1e-8),out(end,2)=in(j,2);else,out(end+1,:)=in(j,:);end;end
end
function t=interval_text(in)
if isempty(in),t="none";return;end
s=strings(size(in,1),1);for j=1:size(in,1),s(j)=sprintf('(%.6f,%.6f)',in(j,:));end;t=strjoin(s,' U ');
end
