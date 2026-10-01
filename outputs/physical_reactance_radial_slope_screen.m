%% Physical downstream-reactor variants: no metric, plots, EMT, or file writes.
% Reactor is between capacitor node and transformer primary series branch.
% External SCR5/XR10 RL, E, H/K/Kd, converter RL and AD remain unchanged.
% Reconstruct the transformer/shunt/breaker network for every candidate.
% Positive-angle domain is the principal interval 0..180 degrees.
% Derivatives are per RADIAN. At the switch use separate one-sided limits.
root=fileparts(mfilename('fullpath'));
a=load(fullfile(root,'current_limiter_SCR5_XR10','current_limiter_SCR5_XR10_results.mat'),'results');
p=a.results.parameters;n0=a.results.network;
assert(p.Imax==1.2&&p.Vmax==1.3&&p.AD==.2&&p.H==3.4&&p.K==20&&p.Kd==5);
nb=rebuild(0,p,n0);assert(abs(nb.Zeq-n0.Zeq)<1e-12&&abs(nb.Vth-n0.Vth)<1e-12);
targets=[.325 .350 .400];designs=cell(3,1);rows=cell(3,1);
for k=1:3
 x=fzero(@(z)imag(rebuild(z,p,n0).Zeq)-targets(k),[0 .6],optimset('TolX',1e-13));
 n=rebuild(x,p,n0);assert(abs(imag(n.Zeq)-targets(k))<1e-11&&n.Zgrid==n0.Zgrid);
 sw=fzero(@(g)abs(electrical(deg2rad(g),p,n,0).q)-p.Imax,[0 180],optimset('TolX',1e-12));
 left=electrical(deg2rad(sw),p,n,1);right=electrical(deg2rad(sw),p,n,2);
 assert(abs(left.Pe-right.Pe)<1e-11&&abs(left.vc-right.vc)<1e-11);
 % Locate slope extrema first, then bracket every slope root on monotone
 % pieces. This avoids mistaking a derivative jump for a smooth zero.
 gd=linspace(sw,180,180001)';d=electrical(deg2rad(gd),p,n,2);
 extrema=find_roots(gd,d.Pgg,@(g)electrical(deg2rad(g),p,n,2).Pgg);
 cuts=unique([sw;extrema;180]);slopes=electrical(deg2rad(cuts),p,n,2).Pg;
 slopeRoots=find_roots(cuts,slopes,@(g)electrical(deg2rad(g),p,n,2).Pg);
 edges=unique([sw;slopeRoots;180]);positive=zeros(0,2);
 for j=1:numel(edges)-1
  if electrical(deg2rad(mean(edges(j:j+1))),p,n,2).Pg>0,positive(end+1,:)=edges(j:j+1)';end %#ok<AGROW>
 end
 maxActiveSlope=max(slopes);assert(all(abs(d.q)>=p.Imax-1e-12));
 % Solve the actual piecewise Pe=1 map, not a forced-current/unlimited-power mix.
 allg=(0:.001:180)';fullscan=electrical(deg2rad(allg),p,n,0);
 er=find_roots(allg,fullscan.Pe-1,@(g)electrical(deg2rad(g),p,n,0).Pe-1);
 es=electrical(deg2rad(er),p,n,0);jj=find(es.Pg>0,1);
 nominal=NaN;in=NaN;un=NaN;feasible=false;
 if ~isempty(jj)
  nominal=er(jj);en=electrical(deg2rad(nominal),p,n,0);in=abs(en.ic);un=abs(en.u);
  feasible=~en.active&&in<p.Imax&&un<=p.Vmax;assert(abs(en.Pe-1)<1e-10);
 end
 % Analytic derivatives validated separately on both smooth branches.
 for mode=1:2
  tests=[sw-.1,sw+.1,40,90,150];if mode==1,tests=[0,sw/2,sw-.1];end
  for z=tests
   g=deg2rad(z);h=1e-6;v=electrical(g,p,n,mode);vp=electrical(g+h,p,n,mode);vm=electrical(g-h,p,n,mode);
   assert(abs((vp.Pe-vm.Pe)/(2*h)-v.Pg)<1e-7&&abs((vp.Pg-vm.Pg)/(2*h)-v.Pgg)<1e-6);
   nodes=n.nodeCoefficients*[v.vc;p.U];ic=(v.vc-nodes(1))/n.Zt1+n.Yc*v.vc;assert(abs(ic-v.ic)<1e-11);
  end
 end
 for mode=1:2
  h=1e-6;sgn=2*mode-3;g=deg2rad(sw);v=electrical(g,p,n,mode);v1=electrical(g+sgn*h,p,n,mode);v2=electrical(g+sgn*2*h,p,n,mode);
  assert(abs(sgn*(-3*v.Pe+4*v1.Pe-v2.Pe)/(2*h)-v.Pg)<1e-7);
 end
 extraL=x*(p.Vlv^2/p.S)/(2*pi*p.f)*1e3;
 txt="none";if ~isempty(positive),parts=strings(size(positive,1),1);for j=1:size(positive,1),parts(j)=sprintf('(%.6f,%.6f)',positive(j,:));end;txt=strjoin(parts,' U ');end
 rows{k}={targets(k),extraL,sw,left.Pe,left.Pg,right.Pg,maxActiveSlope,txt,nominal,in,un,feasible};
 designs{k}=struct('network',n,'addedX_LV_pu',x,'addedL_LV_mH',extraL,'gamma_deg',allg,'Pe',fullscan.Pe,'Pe_gamma',fullscan.Pg, ...
  'switch_deg',sw,'inactiveLimit',left,'activeLimit',right,'activePositiveIntervals_deg',positive,'activeSlopeRoots_deg',slopeRoots,'nominal_deg',nominal,'nominalFeasible',feasible);
end
fprintf(' Xeq   Ladd_LV   gamma_sw  Pe_sw    Pe_g^-     Pe_g^+     max Pe_g active  Active Pe_g>0 (deg)*  gamma_nom   I_nom   u_nom  NomOK\n');
fprintf(' (pu)   (mH)      (deg)    (pu)   (pu/rad)   (pu/rad)     (pu/rad)         *sw<gamma<=180       (deg)      (pu)    (pu)   (0/1)\n');
for k=1:3,r=rows{k};fprintf(' %.3f  %.6f  %9.6f %.6f %+10.6f %+10.6f  %+12.6f    %-20s %10.6f  %.6f %.6f   %d\n',r{:});end

function n=rebuild(x,p,n0)
n=n0;z1=n0.Zt1+1i*x;z2=n0.Zt2;zg=n0.Zgrid;zb=n0.Zbreaker;
y=[1/z1+1/z2+n0.Ym,-1/z2;-1/z2,1/z2+1/(zg+zb)];co=y\diag([1/z1,1/(zg+zb)]);
iv=co(2,1)/(zg+zb);iu=(co(2,2)-1)/(zg+zb);av=co(2,1)-zb*iv;au=co(2,2)-zb*iu;
icv=(1-co(1,1))/z1+n0.Yc;icu=-co(1,2)/z1;
n.Zt1=z1;n.nodalY=y;n.nodeCoefficients=co;n.vgE=av;n.vgU=au;n.igE=iv;n.igU=iu;n.icE=icv;n.icU=icu;n.Zeq=1/icv;n.Vth=-icu/icv*p.U;
end
function d=electrical(g,p,n,mode)
v=p.E*exp(1i*(g+angle(n.Vth)));q=(v-n.Vth)/n.Zeq;r=abs(q);qg=1i*v/n.Zeq;qgg=-v/n.Zeq;
active=r>p.Imax;if mode==1,active=false(size(g));elseif mode==2,active=true(size(g));end
ic=q;ig1=qg;ig2=qgg;
if any(active(:))
 a=active;qa=q(a);ra=r(a);ga=qg(a);gga=qgg(a);rg=real(conj(qa).*ga)./ra;rgg=(abs(ga).^2+real(conj(qa).*gga))./ra-rg.^2./ra;
 ic(a)=p.Imax*qa./ra;ig1(a)=p.Imax*(ga./ra-qa.*rg./ra.^2);
 ig2(a)=p.Imax*(gga./ra-2*ga.*rg./ra.^2-qa.*rgg./ra.^2+2*qa.*rg.^2./ra.^3);
end
vc=n.Vth+n.Zeq*ic;vg1=n.Zeq*ig1;vg2=n.Zeq*ig2;vc(~active)=v(~active);
vg=n.vgE*vc+n.vgU*p.U;ig=n.igE*vc+n.igU*p.U;a=n.vgE*vg1;b=n.igE*vg1;aa=n.vgE*vg2;bb=n.igE*vg2;
d=struct('q',q,'ic',ic,'vc',vc,'active',active,'Pe',real(vg.*conj(ig)), ...
 'Pg',real(a.*conj(ig)+vg.*conj(b)),'Pgg',real(aa.*conj(ig)+2*a.*conj(b)+vg.*conj(bb)), ...
 'u',vc+(n.Zconverter+p.AD)*ic);
end
function r=find_roots(g,y,f)
ix=find(y(1:end-1).*y(2:end)<0);r=g(abs(y)<1e-12);for j=ix',r(end+1,1)=fzero(f,g(j+[0 1]),optimset('TolX',1e-12));end;r=uniquetol(r,1e-9,'DataScale',1);
end
