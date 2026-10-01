%% Three-design physically consistent equilibrium-point screening only.
% Reference: baseline_gfm_contraction_SCR5_XR10.m; no existing files modified.
% Add an IDEAL series reactor between the physical capacitor node and the
% transformer primary series branch. Its impedance adds to Zt1 in the nodal
% reduction; Yc remains at the capacitor node, Ym at the transformer T node.
% External Zgrid, breaker, transformer ratio and converter RL stay unchanged.
% This is an analytical network design, NOT a saved Simulink modification.
% Keep physical capacitor-voltage magnitude E at its saved baseline value.
% Repeat the original parameter-dependent beta/c/Theta construction at each
% equilibrium and design; do not freeze old numeric matrices or optimize M.
% Pointwise M>0, Sym(F)<0 tests do NOT establish a uniform nonlinear or
% switched-system contraction certificate over the reported angle intervals.
root=fileparts(mfilename('fullpath'));
a=load(fullfile(root,'baseline_SCR5_XR10','SCR5_XR10_baseline_results.mat'),'results');
b=load(fullfile(root,'current_limiter_SCR5_XR10','current_limiter_SCR5_XR10_results.mat'),'results');
p0=a.results.parameters;n0=a.results.network;
assert(p0.H==3.4&&p0.K==20&&p0.Kd==5&&p0.D==25);
Imax=b.results.parameters.Imax;Vmax=b.results.parameters.Vmax;AD=b.results.parameters.AD;
assert(Imax==1.2&&Vmax==1.3&&AD==.2);
targets=[imag(n0.Zeq),.325,.350];designs=cell(3,1);
gd=(0:.0005:180)';comparison=table();
for k=1:3
    extraX=0;
    if k>1
        extraX=fzero(@(x)imag(reduce_network(x,p0,n0).Zeq)-targets(k),[0,.5],optimset('TolX',1e-13));
    end
    n=reduce_network(extraX,p0,n0);p=power_parameters(p0,n);
    assert(abs(imag(n.Zeq)-targets(k))<1e-11);
    assert(n.Zgrid==n0.Zgrid&&n.Zconverter==n0.Zconverter);
    if k==1
        assert(norm(n.nodeCoefficients-n0.nodeCoefficients,inf)<1e-12);
        assert(norm(p.pc-p0.pc,inf)<1e-12&&abs(n.Vth-n0.Vth)<1e-12);
    end
    current=@(g)abs((p.E*exp(1i*(deg2rad(g)+p.phi))-n.Vth)/n.Zeq);
    gammaSwitch=fzero(@(g)current(g)-Imax,[0 180],optimset('TolX',1e-12));
    e=electric(deg2rad(gd),p.U,p);roots=zero_roots(gd,e.Pe-1,@(g)electric(deg2rad(g),p.U,p).Pe-1);
    rootSlope=electric(deg2rad(roots),p.U,p).Pg;
    nominal=roots(find(rootSlope>0,1));assert(~isempty(nominal));
    ns=screen(nominal,p,n,Imax,Vmax,AD);
    assert(ns.currentMargin>0&&ns.voltageMargin>0,'Nominal Pstar=1 is not an unsaturated voltage-feasible equilibrium.');
    s=screen(gd,p,n,Imax,Vmax,AD);
    regions=intervals(gd,s.certificateMargin,@(g)screen(g,p,n,Imax,Vmax,AD).certificateMargin);
    usable=intervals(gd,s.usableMargin,@(g)screen(g,p,n,Imax,Vmax,AD).usableMargin);
    voltage=intervals(gd,s.voltageMargin,@(g)screen(g,p,n,Imax,Vmax,AD).voltageMargin);
    if k==1
        saved=a.results.regions_deg;saved=saved(saved(:,1)>0 & saved(:,2)<180,:);
        assert(isequal(size(regions),size(saved))&&max(abs(regions-saved),[],'all')<1e-7,'Baseline intervals were not reproduced.');
    end
    % Direct phasor, KCL and matrix checks against the compact expressions.
    verify_design(p,n,n0,extraX,Imax,Vmax,AD);
    R1="none";R2="none";if size(regions,1)>=1,R1=angle_intervals(regions(1,:));end
    if size(regions,1)>=2,R2=angle_intervals(regions(2,:));end
    assert(size(regions,1)<=2,'More than two regions: expand the console table.');
    extraL_mH=extraX*(p.Vlv^2/p.S)/(2*pi*p.f)*1e3;
    row=table(imag(n.Zeq),extraX,extraL_mH,string(sprintf('%.6f%+.6fj',real(n.Zeq),imag(n.Zeq))), ...
        string(sprintf('%.6f%+.6fj',real(n.Vth),imag(n.Vth))),nominal,gammaSwitch,R1,R2,angle_intervals(usable), ...
        'VariableNames',{'Xeq_pu','Added_X_LV_pu','Added_L_LV_mH','Zeq_pu','Vth_pu','Nominal_deg','Switch_deg','R1_deg','R2_deg','Usable_deg'});
    comparison=[comparison;row]; %#ok<AGROW>
    designs{k}=struct('parameters',p,'network',n,'addedX_LV_pu',extraX,'addedL_LV_mH',extraL_mH, ...
        'powerCoefficients',p.pc,'nominal_deg',nominal,'nominalCurrent_pu',current(nominal), ...
        'nominalVoltageCommand_pu',Vmax-ns.voltageMargin,'switch_deg',gammaSwitch, ...
        'regions_deg',regions,'voltageFeasible_deg',voltage,'usable_deg',usable);
end
% One compact console table only. No files, figures or simulations generated.
fprintf(' Xeq     Added X / L(mH), LV      Zeq (pu)              Vth (pu)              Nominal   Switch    R1 (deg)             R2 (deg)             Usable (deg)\n');
for k=1:height(comparison)
    t=comparison(k,:);
    fprintf(' %.6f  %.6f / %.6f   %-20s %-20s %8.4f  %8.4f  %-20s %-20s %s\n', ...
        t.Xeq_pu,t.Added_X_LV_pu,t.Added_L_LV_mH,t.Zeq_pu,t.Vth_pu,t.Nominal_deg,t.Switch_deg,t.R1_deg,t.R2_deg,t.Usable_deg);
end

function n=reduce_network(extraX,p,n0)
n=n0;z1=n0.Zt1+1i*extraX;z2=n0.Zt2;ym=n0.Ym;yc=n0.Yc;zg=n0.Zgrid;zb=n0.Zbreaker;
y=[1/z1+1/z2+ym,-1/z2;-1/z2,1/z2+1/(zg+zb)];co=y\diag([1/z1,1/(zg+zb)]);
iv=co(2,1)/(zg+zb);iu=(co(2,2)-1)/(zg+zb);av=co(2,1)-zb*iv;au=co(2,2)-zb*iu;
icv=(1-co(1,1))/z1+yc;icu=-co(1,2)/z1;
n.Zt1=z1;n.nodalY=y;n.nodeCoefficients=co;n.vgE=av;n.vgU=au;n.igE=iv;n.igU=iu;n.icE=icv;n.icU=icu;
n.Zeq=1/icv;n.Vth=-icu/icv*p.U;n.addedSeriesReactor=1i*extraX;
end
function p=power_parameters(p,n)
p.phi=angle(n.Vth);av=n.vgE;au=n.vgU;iv=n.igE;iu=n.igU;
A=av*conj(iv);B=au*conj(iu);C=av*conj(iu)+au*conj(iv);D=1i*(av*conj(iu)-au*conj(iv));
c=real([A B C D]);q=imag([A B C D]);
p.pc=[c(1:2),c(3)*cos(p.phi)+c(4)*sin(p.phi),-c(3)*sin(p.phi)+c(4)*cos(p.phi)];
p.qc=[q(1:2),q(3)*cos(p.phi)+q(4)*sin(p.phi),-q(3)*sin(p.phi)+q(4)*cos(p.phi)];
end
function d=electric(g,U,p)
c=p.pc;h=c(3)*cos(g)+c(4)*sin(g);hg=-c(3)*sin(g)+c(4)*cos(g);
d=struct('Pe',c(1)*p.E^2+c(2)*U.^2+p.E*U.*h,'Pg',p.E*U.*hg,'PU',2*c(2)*U+p.E*h,'PgU',p.E*hg);
end
function s=screen(gd,p,n,Imax,Vmax,AD)
g=deg2rad(gd);d=electric(g,p.U,p);beta=-d.PU/p.mass;c=-p.omega0_pu*d.PgU/p.mass;r=c./beta;a=-d.Pg/p.mass;b=p.D/p.mass;
off=(1+a+b*r-r.^2)/2;lm=-b/2+hypot(b/2-r,off);
m11=beta.^2+c.^2;m12=beta.*c;m22=beta.^2;mm=(m11+m22+hypot(m11-m22,2*m12))/2;mi=beta.^4./mm;
cert=min(mi-p.tol,-lm-p.tol);cert(abs(beta)<=1e-12)=-1;
vc=p.E*exp(1i*(g+p.phi));ic=(vc-n.Vth)/n.Zeq;u=vc+(n.Zconverter+AD)*ic;
cm=Imax-abs(ic);vm=Vmax-abs(u);
% Strict interiors are reported; current-switch and metric boundaries are
% not claimed as smooth certified points. Voltage feasibility is <=Vmax.
s=struct('certificateMargin',cert,'currentMargin',cm,'voltageMargin',vm,'usableMargin',min(min(cert,cm),vm));
end
function I=intervals(g,y,f)
ok=y>0;beg=find(diff([false;ok])==1);fin=find(diff([ok;false])==-1);I=zeros(numel(beg),2);
for j=1:numel(beg)
    lo=g(beg(j));hi=g(fin(j));
    if beg(j)>1,lo=fzero(f,g(beg(j)+[-1 0]),optimset('TolX',1e-11));end
    if fin(j)<numel(g),hi=fzero(f,g(fin(j)+[0 1]),optimset('TolX',1e-11));end
    I(j,:)=[lo hi];
end
end
function r=zero_roots(g,y,f)
ix=find(y(1:end-1).*y(2:end)<0);r=g(abs(y)<1e-12);for k=ix',r(end+1,1)=fzero(f,g(k+[0 1]),optimset('TolX',1e-12));end;r=uniquetol(r,1e-9,'DataScale',1);
end
function txt=angle_intervals(I)
if isempty(I),txt="none";return;end
parts=strings(size(I,1),1);for j=1:size(I,1),parts(j)=sprintf('(%.6f,%.6f)',I(j,1),I(j,2));end;txt=strjoin(parts,' U ');
end
function verify_design(p,n,n0,extraX,Imax,Vmax,AD)
for g=deg2rad([0 10 20 30 50 70 90 120])
    vc=p.E*exp(1i*(g+p.phi));nodes=n.nodeCoefficients*[vc;p.U];vg=n.vgE*vc+n.vgU*p.U;ig=n.igE*vc+n.igU*p.U;
    ic=(vc-nodes(1))/(n0.Zt1+1i*extraX)+n.Yc*vc;q=(vc-n.Vth)/n.Zeq;d=electric(g,p.U,p);
    assert(abs(ic-q)<1e-11&&abs(real(vg*conj(ig))-d.Pe)<1e-11);
    assert(norm(n.nodalY*nodes-[vc/(n0.Zt1+1i*extraX);p.U/(n.Zgrid+n.Zbreaker)],inf)<1e-10);
    beta=-d.PU/p.mass;c=-p.omega0_pu*d.PgU/p.mass;T=[-beta 0;-c -beta];J=[0 1;-d.Pg/p.mass -p.D/p.mass];
    if abs(beta)>1e-7
        F=T*J/T;direct=min(min(eig(T'*T))-p.tol,-max(eig((F+F')/2))-p.tol);
        fast=screen(rad2deg(g),p,n,Imax,Vmax,AD).certificateMargin;assert(abs(direct-fast)<1e-8*max(1,abs(direct)));
    end
    h=1e-6;dp=electric(g+h,p.U,p);dm=electric(g-h,p.U,p);assert(abs((dp.Pe-dm.Pe)/(2*h)-d.Pg)<1e-7);
    up=electric(g,p.U+h,p);um=electric(g,p.U-h,p);assert(abs((up.Pe-um.Pe)/(2*h)-d.PU)<1e-7&&abs((up.Pg-um.Pg)/(2*h)-d.PgU)<1e-7);
end
end
