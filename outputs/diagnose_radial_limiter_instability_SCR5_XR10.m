%% Diagnose existing radial benchmark only; no redesign or optimization.
clearvars;clc;close all;format long g
root=fileparts(mfilename('fullpath'));out=fullfile(root,'diagnose_radial_limiter_instability_SCR5_XR10');if ~isfolder(out),mkdir(out);end
f=fopen(fullfile(out,'diagnostic_log.txt'),'w');fclose(f);diary(fullfile(out,'diagnostic_log.txt'));guard=onCleanup(@()diary('off'));
a=load(fullfile(root,'current_limiter_SCR5_XR10','current_limiter_SCR5_XR10_results.mat'),'results');ref=a.results;p=ref.parameters;n=ref.network;
assert(p.Imax==1.2&&p.Vmax==1.3&&p.H==3.4&&p.K==20&&p.Kd==5&&p.D==25&&p.AD==.2);
sw=fzero(@(g)electrical(deg2rad(g),p.U,p,n,0).Ireq-p.Imax,[15 18],optimset('TolX',1e-12));
gd=unique([linspace(0,90,90001).';sw]);[scan,phasors]=diagnostics(gd,p,n,0);[baseline,baselinePhasors]=diagnostics(gd,p,n,1);
% Ordinary derivatives at the exact nonsmooth point are not defined.
at=abs(gd-sw)<1e-9;cols={'Pe_gamma','Qe_gamma','C_V','C_I','C_Vmag','C_Imag','C_phase','max_real_eigJ','min_real_eigJ','lambda_max_SymF'};for k=1:numel(cols),scan.(cols{k})(at)=NaN;end
derivativeColumns=scan.Properties.VariableNames(contains(scan.Properties.VariableNames,'gamma') & ~strcmp(scan.Properties.VariableNames,'gamma_deg'));for k=1:numel(derivativeColumns),scan.(derivativeColumns{k})(at)=NaN;end;scan.detJ(at)=NaN;scan.smoothBranch=~at;
epsDeg=1e-5;sampleG=sw+[-2;-.1;-epsDeg;epsDeg;.1;1;2;5];[points,pointPhasors]=diagnostics(sampleG,p,n,0);
[left,leftPhasors]=diagnostics(sw,p,n,1);[right,rightPhasors]=diagnostics(sw,p,n,2);sides=[left;right];sides.side=["inactive limit";"active limit"];
roots=[];for mode=1:2,if mode==1,g=linspace(0,sw,10001).';else,g=linspace(sw,180,50001).';end;d=electrical(deg2rad(g),p.U,p,n,mode);rr=crossings(g,d.Pg,@(v)electrical(deg2rad(v),p.U,p,n,mode).Pg);roots=[roots;rr];end %#ok<AGROW>
[widePoints,widePhasors]=diagnostics([sw;20;40;60;90],p,n,2);qZero=fzero(@(g)qvalue(g,p,n),[sw 40]);
validation=verify_diagnostics(p,n,sw,scan,phasors,leftPhasors,rightPhasors);
writetable(scan,fullfile(out,'radial_limiter_instability_diagnostics.csv'));writetable(points,fullfile(out,'radial_limiter_switch_diagnostics.csv'));writetable(sides,fullfile(out,'exact_one_sided_switch_diagnostics.csv'));
writetable(widePoints,fullfile(out,'active_branch_mechanism_context.csv'));
figures(scan,baseline,sw,out);
results=struct('parameters',p,'network',n,'switch_deg',sw,'smoothSlopeRoots_deg',roots,'scan',scan,'phasors',phasors,'baseline',baseline,'baselinePhasors',baselinePhasors,'points',points,'pointPhasors',pointPhasors,'oneSided',sides,'oneSidedPhasors',{{leftPhasors,rightPhasors}},'validation',validation,'referenceIntervals',ref.intervals,'widePoints',widePoints,'widePhasors',widePhasors,'activeQZero_deg',qZero);
save(fullfile(out,'radial_limiter_diagnostics.mat'),'results','-v7.3');write_summary(results,out);diary off

function [t,d]=diagnostics(gd,p,n,mode)
g=deg2rad(gd(:));d=electrical(g,p.U,p,n,mode);d.Vtarget=p.E*exp(1i*(g+p.phi));d.q_gamma=1i*d.Vtarget/p.Z;d.Iconv_gamma=d.Ic_g;d.Vc_gamma=p.Z*d.Ic_g;d.Vgrid_gamma=n.vgE*d.Vc_gamma;d.Igrid_gamma=n.igE*d.Vc_gamma;
d.S=d.Vg.*conj(d.Ig);d.Qe=imag(d.S);d.C_V=real(d.Vgrid_gamma.*conj(d.Ig));d.C_I=real(d.Vg.*conj(d.Igrid_gamma));d.Qe_gamma=imag(d.Vgrid_gamma.*conj(d.Ig)+d.Vg.*conj(d.Igrid_gamma));
t=table(gd(:),d.active,'VariableNames',{'gamma_deg','limiter_active'});
zn={'Vtarget','q','Iconv','Vc','Vg','Ig'};dispnames={'Vtarget','q','Iconv','Vc','Vgrid','Igrid'};dz={1i*d.Vtarget,d.q_gamma,d.Iconv_gamma,d.Vc_gamma,d.Vgrid_gamma,d.Igrid_gamma};
for k=1:numel(zn)
 v=d.(zn{k});dv=dz{k};a=dispnames{k};ang=unwrap(angle(v));mag=abs(v);md=real(conj(v).*dv)./mag;ad=imag(conj(v).*dv)./mag.^2;
 t.([a '_real'])=real(v);t.([a '_imag'])=imag(v);t.([a '_abs'])=mag;t.([a '_angle_deg'])=rad2deg(ang);t.([a '_gamma_real'])=real(dv);t.([a '_gamma_imag'])=imag(dv);t.([a '_mag_gamma'])=md;t.([a '_angle_gamma'])=ad;
end
delta=deg2rad(t.Vgrid_angle_deg-t.Igrid_angle_deg);dd=t.Vgrid_angle_gamma-t.Igrid_angle_gamma;
t.delta_VI_deg=rad2deg(delta);t.delta_VI_gamma=dd;t.cos_delta_VI=cos(delta);t.Pe=d.Pe;t.Qe=d.Qe;t.Pe_gamma=d.Pg;t.Qe_gamma=d.Qe_gamma;t.C_V=d.C_V;t.C_I=d.C_I;
t.C_Vmag=t.Vgrid_mag_gamma.*t.Igrid_abs.*cos(delta);t.C_Imag=t.Vgrid_abs.*t.Igrid_mag_gamma.*cos(delta);t.C_phase=-t.Vgrid_abs.*t.Igrid_abs.*sin(delta).*dd;
t.Pe_polar=t.Vgrid_abs.*t.Igrid_abs.*cos(delta);t.u_required=abs(d.u);t.voltage_valid=abs(d.u)<=p.Vmax;
b=p.D/p.mass;disc=complex(b^2-4*d.Pg/p.mass);l1=(-b+sqrt(disc))/2;l2=(-b-sqrt(disc))/2;t.max_real_eigJ=max(real(l1),real(l2));t.min_real_eigJ=min(real(l1),real(l2));t.traceJ=-b+zeros(size(g));t.detJ=d.Pg/p.mass;t.Pstar=d.Pe;t.equilibriumResidual=zeros(size(g));
beta=-d.PU/p.mass;c=-d.PgU/p.mass;r=c./beta;off=(1-d.Pg/p.mass+b*r-r.^2)/2;lm=-b/2+hypot(b/2-r,off);lm(abs(beta)<1e-12)=NaN;t.lambda_max_SymF=lm;
end
function v=verify_diagnostics(p,n,sw,t,d,l,r)
v=struct('activeCurrentMagnitude',max(abs(abs(d.Iconv(d.active))-p.Imax)),'radialAngleDifference_rad',max(abs(angle(d.Iconv.*conj(d.q)))),'complexPower',max(abs(t.Pe-real(d.Vg.*conj(d.Ig)))),'polarPower',max(abs(t.Pe-t.Pe_polar)),'complexDecomposition',max(abs(d.Pg-t.C_V-t.C_I)),'polarDecomposition',max(abs(d.Pg-t.C_Vmag-t.C_Imag-t.C_phase)),'dqQSign',max(abs(t.Qe-(imag(d.Vg).*real(d.Ig)-real(d.Vg).*imag(d.Ig)))),'activeCurrentMagnitudeDerivative',max(abs(t.Iconv_mag_gamma(d.active))),'complexDerivatives',0,'magnitudeDerivatives',0,'phaseDerivatives',0,'PeDerivative',0,'QeDerivative',0,'Jacobian',0,'equilibriumResidual',0,'switchPowerContinuity',abs(l.Pe-r.Pe),'switchVoltageContinuity',abs(l.Vc-r.Vc),'oneSidedPeDerivative',0);
h=1e-6;
for gd=[.5:1:89.5,sw-.01,sw+.01]
 g=deg2rad(gd);[a,z]=diagnostics(gd,p,n,0);mode=1+double(z.active);[ap,zp]=diagnostics(rad2deg(g+h),p,n,mode);[am,zm]=diagnostics(rad2deg(g-h),p,n,mode);
 names={'q','Iconv','Vc','Vg','Ig'};ds={'q_gamma','Iconv_gamma','Vc_gamma','Vgrid_gamma','Igrid_gamma'};tn={'q','Iconv','Vc','Vgrid','Igrid'};
 for k=1:5
  v.complexDerivatives=max(v.complexDerivatives,abs((zp.(names{k})-zm.(names{k}))/(2*h)-z.(ds{k})));
  v.magnitudeDerivatives=max(v.magnitudeDerivatives,abs((abs(zp.(names{k}))-abs(zm.(names{k})))/(2*h)-a.([tn{k} '_mag_gamma'])));
  v.phaseDerivatives=max(v.phaseDerivatives,abs(angle(zp.(names{k})*conj(zm.(names{k})))/(2*h)-a.([tn{k} '_angle_gamma'])));
 end
 v.PeDerivative=max(v.PeDerivative,abs((zp.Pe-zm.Pe)/(2*h)-z.Pg));v.QeDerivative=max(v.QeDerivative,abs((zp.Qe-zm.Qe)/(2*h)-z.Qe_gamma));
 x=[g;1];ps=z.Pe;J=[0 1;-z.Pg/p.mass -p.D/p.mass];Jn=zeros(2);for k=1:2,eh=zeros(2,1);eh(k)=h;Jn(:,k)=(rhs_diag(x+eh,ps,p,n,mode)-rhs_diag(x-eh,ps,p,n,mode))/(2*h);end;v.Jacobian=max(v.Jacobian,norm(J-Jn,inf));v.equilibriumResidual=max(v.equilibriumResidual,norm(rhs_diag(x,ps,p,n,mode),inf));
 assert((z.Pg>0 && max(real(eig(J)))<0)||(z.Pg<0 && det(J)<0));
end
for mode=1:2
 g=deg2rad(sw);sgn=2*mode-3;d0=electrical(g,p.U,p,n,mode);d1=electrical(g+sgn*h,p.U,p,n,mode);d2=electrical(g+sgn*2*h,p.U,p,n,mode);fd=sgn*(-3*d0.Pe+4*d1.Pe-d2.Pe)/(2*h);v.oneSidedPeDerivative=max(v.oneSidedPeDerivative,abs(fd-d0.Pg));
end
assert(v.complexDerivatives<1e-6&&v.PeDerivative<1e-6&&v.QeDerivative<1e-6&&v.magnitudeDerivatives<1e-6&&v.phaseDerivatives<1e-4&&v.Jacobian<1e-6);
assert(v.activeCurrentMagnitude<1e-12&&v.radialAngleDifference_rad<1e-12&&v.polarPower<1e-12&&v.polarDecomposition<1e-11&&v.complexDecomposition<1e-11&&v.dqQSign<1e-12&&v.oneSidedPeDerivative<1e-6);
fprintf('NUMERICAL DIAGNOSTIC VALIDATION PASSED\n');disp(v);
end
function dx=rhs_diag(x,ps,p,n,mode)
d=electrical(x(1),p.U,p,n,mode);dx=[x(2)-1;((ps-d.Pe)/x(2)-p.D*(x(2)-1))/p.mass];
end
function figures(t,b,sw,out)
spec={ {'Pe'}, {'Pe_gamma'}, {'q_angle_deg','Iconv_angle_deg'}, {'Vc_abs','Vc_angle_deg'}, {'Pe','Qe'}, {'Pe_gamma','C_V','C_I'}, {'Pe_gamma','C_Vmag','C_Imag','C_phase'}, {'Iconv_abs','Iconv_angle_deg'}, {'delta_VI_deg','cos_delta_VI'}, {'max_real_eigJ'}, {'lambda_max_SymF'} };
names={'01_power','02_synchronizing_slope','03_current_direction','04_capacitor_voltage','05_PQ_transfer','06_complex_power_slope','07_polar_power_slope','08_current_clamping','09_voltage_current_phase','10_local_eigenvalue','11_metric_context'};
labels={'Grid Pe (pu)','dPe/dgamma (pu/rad)','Global current angle (deg)','Capacitor voltage','Grid P and Q (pu)','Power slope contributions (pu/rad)','Power slope contributions (pu/rad)','Converter current','Grid V-I phase relation','max Re eig(J)','max eig(Sym(F)) (not local stability)'};
for k=1:11
 sep=ismember(k,[4 8 9]);nr=1+sep;fig=figure('Color','w','Position',[40 40 1500 450*nr]);tl=tiledlayout(fig,nr,2,'TileSpacing','compact');fields=spec{k};
 for view=1:2
  limits=[0 90];if view==2,limits=[10 40];end
  for row=1:nr
   ax=nexttile(tl,(row-1)*2+view);hold(ax,'on');ff=fields;if sep,ff=fields(row);end
   if k==1,plot(ax,b.gamma_deg,b.Pe,'--','Color',[.5 .5 .5],'DisplayName','Unlimited baseline','LineWidth',1.4);end
   for j=1:numel(ff),y=t.(ff{j});if contains(ff{j},'gamma')||ismember(k,[6 7 10 11]),y(abs(t.gamma_deg-sw)<1e-9)=NaN;end;plot(ax,t.gamma_deg,y,'LineWidth',1.5,'DisplayName',strrep(ff{j},'_',' '));end
   xline(ax,sw,'k--','Current switch','HandleVisibility','off');if ismember(k,[2 6 7 10 11]),yline(ax,0,'k:','HandleVisibility','off');end
   if k==8&&row==1,yline(ax,1.2,'r:','Imax','HandleVisibility','off');end
   if k==11,ylim(ax,[-1 30]);end
   xlim(ax,limits);grid(ax,'on');xlabel(ax,'VSM voltage-target gamma (deg)');ylabel(ax,labels{k});if sep,ylabel(ax,strrep(ff{1},'_',' '));if k==4&&row==1,ylabel(ax,'|Vc| (pu)');elseif k==4,ylabel(ax,'Global angle Vc (deg)');elseif k==8&&row==1,ylabel(ax,'|Iconv| (pu)');elseif k==8,ylabel(ax,'Global angle Iconv (deg)');elseif row==1,ylabel(ax,'Relative delta VI (deg)');else,ylabel(ax,'cos(delta VI)');end;end;legend(ax,'Location','best','FontSize',9);
   if ismember(k,[6 7])&&view==2,xlim(ax,[sw 40]);idx=t.gamma_deg>sw+1e-8 & t.gamma_deg<=40;yy=t{idx,ff};lo=min(yy,[],'all');hi=max(yy,[],'all');pad=.1*(hi-lo);ylim(ax,[lo-pad hi+pad]);title(ax,'Limiter-active contribution close-up');end
   if k==2,title(ax,'Left of switch: locally stable; right: saddle (smooth branches)');end
  end
 end
 title(tl,{'SCR5/XR10 radial-current benchmark diagnostics | Imax=1.2, H=3.4, D=25',strrep(names{k},'_',' ')});axs=findall(fig,'Type','axes');for j=1:numel(axs),axs(j).Toolbar.Visible='off';end;drawnow;exportgraphics(fig,fullfile(out,[names{k},'.png']),'Resolution',170);savefig(fig,fullfile(out,[names{k},'.fig']));close(fig);
end
end
function write_summary(r,out)
s=r.oneSided;fprintf('\nSWITCH %.12f deg\n',r.switch_deg);disp(s(:,{'Pe','Qe','Pe_gamma','Qe_gamma','Iconv_mag_gamma','Iconv_angle_gamma','Vc_mag_gamma','Vc_angle_gamma','C_V','C_I','C_Vmag','C_Imag','C_phase','max_real_eigJ'}));fprintf('Smooth Pe_gamma zero roots over 0..180 deg: ');disp(r.smoothSlopeRoots_deg);
fid=fopen(fullfile(out,'radial_limiter_instability_summary.md'),'w','n','UTF-8');guard=onCleanup(@()fclose(fid));
fprintf(fid,'# Radial-current limiter: instability diagnosis\n\nSame saved SCR5/XR10 network and parameters; no redesign, tuning, new angle parameter, metric optimization or EMT. This is the ideal two-state radial benchmark, not the complete PI/freeze controller.\n\n');
fprintf(fid,'Current switch: **%.12f deg**. At this point Pe and Vc are continuous, but their derivatives need not be. No smooth Pe_gamma=0 root was found on either positive-angle branch over 0..180 deg. The slope changes SIGN BY A JUMP at the switch, not by continuously crossing zero. Therefore there is no positive-width stable active interval.\n\n',r.switch_deg);
fprintf(fid,'## Exact one-sided values\n\nAll magnitude/power derivatives are per RADIAN of target gamma. Phase derivatives are rad/rad (numerically also deg/deg). Global phasors use the referred physical grid bus as zero angle; gamma itself references Vth.\n\n| Quantity | Inactive side | Active side |\n|---|---:|---:|\n');fields={'Pe','Qe','Pe_gamma','Qe_gamma','Iconv_mag_gamma','Iconv_angle_gamma','Vc_mag_gamma','Vc_angle_gamma','Vgrid_mag_gamma','Igrid_mag_gamma','delta_VI_gamma','C_V','C_I','C_Vmag','C_Imag','C_phase','max_real_eigJ'};for k=1:numel(fields),fprintf(fid,'| %s | %.12g | %.12g |\n',fields{k},s.(fields{k})(1),s.(fields{k})(2));end
fprintf(fid,'\n## What removes synchronizing power?\n\nBefore limiting both the magnitude and direction of converter current can change. The prescribed capacitor voltage has fixed magnitude E and follows the target angle with unit angle derivative. At saturation the radial component of the current derivative is removed: dIconv/dgamma=(Imax/abs(q)) [q_gamma-q real(conj(q)q_gamma)/abs(q)^2]. The surviving derivative is tangent to the current circle. Current direction is exactly that of q; there is no independent limiter-angle control.\n\nThe capacitor voltage must then move according to Vc=Vth+Zeq Iconv rather than remain equal to Vtarget. Its magnitude decreases and its angle barely follows the target near switching. Grid current is not the converter current because shunt/magnetizing currents and transformer referral remain in the network.\n\n');
fprintf(fid,'At the active-side limit the COMPLEX decomposition is C_V=%.9g and C_I=%.9g, totaling %.9g pu/rad. The POLAR decomposition is voltage magnitude %.9g, grid-current magnitude %.9g, and V-I phase contribution %.9g pu/rad. The signed terms, rather than current magnitude clamping alone, explain the negative synchronizing slope.\n\n',s.C_V(2),s.C_I(2),s.Pe_gamma(2),s.C_Vmag(2),s.C_Imag(2),s.C_phase(2));
fprintf(fid,'Reactive-power slope changes from %.9g to %.9g pu/rad while active-power slope changes from %.9g to %.9g. Positive Q is verified against Q=v_gq*i_gd-v_gd*i_gq. The stored P/Q curves distinguish reactive increase from changes in apparent-power magnitude; no conservation of P+Q or fixed grid-current magnitude is assumed.\n\n',s.Qe_gamma(1),s.Qe_gamma(2),s.Pe_gamma(1),s.Pe_gamma(2));
fprintf(fid,'**Near the switch the dominant negative term is grid-voltage magnitude reduction**, not V-I phase rotation: the phase contribution is initially positive. The current-magnitude growth term that previously supported P is removed and the reconstructed Vc/grid voltage can sag. At larger angles the phase term becomes negative and its importance changes, as quantified below.\n\n| gamma deg | Pe | Qe | Voltage-magnitude contribution | Grid-current-magnitude contribution | Phase contribution |\n|---|---:|---:|---:|---:|---:|\n');
for k=1:height(r.widePoints),w=r.widePoints(k,:);fprintf(fid,'| %.9f | %.9g | %.9g | %.9g | %.9g | %.9g |\n',w.gamma_deg,w.Pe,w.Qe,w.C_Vmag,w.C_Imag,w.C_phase);end
fprintf(fid,'\nOn the active branch Q first DECREASES from positive injection, crosses zero at %.12f deg, and subsequently becomes negative (reactive absorption under this convention). Thus the initial loss of P slope is NOT explained by an increase in positive reactive-power injection. The power-factor cosine initially improves; farther along the active branch growing phase separation contributes to declining P.\n\n',r.activeQZero_deg);
fprintf(fid,'## Stability versus certificate\n\nAt each point Pstar=Pe and nu=1. trace(J)=-25/6.8 and det(J)=Pe_gamma/6.8. Thus a negative slope gives a saddle regardless of the chosen contraction metric. The ordinary Jacobian is undefined at the exact switch; the one-sided stable and saddle classifications meet there. On the positive principal interval 0..180 deg: inactive locally stable 0<=gamma<%.12f; active locally stable EMPTY; active saddle %.12f<gamma<=180. Displayed plots cover 0..90 and 10..40 deg.\n\n',r.switch_deg,r.switch_deg);
fprintf(fid,'The unchanged metric already fails to certify the locally stable pre-switch branch (including nominal): its positive maximum symmetric eigenvalue must not be confused with local instability. At the switch genuine local instability occurs in addition. The negative limited power slope is a property of this electrical equilibrium family; the earlier pre-switch certificate failure is metric conservativeness.\n\n');
fprintf(fid,'## Validation\n\n');fn=fieldnames(r.validation);for k=1:numel(fn),fprintf(fid,'- %s: %.9g\n',fn{k},r.validation.(fn{k}));end
fprintf(fid,'\nCSV includes Cartesian phasors and derivatives, continuous unwrapped global angles, relative V-I angle, polar derivatives, both power-slope decompositions, grid P/Q, voltage capability, Jacobian spectrum and metric context. Exact one-sided diagnostics are separately exported; smooth derivative columns at the exact switch are marked NaN. MAT retains all figure source data.\n');
end
function d=electrical(g,U,p,n,mode)
% mode 0 automatic; 1 force inactive; 2 force active (one-sided derivatives).
v=p.E*exp(1i*(g+p.phi));q=(v-p.t*U)/p.Z;qg=1i*v/p.Z;qU=-p.t/p.Z+zeros(size(q));r=abs(q);
active=r>p.Imax;if mode==1,active=false(size(q));elseif mode==2,active=true(size(q));end
ic=q;icg=qg;icU=qU;icgU=zeros(size(q));
if any(active(:))
 a=active;qa=q(a);ra=r(a);ga=qg(a);ua=qU(a);rg=real(conj(qa).*ga)./ra;rU=real(conj(qa).*ua)./ra;rgU=real(conj(ua).*ga)./ra-rg.*rU./ra;
 ic(a)=p.Imax*qa./ra;icg(a)=p.Imax*(ga./ra-qa.*rg./ra.^2);icU(a)=p.Imax*(ua./ra-qa.*rU./ra.^2);
 icgU(a)=p.Imax*(-ga.*rU./ra.^2-ua.*rg./ra.^2-qa.*rgU./ra.^2+2*qa.*rg.*rU./ra.^3);
end
vc=p.t*U+p.Z*ic;vcg=p.Z*icg;vcU=p.t+p.Z*icU;vcgU=p.Z*icgU;
% On the inactive branch set exact target to avoid cancellation roundoff.
vc(~active)=v(~active);vcg(~active)=1i*v(~active);vcU(~active)=0;vcgU(~active)=0;
vg=n.vgE*vc+n.vgU*U;ig=n.igE*vc+n.igU*U;
vgg=n.vgE*vcg;igg=n.igE*vcg;vgU=n.vgE*vcU+n.vgU;igU=n.igE*vcU+n.igU;
vgm=n.vgE*vcgU;igm=n.igE*vcgU;
Pe=real(vg.*conj(ig));Pg=real(vgg.*conj(ig)+vg.*conj(igg));PU=real(vgU.*conj(ig)+vg.*conj(igU));PgU=real(vgm.*conj(ig)+vgg.*conj(igU)+vgU.*conj(igg)+vg.*conj(igm));
vb=vc+n.Zconverter*ic;u=vb+p.AD*ic;
Uvec=U+zeros(size(vc));
d=struct('Pe',Pe,'Pg',Pg,'PU',PU,'PgU',PgU,'q',q,'Ireq',r,'Iconv',ic,'Vc',vc,'Vg',vg,'Ig',ig,'u',u,'Vbridge',vb,'active',active,'Ic_g',icg,'Ic_U',icU,'nodes',n.nodeCoefficients*[vc(:).';Uvec(:).']);
end
function roots=crossings(g,y,fun)
idx=find(y(1:end-1).*y(2:end)<0);roots=g(y==0);for k=idx(:).',roots(end+1,1)=fzero(fun,g(k+[0 1]),optimset('TolX',1e-11));end;roots=unique(sort(roots));
end

function q=qvalue(g,p,n)
d=electrical(deg2rad(g),p.U,p,n,2);q=imag(d.Vg*conj(d.Ig));
end


