root=fileparts(mfilename('fullpath'));out=fullfile(root,'radial_switch_connection_20261001');a=load(fullfile(root,'radial_active_transient_20261001_152918','results.mat'),'p','n','primary');p=a.p;n=a.n;
ge=0.36012957762291908;g=unique(sort([linspace(16.426345885756,33.858410380369,2501)*pi/180 ge]));d=coupled(g,1,p,n,0);W=cumtrapz(g,d.P-1);W=W-interp1(g,W,ge);ix=find(W<=.0063);vvlo=zeros(size(ix));vvhi=vvlo;for j=1:numel(ix),en=.0063-W(ix(j));vvlo(j)=fzero(@(v)6.8*100*pi*((v-1)^2/2+(v-1)^3/3)-en,[.99 1]);vvhi(j)=fzero(@(v)6.8*100*pi*((v-1)^2/2+(v-1)^3/3)-en,[1 1.01]);end
f=figure('Visible','off','Color','w','Position',[50 50 1450 800]);ax=axes(f);hold(ax,'on');
h1=patch([16.426345885756 24.323430419712 24.323430419712 16.426345885756],[.996527139105 .996527139105 1.002060731220 1.002060731220],[.35 .65 1],'FaceAlpha',.22,'EdgeColor',[.1 .4 .8]);
h2=patch([24.423430419712 33.858410380369 33.858410380369 24.423430419712],[.996174514620 .996174514620 1.005178107392 1.005178107392],[1 .65 .2],'FaceAlpha',.22,'EdgeColor',[.8 .4 .1]);
h3=plot([16.426345885756 24.373430419712 24.373430419712 16.426345885756 16.426345885756],[.996527139105 .996527139105 1.0025 1.0025 .996527139105],'--','Color',[.1 .3 .8],'LineWidth',1.5);
h4=patch([24.323430419712 24.423430419712 24.423430419712 24.323430419712],[.996527139105 .996527139105 1.0025 1.0025],[.7 .25 .75],'FaceAlpha',.55,'EdgeColor','none');
x=g(ix)*180/pi;h5=patch([x fliplr(x)],[vvlo fliplr(vvhi)],[.15 .7 .3],'FaceAlpha',.35,'EdgeColor',[0 .5 .2]);
yline(1,':','HandleVisibility','off');xline(24.373430419712,'-.','HandleVisibility','off');
h6=plot(a.primary.xPost(:,1)*180/pi,a.primary.xPost(:,2),':','Color',[.35 .35 .35],'LineWidth',1.3);
h7=plot(a.primary.clearState(1)*180/pi,a.primary.clearState(2),'rx','LineWidth',2,'MarkerSize',10);
h8=patch([24.45 24.50005 24.50005 24.45],[.99995 .99995 1.00005 1.00005],[0 .45 .15],'FaceAlpha',1,'EdgeColor','k');
e=jsondecode(fileread(fullfile(root,'radial_inactive_metric_20261001_170815','invariant_subset.json')));P=e.P;rho=str2double(e.rho);th=linspace(0,2*pi,100);ell=chol(P)\(sqrt(rho)*[cos(th);sin(th)]);h9=plot((ge+ell(1,:))*180/pi,1+ell(2,:),'k-','LineWidth',2);plot(ge*180/pi,1,'k.','MarkerSize',12,'HandleVisibility','off');
xlim([16 34.2]);ylim([.9959 1.0056]);grid on;box on;xlabel('gamma (deg)');ylabel('nu (pu)');title({'Nominal-voltage post-fault recovery: fixed Pstar = 1','Green energy sublevel is invariant and convergent; rectangles alone are not invariant'});
legend([h1 h2 h3 h4 h5 h6 h7 h8 h9],{'Original inactive contraction rectangle','Original active contraction rectangle','Validated inactive extension (same metric)','Validated one-sided switching layer','Proven recovery set: H <= 0.0063','Saved trajectory: numerical guide only','Primary clearing point: outside proven set','Certified active clearing subset','Terminal invariant ellipsoid'},'Location','northeastoutside');
text(17,.99599,'Energy contour drawn numerically for visualization; proof uses interval bounds.','FontSize',10);exportgraphics(f,fullfile(out,'recovery_state_space.png'),'Resolution',180);savefig(f,fullfile(out,'recovery_state_space.fig'));close(f);
function [M,Mdot]=metric(C,z,dz)
 b=[1 z z(1)^2 z(1)*z(2) z(2)^2];bd=[0 dz 2*z(1)*dz(1) z(2)*dz(1)+z(1)*dz(2) 2*z(2)*dz(2)];M=zeros(2);Mdot=M;for j=1:6,M=M+b(j)*C(:,:,j);Mdot=Mdot+bd(j)*C(:,:,j);end
end
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







