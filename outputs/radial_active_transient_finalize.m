root=fileparts(mfilename('fullpath'));out=fullfile(root,'radial_active_transient_20261001_152918');a=load(fullfile(out,'results.mat'));p=a.p;n=a.n;pf=a.pf;nf=a.nf;r=a.primary;
diary(fullfile(out,'final_checks.txt'));
% Keep existing simulations; refine only state-space candidate presentation.
postSwitch=fzero(@(g)coupled(deg2rad(g),1,p,n,0).rawI-p.Imax,[23 26]);
boundsFault=a.boxes{1}.bounds;boundsPost=a.boxes{2}.bounds;boundsPost(1,1)=deg2rad(postSwitch+.05);
summary=cell(2,12);domains=cell(2,1);
for mode=1:2
 pp=pf;nn=nf;bb=boundsFault;if mode==2,pp=p;nn=n;bb=boundsPost;end
 [G,V]=ndgrid(linspace(bb(1,1),bb(2,1),501),linspace(bb(1,2),bb(2,2),101));z=coupled(G,1,pp,nn,0);
 A=-z.totalPg./(p.mass*V);B=(-(1-z.P)./V.^2-p.D)/p.mass;
 lam=real((B+sqrt(complex(B.^2+4*p.omega*A)))/2);
 assert(all(z.active,'all')&&max(abs(z.u),[],'all')<1.3&&max(abs(z.ic),[],'all')<1.2+1e-12);
 err=0;
 for g=linspace(bb(1,1),bb(2,1),11)
  for v=linspace(bb(1,2),bb(2,2),7)
   q=coupled(g,1,pp,nn,0);J=[0 p.omega;-q.totalPg/(p.mass*v),(-(1-q.P)/v^2-p.D)/p.mass];JJ=zeros(2);h=1e-6;x=[g;v];
   for col=1:2,ee=zeros(2,1);ee(col)=h;JJ(:,col)=(field(x+ee,pp,nn)-field(x-ee,pp,nn))/(2*h);end
   err=max(err,max(abs(J-JJ),[],'all'));
  end
 end
 assert(err<1e-6);
 roots=a.rootsByMode{1};if mode==1,roots=a.rootsByMode{2};end
 for j=1:numel(roots),rr=roots{j};assert(~(deg2rad(rr.gamma)>=bb(1,1)&&deg2rad(rr.gamma)<=bb(2,1)&&bb(1,2)<=1&&bb(2,2)>=1));end
 ranges=@(y)[min(y,[],'all') max(y,[],'all')];
 domains{mode}=struct('bounds',bb,'rawMargin',min(z.rawI-1.2,[],'all'),'uMargin',1.3-max(abs(z.u),[],'all'),'J21',ranges(A),'J22',ranges(B),'Pg',ranges(z.totalPg),'Pe',ranges(z.P),'maxRealEig',ranges(lam),'JacobianError',err);
 summary(mode,:)={mode,rad2deg(bb(:,1))',bb(:,2)',ranges(z.P),ranges(z.totalPg),ranges(z.rawI),ranges(abs(z.u)),ranges(lam),min(z.rawI-1.2,[],'all'),1.3-max(abs(z.u),[],'all'),ranges(A),ranges(B)};
end
D=cell2table(summary,'VariableNames',{'Mode','gamma_deg','nu','Pe','TotalPg','RawI','uRequired','MaxRealEigenvalue','CurrentBranchMargin','VoltageMargin','J21','J22'});disp(D);writetable(D,fullfile(out,'candidate_domains.csv'));
fprintf('Primary clearing state: gamma %.12f deg nu %.12f\n',rad2deg(r.clearState(1)),r.clearState(2));
fprintf('Current on initially: %d; post-clearing flag %d; smooth release events %s\n',r.dFault.active(1),r.dPost.active(1),mat2str(r.currentEventsPost',12));
fprintf('Primary ranges gamma=[%.12f %.12f], nu=[%.12f %.12f], Pe=[%.12f %.12f], rawI=[%.12f %.12f], u=[%.12f %.12f]\n',min(rad2deg([r.xFault(:,1);r.xPost(:,1)])),max(rad2deg([r.xFault(:,1);r.xPost(:,1)])),min([r.xFault(:,2);r.xPost(:,2)]),max([r.xFault(:,2);r.xPost(:,2)]),min([r.dFault.P;r.dPost.P]),max([r.dFault.P;r.dPost.P]),min([r.dFault.rawI;r.dPost.rawI]),max([r.dFault.rawI;r.dPost.rawI]),min(abs([r.dFault.u;r.dPost.u])),max(abs([r.dFault.u;r.dPost.u])));
fprintf('Primary Pe / Eref at sag application: %.12f / %.12f; clearing before %.12f / %.12f; after %.12f / %.12f\n',r.dFault.P(1),r.dFault.E(1),r.dFault.P(end),r.dFault.E(end),r.dPost.P(1),r.dPost.E(1));
fprintf('Final state gamma %.12f nu %.12f; final power residual %.12g\n',rad2deg(r.finalState(1)),r.finalState(2),r.dPost.P(end)-1);
for m=1:2,for j=1:numel(a.rootsByMode{m}),rr=a.rootsByMode{m}{j};fprintf('Mode %d fixed-P equilibrium %.12f deg active %d u %.12f eig %s\n',m,rr.gamma,rr.active,rr.u,mat2str(rr.eig',12));end,end
% Assemble actual active time intervals, including instantaneous sag/clearing switches.
Arows=cell(0,5);
for k=1:numel(a.runs)
 rr=a.runs{k};
 for mode=1:2
  tt=rr.tFault;xx=rr.xFault;dd=rr.dFault;ee=rr.currentEventsFault;pp=pf;nn=nf;
  if mode==2,tt=rr.tPost;xx=rr.xPost;dd=rr.dPost;ee=rr.currentEventsPost;pp=p;nn=n;end
  cuts=unique([tt(1);ee(:);tt(end)]);
  for j=1:numel(cuts)-1
   mid=mean(cuts(j:j+1));x=interp1(tt,xx,mid,'pchip');q=coupled(x(1),1,pp,nn,0);
   if q.active,Arows(end+1,:)={k,mode,cuts(j),cuts(j+1),cuts(j+1)-cuts(j)};end
  end
 end
end
AT=cell2table(Arows,'VariableNames',{'Run','Mode','Start_s','End_s','Duration_s'});writetable(AT,fullfile(out,'active_time_intervals.csv'));
save(fullfile(out,'candidate_domains.mat'),'domains','D','postSwitch','AT');diary off;
function f=field(x,p,n)
 z=coupled(x(1),1,p,n,0);f=[p.omega*(x(2)-1);((1-z.P)/x(2)-p.D*(x(2)-1))/p.mass];
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




