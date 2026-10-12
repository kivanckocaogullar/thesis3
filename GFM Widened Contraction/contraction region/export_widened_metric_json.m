function export_widened_metric_json(src,out)
% Export saved widened metrics (no refit) to the JSON format read by the interval checker.
% src: a search_runs/widened_nu_* folder (inactive_metric.mat, active_metric.mat); out: e.g. metrics.
% The JSON in metrics/ was exported from the certified run; re-exporting rewrites metric.json.
% Also records sampled (non-rigorous) ranges of E, |q| and |u| to choose the root bracket.
root=fileparts(mfilename('fullpath'));a=load(fullfile(root,'inputs','model_inputs.mat'),'p','n');p=a.p;n=a.n;
if ~exist(out,'dir'),mkdir(out);end
names={'Zeq','Vth','vgE','vgU','igE','igU','Zconverter'};net=struct();for j=1:numel(names),zz=n.(names{j});net.(names{j})=[real(zz),imag(zz)];end
labels={'inactive','active'};mode=[1 0];
for j=1:2
 m=load(fullfile(src,[labels{j} '_metric.mat']));coef=cell(6,1);for k=1:6,coef{k}=m.C(:,:,k);end
 export=struct('coefficients',{coef},'lambda',m.lambda,'bounds',m.bounds,'center',m.center,'scale',m.scale,'epsilon',m.epsilon,'upper',m.upper,'network',net,'U',p.U,'AD',p.AD,'omega',p.omega,'mass',p.mass,'D',p.D,'Imax',p.Imax,'umax',p.Vmax);
 d=fullfile(out,labels{j});if ~exist(d,'dir'),mkdir(d);end
 fid=fopen(fullfile(d,'metric.json'),'w');fwrite(fid,jsonencode(export),'char');fclose(fid);
 g=linspace(m.bounds(1,1),m.bounds(2,1),2001);dd=coupled(g,1,p,n,mode(j));
 fprintf('%s gamma [%.12f %.12f] deg nu [%.6f %.6f] | E [%.6f %.6f] |q| [%.6f %.6f] |u| max %.6f | active all=%d any=%d\n',labels{j},rad2deg(m.bounds(:,1)),m.bounds(:,2),min(dd.E),max(dd.E),min(dd.rawI),max(dd.rawI),max(abs(dd.u)),all(dd.active),any(dd.active));
 % Float point samples for cross-checking the interval code (M, Mdot, Q eigenvalues).
 pts=[m.bounds(1,:);m.center;m.bounds(2,:);[m.bounds(1,1) m.bounds(2,2)];[m.bounds(2,1) m.bounds(1,2)]];chk=zeros(size(pts,1),5);
 for k=1:size(pts,1)
  x=pts(k,:);e=coupled(x(1),1,p,n,mode(j));f=[p.omega*(x(2)-1),((1-e.P)/x(2)-p.D*(x(2)-1))/p.mass];J=[0 p.omega;-e.totalPg/(p.mass*x(2)),(-(1-e.P)/x(2)^2-p.D)/p.mass];
  z=(x-m.center)./m.scale;dz=f./m.scale;b=[1 z z(1)^2 z(1)*z(2) z(2)^2];bd=[0 dz 2*z(1)*dz(1) z(2)*dz(1)+z(1)*dz(2) 2*z(2)*dz(2)];M=zeros(2);MD=M;for i=1:6,M=M+b(i)*m.C(:,:,i);MD=MD+bd(i)*m.C(:,:,i);end
  Q=J'*M+M*J+MD+2*m.lambda*M;chk(k,:)=[x,min(eig(M)),max(eig(M)),max(eig(Q))];
 end
 fid=fopen(fullfile(d,'matlab_point_checks.json'),'w');fwrite(fid,jsonencode(struct('gamma_rad',chk(:,1),'nu',chk(:,2),'minM',chk(:,3),'maxM',chk(:,4),'maxQ',chk(:,5))),'char');fclose(fid);
end
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
