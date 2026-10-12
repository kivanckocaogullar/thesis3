function hd_export_metric_json(R,p,n,file)
%HD_EXPORT_METRIC_JSON Write a rectangle/metric in the format read by widened_interval_checker.py.
% Same fields as export_widened_metric_json.m; mass and D carry the swept H and Kd.
names={'Zeq','Vth','vgE','vgU','igE','igU','Zconverter'};net=struct();for j=1:numel(names),zz=n.(names{j});net.(names{j})=[real(zz),imag(zz)];end
coef=cell(6,1);for k=1:6,coef{k}=R.C(:,:,k);end
ex=struct('coefficients',{coef},'lambda',R.lambda,'bounds',R.bounds,'center',R.center,'scale',R.scale,'epsilon',R.epsilon,'upper',R.upper,'network',net,'U',p.U,'AD',p.AD,'omega',p.omega,'mass',p.mass,'D',p.D,'Imax',p.Imax,'umax',p.Vmax, ...
 'H',p.H,'K',p.K,'Kd',p.Kd,'status_note','Numerically grid-checked - not interval-certified');
d=fileparts(file);if ~isfolder(d),mkdir(d);end
fid=fopen(file,'w');fwrite(fid,jsonencode(ex),'char');fclose(fid);
end
