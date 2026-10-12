function hd_quick_test(root)
addpath(fullfile(root,'scripts'));cfg=hd_default_config(root);a=load(cfg.modelInputs,'p','n');p0=a.p;n=a.n;
dom=hd_domain_limits(p0,n,cfg);disp(dom);fprintf('floor %.3f deg (%s), ceiling %.3f deg (%s), eq %.6f, switch %.6f, gamma_u %.4f\n',rad2deg(dom.inactiveFloor),dom.inactiveFloorReason,rad2deg(dom.activeCeiling),dom.activeCeilingReason,rad2deg(dom.gamma_eq),rad2deg(dom.gamma_switch),rad2deg(dom.gamma_u));
p=hd_set_HD(p0,3.4,5);fprintf('baseline mass %.4f D %.4f equal to stored: %d\n',p.mass,p.D,isequal(p.mass,p0.mass)&&isequal(p.D,p0.D));
for b=1:2
 R=check_contraction_HD(cfg.baseBounds{b},b,p,n,cfg);fprintf('branch %d: %s | %s | tmin %.3g rounds %d minM %.3g maxQ %.3g cur %.4f volt %.4f | %.1f s\n',b,R.status,R.detail,R.tmin,R.rounds,R.minM,R.maxQ,R.current_margin,R.voltage_margin,R.runtime_s);
end
R=check_contraction_HD([cfg.baseBounds{1}(1,1) .98;cfg.baseBounds{1}(2,1) 1.02],1,hd_set_HD(p0,8,0),n,cfg);fprintf('H=8,Kd=0 inactive: %s | %s | %.1f s\n',R.status,R.detail,R.runtime_s);
end
