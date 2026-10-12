function p=hd_set_HD(p0,H,Kd)
%HD_SET_HD Model parameter struct for virtual inertia H and damping gain Kd (K fixed).
% Reduced two-state model (physical time):
%   gamma_dot = omega_b*(nu-1),  nu_dot = ((Pstar-Pe)/nu - D*(nu-1))/m,  m = 2H.
% From the VSM block (tdmsquared Subsys_GFormConv_OuterCtrl_VSM.slx):
%   2H*omega_dot = (P*-P)/omega + K*(omega_ref-omega) - K_d*(omega-omega_PLL),
% with omega_ref = 1 and omega_PLL approximated by nominal 1:  D = K + K_d.
% The sweep keeps the droop gain K = p0.K (20) fixed and varies K_d, so D = K + K_d.
% H enters only through p.mass; D only through p.D. The electrical closure
% (hd_coupled) uses neither, so the equilibrium, Pe(gamma), the limiter switch and
% the voltage-command limit are independent of H and K_d.
assert(abs(p0.mass-2*p0.H)<1e-12,'Baseline mass is not 2H');
assert(abs(p0.D-(p0.K+p0.Kd))<1e-12,'Baseline D is not K+Kd');
assert(Kd>=0,'Negative K_d requested');
p=p0;p.H=H;p.mass=2*H;p.Kd=Kd;p.D=p0.K+Kd;
end
