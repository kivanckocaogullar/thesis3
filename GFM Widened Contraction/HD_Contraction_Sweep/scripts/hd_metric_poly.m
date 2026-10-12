function [M,Mdot]=hd_metric_poly(C,z,dz)
% Verbatim copy of the local function 'metric' from
% GFM Contraction Analsys/limiter region/widen_nu_contraction.m (itself copied from
% radial_inactive_metric_rate_trial.m). Only function names were changed.
% H and D enter exclusively through p.mass (=2H) and p.D.
 b=[1 z z(1)^2 z(1)*z(2) z(2)^2];bd=[0 dz 2*z(1)*dz(1) z(2)*dz(1)+z(1)*dz(2) 2*z(2)*dz(2)];M=zeros(2);Mdot=M;for j=1:6,M=M+b(j)*C(:,:,j);Mdot=Mdot+bd(j)*C(:,:,j);end
end
