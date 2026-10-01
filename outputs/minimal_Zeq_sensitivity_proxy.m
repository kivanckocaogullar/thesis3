%% Minimal equivalent-impedance proxy: limiter activation angle only.
% Reference: baseline_gfm_contraction_SCR5_XR10.m (left unchanged).
% External SCR=5, X/R=10; H=3.4, K=20, Kd=5 remain fixed.
% These dynamics parameters do not enter the algebraic switching condition.
% Varying Zeq here is NOT a physically realized network change. No original
% nodal coefficients, power-flow calculation, or contraction test is used.
% Old R1 lower boundary is a comparison threshold, NOT a recalculated region.
E = 0.997477277664;
Vth = complex(0.990384568064, -0.001741046188);
Req = 0.0298789620445;
Imax = 1.2;
Xeq = [0.2335966; 0.250; 0.275; 0.300; 0.325];
old_R1_lower_deg = 22.3523784;
Zeq = complex(Req, Xeq);

% Alignment to angle(Vth) gives:
% |q|^2 = (E^2 + |Vth|^2 - 2*E*|Vth|*cos(gamma))/|Zeq|^2.
cos_gamma = (E^2 + abs(Vth)^2 - (Imax*abs(Zeq)).^2) ...
    ./ (2*E*abs(Vth));
assert(all(cos_gamma > -1 & cos_gamma < 1), ...
    'A positive switching root in (0,180) deg is required for each candidate.');
gamma_activation_deg = acosd(cos_gamma);
q_at_switch = (E*exp(1i*(deg2rad(gamma_activation_deg)+angle(Vth))) ...
    - Vth)./Zeq;
assert(max(abs(abs(q_at_switch)-Imax)) < 1e-12, ...
    'Switching-current residual exceeds tolerance.');
exceeds_old_R1 = gamma_activation_deg > old_R1_lower_deg;

% Console table only: no plots, result files, or baseline modifications.
fprintf(' Xeq (pu)   Positive activation (deg)   Exceeds old R1 lower?\n');
fprintf(' --------   -------------------------   ---------------------\n');
for k = 1:numel(Xeq)
    answer = 'No';
    if exceeds_old_R1(k), answer = 'Yes'; end
    fprintf(' %.7f          %12.8f                  %s\n', ...
        Xeq(k), gamma_activation_deg(k), answer);
end
