function s=hd_style()
%HD_STYLE Shared colours and typography for all H-D sweep figures.
s.inactive=[0.129 0.400 0.674];   % blue  (limiter inactive), same hue as earlier region plots
s.active=[0.843 0.373 0.008];     % orange (limiter active)
s.inactiveLight=[0.62 0.76 0.89];s.activeLight=[0.99 0.78 0.58];
s.fixed=[0.45 0.45 0.45];         % fixed (H-D independent) boundaries
s.feasible=[0.40 0.65 0.52];      % numerically feasible
s.nofeasible=[0.80 0.80 0.80];    % no feasible rectangle found
s.inconclusive=[0.93 0.69 0.13];  % solver inconclusive / failed
s.noteval=[1 1 1];
s.font='Arial';s.fs=11;s.fsSmall=9;s.fsTitle=12;
s.note='Numerically grid-checked (301\times301) — not interval-certified';
end
