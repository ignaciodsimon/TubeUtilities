% Data obtained from the 12AX7 curves using the TubeUtilities program. To
% be used in the circuit emulation.
%
% Joe.

Gm = [-0.50, 1.91
      -1.00, 1.82
      -1.50, 1.68
      -2.00, 1.52
      -2.50, 1.32
      -3.00, 1.12
      -3.50, 0.85
      -4.00, 0.55];   % <-- Originally 0.65 mA

Ri = [ 0.0, 49000.820275714643 % <-- Originally 46299 kOhm
      -0.5, 52157.554473869692
      -1.0, 54285.675148822549
      -1.5, 58108.341103889412
      -2.0, 64772.066932879468
      -2.5, 74548.305694974246
      -3.0, 84624.833194823615
      -3.5, 105193.50339705426
      -4.0, 144579.43919013219
      -4.5, 206860.51243538884];

plot(Gm(:,1), Gm(:,2), 'x-');
hold on
plot(Ri(:,1), Ri(:,2) / 100000, 'x-');
xlim([-5 0.5])
ylim([0 2.5])
xlabel('Vgrid [V]', 'FontSize', 12);
ylabel(sprintf('Ri [ohm] x 10^5\nGm [mA/V]'), 'FontSize', 12)
set(gca, 'FontSize', 12);
legend({'Ri', 'Gm'}, 'FontSize', 12)
grid
