% plots figure 2 for TTW frontal arrest paper and some additional results of interest

close all
addpath("functions\")

% code sections to run

calc_sol = false;
plot_sol = false;
plot_u_balance = false;
plot_sol_vs_TTW = true;

% parameters for particular case studied

Nx = 64;
Nz = 51;

E = 0.1;
Ro = 1.6;
Pr = 1;

b = @(x) erf(x*sqrt(pi)/2);
dbdx = @(x) exp(-(x*sqrt(pi)/2).^2);

err_tol = 1e-5;
err_max = 10;

% get solution for parameters defined above

if calc_sol

    [U, V, W, B, P, x, z, ~, ~, ~, ~, dx, dz] = TTW_State(10, E, Ro, Pr, 0, Nx, Nz, 1, 2, dbdx, err_tol, err_max, 1);

    d_dx = @(F) dx * F;
    d_dz = @(F) (dz * (F'))';

    Iz = dz;
    Iz(1, :) = [1 zeros(1, Nz-1)];
    Iz = Iz^-1 * diag([0 ones(1, Nz-1)]);
    I_z = @(F) (Iz * (F'))';

    Lx = -2*x(1);

end

% plot U, V, W, B and Psi for particular case

if plot_sol

    Psi = I_z(U);     % calculate streamfunction by integrating U in z      

    splot(U, x, z, xlabel = 'x', ylabel = 'z', xlim = [-3 3])
    hold on; contour(x, z, U', [-0.15 -0.1 -0.05 0 0.05 0.1 0.15], 'k'); hold off

    splot(V, x, z, xlabel = 'x', ylabel = 'z', xlim = [-3 3])
    hold on; contour(x, z, V', [-0.15 -0.1 -0.05 0 0.05 0.1 0.15], 'k'); hold off

    splot(W, x, z, xlabel = 'x', ylabel = 'z', xlim = [-3 3])
    hold on; contour(x, z, W', [-0.04 -0.02 0 0.02 0.04], 'k'); hold off

    splot(B, x, z, xlabel = 'x', ylabel = 'z', xlim = [-3 3]); clim([-1.0001 1.0001])
    hold on; contour(x, z, B', [-0.8 -0.6 -0.4 -0.2 0 0.2 0.4 0.6 0.8], 'k'); hold off

    splot(Psi, x, z, xlabel = 'x', ylabel = 'z', xlim = [-3 3]); colormap(cmap(centre = -0.000001))
    hold on; contour(x, z, Psi', [0.02 0.04 0.06], 'k'); hold off
    
end


% Look at surface fields to identify dominant terms in balances

if plot_u_balance

    % calculate terms in cross-front momentum equation

    Adv_b = Ro*(U .* (d_dx(B-2*x/Lx) + 2/Lx) + W .* d_dz(B));
    Diff_b = E/Pr*d_dz(d_dz(B));
    
    Adv_u1 = Ro * U .* d_dx(U);
    Adv_u2 = Ro * W .* d_dz(U);   % zero on surface due to w = 0 BC
    Diff_u = E*d_dz(d_dz(U));
    dPdx = d_dx(P) + dbdx(x)*(z');
    
    splot(Adv_u1 + Adv_u2 - V - Diff_u + dPdx, x, z, xlabel = 'x', ylabel = 'z')   % equation balance, will equal subtracted depth-averaged jet tendency
    
    % plot terms in cross-front momentum equation on the top surface (z = 0.5)

    figure;
    plot(x, Adv_u1(:,end), x, -V(:, end), x, -Diff_u(:, end), x, dPdx(:, end), 'LineWidth', 2)                                                              
    xlabel('x'); legend('$Ro\, u\, \partial u/\partial x$', '$-v$', '$-E\, \partial^2 u/\partial z^2$', '$\partial p/\partial x$', interpreter = 'latex'); grid
    set(gca,'FontSize',12,'linewidth',0.7); xlim([-3 3])

end

% Plot cross-front velocity for one case, compare full result minus TTW
% velocity with order Ro correction

if plot_sol_vs_TTW

    structure_functions     % load definitions of z0, K, K_p, K_pp, P1, and P2
    
    b0 = b(x);
    b0_x = dbdx(x);
    b0_xx = d_dx(b0_x);
    
    u0 = -sqrt(E) * (b0_x * K_pp(z'/sqrt(E), z0));
    v0 = -sqrt(E) * (b0_x * K(z'/sqrt(E), z0));

    u1 = E * (2*(b0_x .* b0_xx) * P1(z'/sqrt(E)));
    v1 = E * (2*(b0_x .* b0_xx) * P2(z'/sqrt(E)));
    b1 = -sqrt(E) * (b0_x .^2 * K(z'/sqrt(E), z0));

    b1_x = -2 * sqrt(E) * (b0_x .* b0_xx * K(z'/sqrt(E), z0));
    
    Bx = d_dx(B-2*x/Lx) + 2/Lx;     % numerically differentiate B, subtract background gradient first to make B periodic

    splot(U - u0 - Ro * u1, x, z, xlabel = 'x', ylabel = 'z')       % compare u between converged solution and TTW O(Ro) prediction
    splot(V - v0 - Ro * v1, x, z, xlabel = 'x', ylabel = 'z')       % compare u between converged solution and TTW O(Ro) prediction
    splot(B - b0 - Ro * b1, x, z, xlabel = 'x', ylabel = 'z')       % compare b between converged solution and TTW O(Ro) prediction
    splot(Bx - b0_x - Ro * b1_x, x, z, xlabel = 'x', ylabel = 'z')  % compare b_x between converged solution and TTW O(Ro) prediction

end

