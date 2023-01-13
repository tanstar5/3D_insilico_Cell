c_p = 0.3;
c_0 = 0;
c_n = -0.3;
c_expected = .3;

total_area = 0.04;
temperature = 3;

k = 10000;
mixing_factor = 1;
hamiltonian = @(x,y) 0.5.*k*total_area*((c_expected - (c_p.*x  + c_n.*y + c_0*(1 - x -y))).^2);
% hamiltonian = @(x,y) 0.5*1*( (.3 - 1*x).^2 + (.3 - 1*y).^2 + (.3 - 1*(1 - x - y)).^2);
entropy = @(x,y) - mixing_factor*temperature*(x.*log(x) + y.*log(y) + (1-x-y).*log(1-x-y));
[X_grid,Y_grid]  = meshgrid((0:0.01:1),(0:0.01:1));

% entropy_2d = @(x) -(x.*log(x) + (1-x).*log(1-x));
hamiltonian_2d = @(x) 0.5*1*(.3 - (.3*x  +  -0.3*(1 - x)  ).^2);

energy_well = (-hamiltonian(X_grid,Y_grid));
probability = exp(-hamiltonian(X_grid,Y_grid)/temperature);
entropy_well = real(entropy(X_grid,Y_grid));

equillibrium_landscape = exp(entropy_well).*probability;

surf(X_grid,Y_grid,equillibrium_landscape);daspect([1,1,1]);
% surf(X_grid,Y_grid,entropy_well);hold off
xlabel('positive curvature');ylabel('negative curvature'); view(2)

[~,ind] = max(equillibrium_landscape(:));

x_fraction = X_grid(ind)
y_fraction = Y_grid(ind)

c_equillibrium = x_fraction*c_p + y_fraction*c_n + (1-x_fraction-y_fraction)*c_0


% plot(0:0.01:1,entropy_2d(0:0.01:1));hold on;
% plot(0:0.01:1,hamiltonian_2d(0:0.01:1));hold off;
