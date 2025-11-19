% Define the function

%f_inter = @(x) 2;
f_inter = @(x) 2.*x.^2;
f = @(x) f_inter(x) .* exp(-x.^2);

% Define the exact integral value using MATLAB's integral function
exact_integral = integral(f, -Inf, Inf); %2*sqrt(pi)

% Trapezoidal Rule
a = -3;
b = 3;
n = 10; % Number of intervals
x_trap = linspace(a, b, n+1);
y_trap = f(x_trap);
trap_integral = (b - a) / (2 * n) * (y_trap(1) + 2 * sum(y_trap(2:end-1)) + y_trap(end));


% Gauss-Hermite Quadrature (6 points)
% x_hermite are the roots of the physicists' version of the Hermite polynomial

H_hermite_5 = [32 0 -160 0 120 0];
H_hermite_6 = [64 0 -480 0 720 0 -120];

% x_hermite = [-2.350604973674492, -1.335849074013697, -0.436077411927617, 0.436077411927617, 1.335849074013697, 2.350604973674492];
% w_hermite = [0.004530009905509, 0.157067320322857, 0.724629595224392, 0.724629595224392, 0.157067320322857, 0.004530009905509];

% Gauss-Hermite Quadrature (10 points)
%H_hermite_9 =
%H_hermite_10 = [1024 0 -23040 0 161280 0 -403200 0 302400 0 -30240];

% x_hermite = [-3.436159118837738, -2.532731674232789, -1.756683649299881, -1.036610829789513, -0.342901327223705, ...
%               0.342901327223705,  1.036610829789513,  1.756683649299881,  2.532731674232789,  3.436159118837738];
% w_hermite = [7.640432855232620e-06, 0.001343645746781232, 0.033874394455481, 0.240138611082314, 0.610862633735325, ...
%              0.610862633735325, 0.240138611082314, 0.033874394455481, 0.001343645746781232, 7.640432855232620e-06];

H_hermite_n = H_hermite_6; 
H_hermite_n_minus = H_hermite_5;
%H_hermite = H_hermite_10;
x_hermite = roots(H_hermite_n)';

nbr_points = size(x_hermite,2); % n

% Computation of the Gauss–Hermite quadrature weights (according to the
% number of points
w_hermite = [];
for i=1:nbr_points
    
    w_hermite_i = 2^(nbr_points-1)*factorial(nbr_points)*sqrt(pi)/(nbr_points^2*polyval(H_hermite_n_minus,x_hermite(:,i))^2);
    w_hermite = [w_hermite,w_hermite_i];
end
%w_hermite

gauss_hermite_integral = sum(w_hermite .* f_inter(x_hermite));

% Display the results
disp(['Exact Integral: ', num2str(exact_integral)]);
disp(['Trapezoidal Rule Integral: ', num2str(trap_integral)]);
disp(['Gauss-Hermite Quadrature Integral: ', num2str(gauss_hermite_integral)]);

% Plot the function and approximations
x = linspace(a, b, 100);
y = f(x);

figure;
h1 = plot(x, y, 'k-', 'LineWidth', 1.5); hold on;

% Plot trapezoidal approximation
for i = 1:n
    h2=fill([x_trap(i) x_trap(i) x_trap(i+1) x_trap(i+1)], [0 f(x_trap(i)) f(x_trap(i+1)) 0], 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
end

% Plot Gaussian quadrature points
h3 = plot(x_hermite, f(x_hermite), 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b');

title('Comparison of Trapezoidal Rule and Gauss-Hermite Quadrature');
xlabel('x');
ylabel('f(x)'); % 2x^2e^{-x^2}
legend([h1, h2,h3], {'f(x) = 2x^2e^{-x^2}', 'Trapezoid integral approximation', 'Gaussian Quadrature Points'}, 'Location', 'NorthWest');
grid on;
hold off;
