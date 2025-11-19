


% Gauss-Laguerre Quadrature (generalized)


alpha = 4.5; % power inside the integral

%inside_f = @(x) 2.*x.^1; 
inside_f = @(x) 2.*x.^(1); %
f = @(x)  inside_f(x) .* x.^(alpha) .* exp(-x);


% Define the exact integral value using MATLAB's integral function
exact_integral = integral(f, 0, Inf); %

% Trapezoidal Rule
a = 0;
b = 20;
n_p = 10; % Number of intervals
x_trap = linspace(a, b, n_p+1);
y_trap = f(x_trap);
trap_integral = (b - a) / (2 * n_p) * (y_trap(1) + 2 * sum(y_trap(2:end-1)) + y_trap(end));

% Gauss-Laguerre Quadrature (1 point: x_1)
m = 1;
% P2 is the number 2 generalized Laguerre polynomial with parameter alpha
P2 = @(x) 0.5*(x^2 - 2*(alpha+2)*x + (alpha+1)*(alpha+2));
% x_1 is a root of P2 
x_1 = alpha+1;


x_GL = [alpha+1];

% Weights are given in terms of the generalized Laguerre polynomials
w_GL = [gamma(m+alpha+1)*x_1/(factorial(m)*(m+1)^2*P2(x_1)^2)];
gauss_L_integral = sum(w_GL .* inside_f(x_GL));

% Display the results
disp(['Exact Integral: ', num2str(exact_integral)]);
disp(['Trapezoidal Rule Integral: ', num2str(trap_integral)]);
disp(['Gauss-Laguerre Quadrature Integral: ', num2str(gauss_L_integral)]);

% Plot the function and approximations
x = linspace(a, b, 100);
y = f(x);

figure;
h1 = plot(x, y, 'k-', 'LineWidth', 1.5); hold on;

% Plot trapezoidal approximation
for i = 1:n_p
    h2=fill([x_trap(i) x_trap(i) x_trap(i+1) x_trap(i+1)], [0 f(x_trap(i)) f(x_trap(i+1)) 0], 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
end

% Plot Gaussian quadrature points
h3 = plot(x_GL, f(x_GL), 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b');

title('Comparison of Trapezoidal Rule and Gauss-Laguerre Quadrature');
xlabel('x');
ylabel('f(x)');
legend([h1, h2,h3], {'f(x)', 'Trapezoid', 'Gaussian Laguerre Quadrature Points'}, 'Location', 'NorthWest');
grid on;
hold off;


% exact for all degree of polynom < 2*m-1


