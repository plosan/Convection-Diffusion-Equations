clear;
close all;
clc;

%% 1. DATA
L = 1;                          % Domain length     [m]
u = @(t,x,y) 2*y*(1-x.^2);      % Velocity field x component
v = @(t,x,y) -2*x*(1-y.^2);     % Velocity field y component

n = 1e6;    % Max size of vectors
h = 1e-3;   % Time step             [s]

x = zeros(1, n);    % Curve x component     [m]
y = zeros(1, n);    % Curve y component     [m]
t = zeros(1, n);    % Time                  [s]

x0 = -0.5*L;        % x initial condition   [m]
y0 = 0;             % y initial condition   [m]
x(1) = x0;
y(1) = y0;

%% 2. RUNGE-KUTTA
i = 1;
finish = false;

while (i < n && finish == false)
    % Coefficients 1
    i1 = u(t(i), x(i), y(i));
    j1 = v(t(i), x(i), y(i));
    % Coefficients 2
    i2 = u(t(i)+0.5*h, x(i)+0.5*h*i1, y(i)+0.5*h*j1);
    j2 = v(t(i)+0.5*h, x(i)+0.5*h*i1, y(i)+0.5*h*j1);
    % Coefficients 3
    i3 = u(t(i)+0.5*h, x(i)+0.5*h*i2, y(i)+0.5*h*j2);
    j3 = v(t(i)+0.5*h, x(i)+0.5*h*i2, y(i)+0.5*h*j2);
    % Coefficients 4
    i4 = u(t(i)+h, x(i)+h*i3, y(i)+h*j3);
    j4 = v(t(i)+h, x(i)+h*i3, y(i)+h*j3);
    % Next step
    t(i+1) = t(i) + h;
    x(i+1) = x(i) + h*(i1 + 2*i2 + 2*i3 + i4)/6;
    y(i+1) = y(i) + h*(j1 + 2*j2 + 2*j3 + j4)/6;
    % End
    i = i + 1;
    if (y(i) < 0 || y(i) > L)
        finish = true;
    end
    if (abs(x(i)) > L)
        finish = true;
    end
end

if finish
    fprintf("i = %d\n", i);
    x(i+1:end) = [];
    y(i+1:end) = [];
    t(i+1:end) = [];
end

%% 3. SAVE
filename = sprintf("output/case_smith_hutton/smith_hutton_streamline_%.2f.dat", abs(x0));

fileID = fopen(filename, 'w');
for k = 1:i
    fprintf(fileID, "%.5f %.5f\n", x(k), y(k));
end
fclose(fileID);

filename = sprintf("output/case_smith_hutton/smith_hutton_streamline_arrow_%.2f.dat", abs(x0));
fileID = fopen(filename, 'w');

for k = 1:i
    posx = x(k);
    posy = y(k);
    velx = u(0, posx, posy);
    vely = v(0, posx, posy);
    vel = sqrt(velx*velx + vely*vely);
    velx = velx / vel;
    vely = vely / vel;
    fprintf(fileID, "%.5f %.5f %.5f %.5f\n", posx, posy, velx, vely);
end
fclose(fileID);

%% 4. PLOT
theta = linspace(0, pi, 1e3);
x_circ = abs(x0) * cos(theta);
y_circ = abs(x0) * sin(theta);

figure(1);
hold on;
plot(x, y, 'b');
plot(x_circ, y_circ, 'r');
plot([-L L], [0 0], 'k', 'LineWidth', 1);
plot([-L L], [L L], 'k', 'LineWidth', 1);
plot([L L], [0 L], 'k', 'LineWidth', 1);
plot([-L -L], [0 L], 'k', 'LineWidth', 1);
xlim([-L L]);
ylim([0 L]);
pbaspect([2 1 1]);
legend("Streamline", "Circle");
hold off;


