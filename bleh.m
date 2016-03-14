syms x z u w theta q u1 u2;
m = 0.435;
g = 9.807;
Iyy = 0.01;
d = 0.25;

f = [
    u * cos(theta) + w * sin(theta);
    - u * sin(theta) + w * cos(theta);
    - g * sin(theta) - q * w;
    -2/m * (u1 + u2) + g * cos(theta) + q * u;
    q;
    d * (u1 - u2) / Iyy;
    ];

dfdx = jacobian(f, [x, z, u, w, theta, q]);
dfdu = jacobian(f, [u1, u2]);

dfdx_0 = subs(dfdx, [u w theta q], [0 0 0 0]);
dfdu_0 = subs(dfdu, [u w theta q], [0 0 0 0]);

F = double(dfdx_0);
G = double(dfdu_0);
sys = ss(F, G, eye(6), zeros(6, 2));

Q = diag([1 1 10 1 1000 1]);
R = 100 * eye(2);
K = lqr(sys, Q, R);

disp('Problem 1')

C = K
[V,D] = eig(F-G*C)

disp('Problem 2')

T = 10;
x0 = [0 10 0 0 0 0]';
sys_lti = ss(F - G * C, [], [], []);
[y, t, x] = initial(sys_lti, x0, T);

%{
subplot(4, 2, 1);
plot(t, x(:,1));
xlabel('t, s');
ylabel('x (m)');

subplot(4, 2, 2);
plot(t, x(:,2));
xlabel('t, s');
ylabel('z (m)');

subplot(4, 2, 3);
plot(t, x(:,3));
xlabel('t, s');
ylabel('u (m/s)');

subplot(4, 2, 4);
plot(t, x(:,4));
xlabel('t, s');
ylabel('w (m/s)');

subplot(4, 2, 5);
plot(t, x(:,5));
xlabel('t, s');
ylabel('\theta (rad)');

subplot(4, 2, 6);
plot(t, x(:,6));
xlabel('t, s');
ylabel('q (rad/s)');

u = zeros(size(x,1), 2);
for i=1:size(x,1)
    [u(i, 1), u(i, 2)] = controls(t(i), x(i,:)');
end
    
subplot(4, 2, 7);
plot(t, u(:,1));
xlabel('t, s');
ylabel('u_1 (N)');

subplot(4, 2, 8);
plot(t, u(:,2));
xlabel('t, s');
ylabel('u_2 (N)');



figure;
%}

disp('Problem 3')

%odeset('RelTol',1e-4,'AbsTol',[1e-5 1e-5 1e-5 1e-5 1e-5 1e-5]);
[T,Y] = ode45(@quadro, [0 T], x0);


subplot(4, 2, 1);
plot(T, Y(:,1));
xlabel('t, s');
ylabel('x (m)');

subplot(4, 2, 2);
plot(T, Y(:,2));
xlabel('t, s');
ylabel('z (m)');

subplot(4, 2, 3);
plot(T, Y(:,3));
xlabel('t, s');
ylabel('u (m/s)');

subplot(4, 2, 4);
plot(T, Y(:,4));
xlabel('t, s');
ylabel('w (m/s)');

subplot(4, 2, 5);
plot(T, Y(:,5));
xlabel('t, s');
ylabel('\theta (rad)');

subplot(4, 2, 6);
plot(T, Y(:,6));
xlabel('t, s');
ylabel('q (rad/s)');

for i=1:size(Y,1)
    [u(i, 1), u(i, 2)] = controls(T(i), Y(i,:)');
end
    
subplot(4, 2, 7);
plot(T, u(:,1));
xlabel('t, s');
ylabel('u_1 (N)');

subplot(4, 2, 8);
plot(T, u(:,2));
xlabel('t, s');
ylabel('u_2 (N)');


%figure;
%plot(T, Y);

disp('Problem 4')


