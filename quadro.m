function dy = quadro(t, y)
d = 0.25;
Iyy = 0.01;
maxT = 2;
minT = 0;
m = 0.435;
g = 9.807;

x = y(1);
z = y(2);
u = y(3);
w = y(4);
theta = y(5);
q = y(6);

[u1, u2] = controls(t, y);
X = 0;
Z = - 2 * (u1 + u2);
M = d * (u1 - u2);

dy = zeros(6, 1);
dy(1) = cos(theta) * u + sin(theta) * w;
dy(2) = -sin(theta) * u + cos(theta) * w;
dy(3) = - g * sin(theta) - q * w;
dy(4) = Z / m + g * cos(theta) + q * u;
dy(5) = q;
dy(6) = M / Iyy;
