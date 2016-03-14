function [u1, u2] = controls(t, y)

% taken from bleh.m
%C = [-0.7071   -0.7071   -1.0458   -0.8086    4.1164    0.8153;
%    0.7071   -0.7071    1.0458   -0.8086   -4.1164   -0.8153];

C = [-0.0707   -0.0707   -0.2965   -0.1428    2.6271    0.3318;
    0.0707   -0.0707    0.2965   -0.1428   -2.6271   -0.3318];


% nominal controls, derived in homework #4

u_star = 1.0665;
u_vec = [u_star; u_star] - C * y;
u1 = u_vec(1);
u2 = u_vec(2);
