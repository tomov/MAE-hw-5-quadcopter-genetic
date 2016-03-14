% stuff from 1, 2, 3

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

F = double(dfdx_0);                % USED IN GA
G = double(dfdu_0);                % USED IN GA
Q = diag([1 1 10 1 1000 1]);
R = 100 * eye(2);
N = zeros(6, 2);
sys = ss(F, G, eye(6), N);
[K, S, e] = lqr(sys, Q, R, N);
C = K;

T = 10;                            % USED IN GA
x0 = [10 0 0 0 0 0]';              % USED IN GA
sys_lti = ss(F - G * C, [], [], []);
[y, t, x] = initial(sys_lti, x0, T);

x_optimal = x;                     % USED IN GA


% stuff from 4
% lots of variables are overwritten lol...

syms c11 c12 c13 c14 c15 c16;
bits = 10;
N = 512;
m = 6;
max_c = 4;
gens = 300;
mutation_rate = 1000;

C_template = [-c11 -c12 -c13 -c14 c15  c16;
              c11  -c12 c13  -c14 -c15 -c16];

conv0 = [2^9 2^8 2^7 2^6 2^5 2^4 2^3 2^2 2^1 2^0 zeros(1, (m - 1) * bits)];
conv = conv0;
for i = 1:m-1
    conv = circshift(conv, [0 bits]);
    conv = vertcat(conv0, conv);
end

 
population = rand(N, bits * m) > 0.5;
Ftotal = zeros(1, gens);
Fbest = Ftotal;
xs = zeros(gens, size(x_optimal, 1), size(x_optimal, 2));

for generation=1:gens
    generation
    % calculate fitness of population
    fitness = zeros(N, 1);
    for i = 1:N
        
        % get C
        chrom = population(i, :);
        c = conv * chrom';
        c = max_c * c / 2^bits;
        C = subs(C_template, [c11 c12 c13 c14 c15 c16], c');
        C = double(C);
        
        % calculate x_GA
        sys_lti = ss(F - G * C, [], [], []);
        [y, t, x] = initial(sys_lti, x0, T);
        
                
        % calculate fitness
        J = 0;
        for j = 1:size(x,1)
            J = J + (x_optimal(j,:) - x(j,:)) * (x_optimal(j,:) - x(j,:))';
        end
        J = J / 2;
        fit = 1 / J;
        if fit > max(fitness)
            xs(generation,:,:) = x;
        end
        fitness(i) = fit;
    end
    
    s = sum(fitness)
    Ftotal(generation) = s;
    Fbest(generation) = max(fitness);
    
    % roulette wheel to get new population
    new_population = zeros(N, bits * m);
    for i = 1:N
        r = rand;
        cumsum = 0;
        for j = 1:N
            cumsum = cumsum + fitness(j) / s;
            if cumsum > r
                new_population(i,:) = population(j,:);
                break;
            end
        end
    end
    
    % crossover
    for k = 1:N/2
        i = k*2 - 1;
        j = k*2;
        pos = int32(ceil(rand * (bits * m - 1)));
        left = new_population(i,:);
        right = new_population(j,:);
        new_population(i,:) = [left(1:pos) right(pos+1:end)];
        new_population(j,:) = [right(1:pos) left(pos+1:end)];
    end
    
    % mutation
    for i = 1:size(population, 1)
        for j = 1:size(population, 2)
            if rand < 1/mutation_rate
                disp('FLIP!');
                population(i, j) = ~population(i, j);
            end
        end
    end
        
    population = new_population;
end

x = squeeze(xs(end,:,:));
plot(t, x);