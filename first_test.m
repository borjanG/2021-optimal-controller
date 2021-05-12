import casadi.*

opti = casadi.Opti();

n = 3;
h = 1/(n-1);
A = -2*eye(n) + diag(ones(n-1,1),-1) + diag(ones(n-1,1),1);
A = 1/(h^2)*A;

lambdas = zeros(n, 1);
for k=1:n
    lambdas(k) = -1/h^2*(sin(pi*k/(2*(n+1))))^2;
end

a = [-lambdas(1)-lambdas(2)-lambdas(3),
    lambdas(1)*lambdas(2)+lambdas(2)*lambdas(3)+lambdas(1)*lambdas(3),
    -lambdas(1)*lambdas(2)*lambdas(3)];

B = zeros(n);
B(2:n, 2:n) = eye(n-1);
B(n,:) = -a;
e3 = zeros(n, 1); e3(n) = 1;
%%
% Symbols/expressions
x = opti.variable(3,3);

constraint1 = A*x-x*B;

cost = sum(sum(x.^2)) + 1e-3*sum(constraint1(:).^2);
opti.minimize(cost);


%opti.subject_to(constraint1(:)==0);
opti.subject_to(sum(sum((x*e3).^2))==1);

p_opts = struct('expand', true);
s_opts = struct('max_iter', 10000);

opti.solver('ipopt', p_opts, s_opts);
tic
sol = opti.solve();
toc