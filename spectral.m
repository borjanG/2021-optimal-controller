clear
%clc
import casadi.*

%% Symbols, expressions
n=3;

% We define the symbolic variable representing the optimal solution $b$
b = SX.sym('b',n);

% We define sybolically the constraint ($|b|^2-1=0$)
constraint = norm(b)^2-1;
tc = constraint(:);

% heat = 0, wave = 1, convect = 2
model = 2; 
if model == 0
    % We initialize the power iteration algorithm with a random vector 
    % (this is just to approximate \lambda_1); f represents the
    % objective, g represents the constraint
    xk = rand(n,1);
    nlp = struct('x', b(:), 'f', fun_heat(b, xk), 'g', tc(:));
elseif model == 1
    xk = rand(2*n,1);
    nlp = struct('x', b(:), 'f', fun_wave(b, xk), 'g', tc(:));
else 
    xk = rand(n,1);
    nlp = struct('x', b(:), 'f', fun_convect(b, xk), 'g', tc(:));
end

%disp(xk);

% We solve the problem using ipopt
S = casadi.nlpsol('S', 'ipopt', nlp);

%% Setup
% We solve for the dico containing the relevant information from the
% problem S (we take an initial guess which is a random vector for b0)

r = S('x0',rand(n,1),'lbg',0,'ubg',0);
%r = S('x0', [0.335503; 0.942039], 'lbg',0,'ubg',0);

% We recover the optimal value b* as r.x
b_opt = r.x;
disp(b_opt)

% Display of the values of the functional (to just check convergence of PI)
disp(r.f)
if model == 0
    disp(fun_heat(b_opt, xk))
elseif model == 1
    disp(fun_wave(b_opt,xk))
else
    disp(fun_convect(b_opt,xk))
end

%% Power Iteration
% This is the power iteration method which can be found on Wikipedia.
% We put many iterations to be sure that we have convergence. For n big,
% more iterations are requried. 
function r = power_iteration(A, xk)
    max_iter = 50000;
    for i=1:max_iter
        x_ = A*xk;
        xk = x_/norm(x_);
    end
    r = (transpose(xk)*(A*xk))/norm(xk)^2;
end

%% Heat
function r = fun_heat(b, xk)

    %% Laplacian
    n = length(b);
    h = 1/(n-1);
    A = -2*eye(n) + diag(ones(n-1,1),-1) + diag(ones(n-1,1),1);
    A = 1/(h^2)*A;
    
    %% Constructing P(b):
    a_ = charpoly(A);
    
    mat = zeros(n,n,'casadi.SX');
    aux = zeros(n,n);
    for k=1:n
        if k==n
            mat(:,k)=b;
        else
            for j=1:(n-k)
                aux = aux + a_(j+1)*mpower(A,n-k-j);
            end
            mat(:,k)=(mpower(A, n-k)+aux)*b;
        end
    end
    
    %% Power iteration to find \lambda_min(P(b)P(b)^*)
    C = mat*transpose(mat);
    %r_ = eigs(C);
    %r = r_(length(r_));
    lbda_max = power_iteration(C, xk);
    D = C-lbda_max*eye(n);
    
    r = -(power_iteration(D, xk)+lbda_max);
end

%% Wave

function r = fun_wave(b, xk)
    
    %%  Wave operator
    n = 2*length(b);
    n_ = n/2;
    h = 1/(n_-1);
    A_ = -2*eye(n_) + diag(ones(n_-1,1),-1) + diag(ones(n_-1,1),1);
    A_ = 1/(h^2)*A_;
    A = zeros(n,n);
    A(1:n_,n_+1:n) = eye(n_); 
    A(n_+1:n,1:n_)=A_;
    
    %%  Constructing P(b):
    a_ = charpoly(A);

    bb = zeros(n, 1, 'casadi.SX');
    bb(n_+1:n, 1) = b;
    mat = zeros(n,n,'casadi.SX');
    aux = zeros(n,n);
    for k=1:n
        if k==n
            mat(:,k)=bb;
        else
            for j=1:(n-k)
                aux = aux + a_(j+1)*mpower(A,n-k-j);
            end
            mat(:,k)=(mpower(A, n-k)+aux)*bb;
        end
    end
    
    %%   Power iteration to find \lambda_min(P(b)P(b)^*)
    C = mat*transpose(mat);
    lbda_max = power_iteration(C, xk);
    D = C-lbda_max*eye(n);
    
    r = -(power_iteration(D, xk)+lbda_max);
end

%% Convection

function r = fun_convect(b, xk)
    
    %%  Convect
    n = length(b);
    h = 1/(n-1);
    A_ = -2*eye(n) + diag(ones(n-1,1),-1) + diag(ones(n-1,1),1);
    A_ = 1/(h^2)*A_;
    A = 1/(2*h)*(zeros(n,n)+diag(ones(n-1,1),-1)-diag(ones(n-1,1),1))+A_;
    
    %% Constructing P(b):
    a_ = charpoly(A);
    
    mat = zeros(n,n,'casadi.SX');
    aux = zeros(n,n);
    for k=1:n
        if k==n
            mat(:,k)=b;
        else
            for j=1:(n-k)
                aux = aux + a_(j+1)*mpower(A,n-k-j);
            end
            mat(:,k)=(mpower(A, n-k)+aux)*b;
        end
    end
    
    %% Power iteration to find \lambda_min(P(b)P(b)^*)
    C = mat*transpose(mat);
    lbda_max = power_iteration(C, xk);
    D = C-lbda_max*eye(n);
    
    r = -(power_iteration(D, xk)+lbda_max);
end
