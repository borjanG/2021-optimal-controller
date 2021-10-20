t = 0:pi/200:10*pi;
r = 1;
n = 2;

xk = rand(n,1);

res = zeros(1, length(t));
for ii = 1:length(t)
    res(ii) = fun_heat(r.*sin(t(ii)), r.*cos(t(ii)), xk);
end

plot3(r.*sin(t), r.*cos(t), zeros(1, length(t)), ['--', 'k']);
hold on;
plot3(r.*sin(t), r.*cos(t), res, 'r', 'LineWidth', 2);
ax = gca;
ax.XGrid = 'on';
ax.YGrid = 'on';
ax.ZGrid = 'on';
zlim([0 0.25]);
exportgraphics(ax,'heat.pdf','ContentType','vector')
max(max(res))

%% Power Iteration
% This is the power iteration method which can be found on Wikipedia.
% We put many iterations to be sure that we have convergence. For n big,
% more iterations are requried. 
function r = power_iteration(A, xk)
    max_iter = 1000;
    for i=1:max_iter
        x_ = A*xk;
        xk = x_/norm(x_);
    end
    r = (transpose(xk)*(A*xk))/norm(xk)^2;
end

function r = p_(k)
  
    %% Laplacian
    n = 2;
    h = 1/(n+1);
    A = -2*eye(n) + diag(ones(n-1,1),-1) + diag(ones(n-1,1),1);
    A = 1/(h^2)*A;
    
    a_ = charpoly(A);
    if k==n
        r = eye(n);
    else
        aux = zeros(n,n);
        for j=1:(n-k)
            aux = aux + a_(j+1)*mpower(A,n-k-j);
        end
        r=(mpower(A, n-k)+aux);
    end
end 


%% Heat
function r = fun_heat(b1, b2, xk)

    mat = zeros(2,2);
    b = [b1; b2];
    for k=1:2
       mat(:, k) = p_(k)*b; 
    end
    
    %% Power iteration to find \lambda_min(P(b)P(b)^*)
    C = mat*transpose(mat);
    r_ = eigs(C);
    r = r_(length(r_));
%     lbda_max = power_iteration(C, xk);
%     D = C-lbda_max*eye(n);
%     
%     r = (power_iteration(D, xk)+lbda_max);
end


