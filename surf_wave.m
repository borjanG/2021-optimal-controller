clc
clear

% % x = -1:0.05:1;  % define range and mesh of x and y which will be shown in figure
% % y = -1:0.05:1;

t = 0:pi/200:10*pi;
r = 1;
n = 2;

%xk = rand(2*n,1);
xk = [0.7900; 0.3185; 0.5341; 0.0900];
% % disp(xk);

% % res = zeros(length(x), length(y));
% % for ii=1:length(x)
% %     for jj=1:length(y)
% %         res(ii, jj) = fun_wave(x(ii),y(jj), xk);
% %     end
% % end

res = zeros(1, length(t));
for ii = 1:length(t)
    res(ii) = fun_wave(r.*sin(t(ii)), r.*cos(t(ii)), xk);
end

plot3(r.*sin(t), r.*cos(t), zeros(1, length(t)), ['--', 'k']);
hold on;
plot3(r.*sin(t), r.*cos(t), res, 'b', 'LineWidth', 2);
ax = gca;
ax.XGrid = 'on';
ax.YGrid = 'on';
ax.ZGrid = 'on';
exportgraphics(ax,'wave.pdf','ContentType','vector')
max(max(res))

% % surf(x,y, res);
% % figure;
% % contourf(x,y, res);




%% Power Iteration
% This is the power iteration method which can be found on Wikipedia.
% We put many iterations to be sure that we have convergence. For n big,
% more iterations are requried. 
function r = power_iteration(A, xk)
    max_iter = 10000;
    for i=1:max_iter
        x_ = A*xk;
        xk = x_/norm(x_);
    end
    r = (transpose(xk)*(A*xk))/norm(xk)^2;
end


function r = fun_wave(b1, b2, xk)
    
    %%  Wave operator
    b = [b1; b2];
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

    bb = zeros(n, 1);
    bb(n_+1:n, 1) = b;
    mat = zeros(n,n);
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
    r_ = eigs(C);
    r = r_(length(r_));
%     lbda_max = power_iteration(C, xk);
%     D = C-lbda_max*eye(n);
%     
%     r = (power_iteration(D, xk)+lbda_max);
end


