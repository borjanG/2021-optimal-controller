tic;

n = 3;

xk = rand(n,1);

%% sphere code
figure(1);
N = 70; % controls the number of faces on the sphere
res = zeros(N+1, N+1);
[X,Y,Z] = sphere(N); % returns a sphere with N x N faces
                     % hence with N+1 x N+1 vertices
for ii=1:N+1
    for jj=1:N+1
        res(ii, jj) =  fun_heat(X(ii, jj),Y(ii, jj), Z(ii, jj), xk);
    end
end

s = surf(X, Y, Z, res); 
%s.EdgeColor = 'interp';
hold on

axis equal      % set the axis ratio so the sphere appear as a sphere

disp(max(max(res)))

ax = gca;
ax.XGrid = 'on';
ax.YGrid = 'on';
ax.ZGrid = 'on';
%light
%lightangle(260,-45)
lighting gouraud % preferred lighting for a curved surface
%lightangle(-15,20)
lightangle(260,-45)
%lightangle(-15,20)
light
s.FaceLighting = 'gouraud';
s.AmbientStrength = 0.3;
s.DiffuseStrength = 0.8;
s.SpecularStrength = 0.9;
s.SpecularExponent = 25;
%s.alpha = 0.85;
set(s,'edgecolor','none')
%s.BackFaceLighting = 'unlit';
colorbar('southoutside')
shading interp   
exportgraphics(ax,'convect_3d.pdf','ContentType','vector')
toc;

%% Power Iteration
% This is the power iteration method which can be found on Wikipedia.
% We put many iterations to be sure that we have convergence. For n big,
% more iterations are requried. 
function r = power_iteration(A, xk)
    max_iter = 100000;
    for i=1:max_iter
        x_ = A*xk;
        xk = x_/norm(x_);
    end
    r = (transpose(xk)*(A*xk))/norm(xk)^2;
end

function r = p_(k)
    %% Laplacian
    n = 3;
    h = 1/(n+1);
    A_ = -2*eye(n) + diag(ones(n-1,1),-1) + diag(ones(n-1,1),1);
    A_ = 1/(h^2)*A_;
    A = 1/(2*h)*(zeros(n,n)+diag(ones(n-1,1),-1)-diag(ones(n-1,1),1))+A_;
    
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
function r = fun_heat(b1, b2, b3, xk)

    b = [b1; b2; b3];
    %% Laplacian
    n = length(b);
    h = 1/(n+1);
    A_ = -2*eye(n) + diag(ones(n-1,1),-1) + diag(ones(n-1,1),1);
    A_ = 1/(h^2)*A_;
    A = 1/(2*h)*(zeros(n,n)-diag(ones(n-1,1),-1)+diag(ones(n-1,1),1))+A_;
    
    %% Constructing P(b):
    a_ = charpoly(A);
    
    mat = zeros(n,n);
    for k=1:n
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


