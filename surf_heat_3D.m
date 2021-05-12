clc
clear

% % x = -1:0.05:1;  % define range and mesh of x and y which will be shown in figure
% % y = -1:0.05:1;

% % t = 0:pi/20:10*pi;
% % r = 1;
n = 3;

%xk = rand(n,1);
xk = [0.6551; 0.1626; 0.1190];

% % disp(xk);

% % res = zeros(length(x), length(y));
% % for ii=1:length(x)
% %     for jj=1:length(y)
% %         res(ii, jj) = fun_wave(x(ii),y(jj), xk);
% %     end
% % end

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

s = surf(X, Y, Z, res); hold on;

axis equal;      % set the axis ratio so the sphere appear as a sphere

disp(max(max(res)))

ax = gca;
ax.XGrid = 'on';
ax.YGrid = 'on';
ax.ZGrid = 'on';
lighting gouraud % preferred lighting for a curved surface
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
exportgraphics(ax,'heat_3d.pdf','ContentType','vector')

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


%% Heat
function r = fun_heat(b1, b2, b3, xk)

    b = [b1; b2; b3];
    %% Laplacian
    n = length(b);
    h = 1/(n-1);
    A = -2*eye(n) + diag(ones(n-1,1),-1) + diag(ones(n-1,1),1);
    A = 1/(h^2)*A;
    
    %% Constructing P(b):
    a_ = charpoly(A);
    
    mat = zeros(n,n);
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
    
    r = (power_iteration(D, xk)+lbda_max);
end


