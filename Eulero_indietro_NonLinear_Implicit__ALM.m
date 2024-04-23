clear all
close all
%function penalty()

% Number of discrete points
%n = 1000;
a=-50;
b=50;
h=0.2;
n=  abs(b-a)/h;
%h= abs(b-a)/n;
dt=0.1;
T=20;
t=0:dt:T;
% Spatial dimension of the problem is 1
mesh = linspace(a, b, n)';
%u0 =zeros(size(mesh)); %condizione iniziale
%u0 =[0; ones(size(mesh) - [2 0]); 0]; %condizione iniziale

% Reaction-diff
% u0 = 0.05 * exp ( -5 * mesh.^2 );
% f= 0*t;
% reaction_term = 1;

% Diff
u0 = zeros(size(mesh));
f = -ones(size(u0));
reaction_term = 0;


A_l=zeros(n,n);
A_l= A_l - 2/h^2 * diag(ones(n,1),0) + 1/h^2 * diag(ones(n-1,1),-1) + 1/h^2 * diag(ones(n-1,1),1);
K=diag(ones(size(A_l(:,1))),0);
%eye - identità
K= speye(n);
A_l=sparse(A_l);



u=u0;

%vincolo u < ub
mu = 1e4;
lambda=zeros (n,1);

ub = 0.5 * ones(n,1);
c = @(x) x - ub;
cg = @(x) eye (n,n) ; %hessiana del vincolo è identità
ch = zeros(n,n); %per vincolo lineare
d = @(x) max ( x - ub , 0 );
epq = @(x)  dot ( d(x),d(x) ) /2;
epg = @(x)  d(x);
eph = @ (x, lambda) spdiags( (mu*c(x)+lambda)>0 , 0, n,n );

%Newton method to solve eq:  u_kk - dt  A * u_kk - dt *  u_kk .* ( 1 - u_kk) - u_k = 0
g = @ (xkk, xk) xkk - dt *  A_l * xkk - dt * reaction_term * (xkk - xkk.^2) -  xk + dt .* f;
h = @ (xkk) K - dt * A_l - dt * (K - 2* sparse(1:n,1:n, xkk,n,n) ) ;
alpha = 1;
niter = 4;
nshift = 60;
s=0;

use_shift = 1;

%With penalty
Qg = @(xkk, xk) g(xkk, xk) + dt * epg(xkk);
Qh = @(x, lambda) h(x) + dt * eph(x, lambda);

%Shifted-penalty
lambda=zeros (n,1);
Lg = @(xkk,xk,lambda) g(xkk,xk) +  cg(xkk) * (mu*c(xkk)+lambda) .* ((mu*c(xkk)+lambda)>0) ;
%cg(xkk) * max(mu * c(xkk) + lambda,0)
Lh = @(x, lambda) h(x) + mu * cg(x) .* eph(x, lambda);


close all;
hold off;

i=1;
violazioni = zeros(length(0:dt:T),1);

for t=0:dt:T
    
    uk = u;
    ukk = 0.2*ones(n,1);
    lambda =zeros(n,1);
    active = u  > ub;
    violazioni (i) = sum(active);
    % while norm(ub-u,2)<=toll;
    q=0;
    for j=1:nshift
        q=q+1;
       
        for k=1:niter            
            pk = - sparse(Lh(ukk, lambda)) \ Lg(ukk, uk, lambda);
            ukk = ukk + alpha * pk;            
        end

        norm_Lg = norm(Lg(ukk, uk,lambda), 2);
        norm_d= norm(d(ukk));
        
        lambda = lambda + c(ukk);
        active = u  > ub;
        lambda= use_shift * active.* lambda;

        %active = ukk -(ub+lambda)  >= 0;
        %lambda =  active .* (ukk - ub - lambda);

        
        if (norm_Lg < 1e-8 && norm_d < 1e-10) || j==nshift
            norm_lambda= norm(lambda,2);
            disp(sprintf('%d %g %g %g', j,norm_Lg,norm_lambda, norm_d));
            break;
    end
        
    end

  

    u = ukk;
    u(1) = 0;
    u(end) = 0;
    
    
    if mod(i, 100) == 1
        nplots = 4;
        subplot(nplots, 1, 1);
        plot(mesh, u);
        % hold on;


        subplot(nplots, 1, 2);
        plot(mesh,lambda)
        %ylim([-1e-3, 1e-3]);
        
        % hold on
        % hold off;
      

        %         subplot(nplots,1,2);
        %         plot(1:n,active, 'o');
        %         hold on
        
        %plot(mesh , f(t));
        
        
        pause(0.05);
    end
    i=i+1;
end
%plot (mesh, ub);