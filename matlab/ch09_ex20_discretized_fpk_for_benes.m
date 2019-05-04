%% Example 9.20: Discretized FPK for the Benes SDE
%
% Copyright: 
%   2018 - Simo Särkkä and Arno Solin
%
% License:
%   This software is provided under the MIT License. See the accompanying 
%   LICENSE file for details.

%%
% This is the true density
%
    x0 = 0.5;
    t  = 5;
    
    trans_dens = @(xx,t) 1./sqrt(2*pi*t).*cosh(xx)./cosh(x0).*exp(-0.5*t).*exp(-1./(2*t).*(xx-x0).^2);

    xx = -15:0.1:15;
    plot(xx,trans_dens(xx,t));

%%
% Dirichlet (Neumann) orthogonal basis
%
    L = 15;

    lam = @(i) i*pi/(2*L);
    phi = @(x,i) 1/sqrt(L)*sin(lam(i)*(x + L));
    
    % Check orthonormality
    M = zeros(4);
    for i=1:size(M,1)
        for j=1:size(M,2)
            M(i,j) = integral(@(x) phi(x,i) .* phi(x,j),-L,L);
        end
    end
    M
    
%%
% Form the delta function with the series
%
    %  int phi_i(x) delta(x - x0) dx
    %  = phi_i(x0)
    %
    N = 20;
        
    c0 = phi(x0,1:N)';
    
    xx = -15:0.1:15;

    del = zeros(size(xx));
    for i=1:N
        del = del + c0(i) * phi(xx,i);
    end
    plot(xx,del);
    
%%
% Form the matrix F
%
    F = zeros(N,N);
    
    for i=1:N
        for j=1:N
            F(i,j) = 1/L*integral(@(x) ((tanh(x).^2-1) .* sin(lam(j).*(x + L)) - tanh(x) .* lam(j) .* cos(lam(j) .* (x + L))) .* sin(lam(i).*(x + L)),-L,L);
        end
    end
    F = F - 0.5*diag(lam(1:N).^2);
    pcolor(F);
    
%%
% Solve the equation using backward Euler
%
    dt = 0.1;
    T = (0:dt:t);
    xx = -15:0.1:15;
    
    c = c0;
    for k=1:length(T)
        c = (eye(size(F)) - F*dt) \ c;
        f = zeros(size(xx));
        for i=1:N
            f = f + c(i) * phi(xx,i);
        end
    end
    plot(xx,trans_dens(xx,t),'--',xx,f)
    
%%
% Solve using matrix exponential
%
    c = expm(F*t)*c0;
    
%%
% Plot the final figure
%
    figure(1); clf; hold on
    
      % Exact density
      fill(xx,trans_dens(xx,t),1,'FaceColor',[.7 .7 .7],'EdgeColor',[.7 .7 .7])
      
      % Series approximation
      plot(xx,f,'-k','LineWidth',1)
    
      ylim([0 0.15])      
      legend('Exact density','Series approximation','Location','northwest');
      xlabel('$x$');
      ylabel('$p(x)$');
    
    
    
    