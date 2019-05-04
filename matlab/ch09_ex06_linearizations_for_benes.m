%% Examples 9.6 and 9.14: Linearizations and approximations for the Benes model
%
% Copyright: 
%   2018 - Simo Särkkä and Arno Solin
%
% License:
%   This software is provided under the MIT License. See the accompanying 
%   LICENSE file for details.

%%
% Reference result
%

    dt = 1;
    T  = (0:dt:5);
    x0 = 0.5;
    t  = 5;
    
    q = 1; % Change this from one only if you know what you are doing
    
    trans_dens = @(xx,t) 1./sqrt(2*pi*t).*cosh(xx)./cosh(x0).*exp(-0.5*t).*exp(-1./(2*t).*(xx-x0).^2);

    % diff(tanh(x))
    % = 1 - tanh(x)^2
    % 
    % diff(tanh(x),2)
    % = 2*tanh(x)*(tanh(x)^2 - 1)
    
    % Lock random seed
    if exist('rng') % Octave doesn't have rng
      rng(1);
    else
      randn('state',2);
    end

    C = chol([dt^3/3 dt^2/2; dt^2/2 dt],'lower');
    XX1 = zeros(1,10000);
    XX2 = zeros(size(XX1));
    XX3 = zeros(size(XX1));
    for n=1:size(XX1,2)
        x = x0;
        y = x0;
        z = x0;
        for k=1:length(T)-1
            % Euler
            x = x + tanh(x) * dt + sqrt(dt)*sqrt(q)*randn;

            % IT1.5
            zb = C * sqrt(q) * randn(2,1);
            fy = tanh(y);
            dfy = 1 - tanh(y)^2;
            d2fy = 2*tanh(y)*(tanh(y)^2 - 1);
            a = dfy * fy + 0.5 * d2fy;
            b = dfy;
            y = y + fy * dt + zb(2) + a * dt^2/2 + b * zb(1);

            % SRK1.5
            xp = z + tanh(z) * dt + sqrt(dt);
            xm = z + tanh(z) * dt - sqrt(dt);
            z = z + tanh(z) * dt + zb(2) + ...
                1/4*(tanh(xp) - 2*tanh(z) + tanh(xm)) * dt + ...
                1/(2*sqrt(dt)) * (tanh(xp) - tanh(xm)) * zb(1);
        end
        XX1(n) = x;
        XX2(n) = y;
        XX3(n) = z;
    end
    
    xx = -15:0.1:15;
    
    bins = linspace(-15,15,100);
    
    subplot(3,2,1);
    nh = hist(XX1,bins);
    hold off;
    bar(bins,nh/sum(nh)/(bins(2)-bins(1)));
    hold on;
    plot(xx,trans_dens(xx,t));
    title('Euler');
    xlim([xx(1) xx(end)]);

    subplot(3,2,2);
    nh = hist(XX2,bins);
    hold off;
    bar(bins,nh/sum(nh)/(bins(2)-bins(1)));
    hold on;
    plot(xx,trans_dens(xx,t));
    title('Ito-Taylor 1.5');
    xlim([xx(1) xx(end)]);

    subplot(3,2,3);
    nh = hist(XX3,bins);
    hold off;
    bar(bins,nh/sum(nh)/(bins(2)-bins(1)));
    hold on;
    plot(xx,trans_dens(xx,t));
    title('SRK 1.5');
    xlim([xx(1) xx(end)]);
    
    
    subplot(3,2,6);
    plot(xx,trans_dens(xx,t));
    title('Exact');
    xlim([xx(1) xx(end)]);
    
    sum(trans_dens(xx,t) .* (xx(2)-xx(1)))

    
%%
% GH weights
%
    n = 10;
    
    %
    % Form Probabilists' Hermite polynomials of
    % order n-1 and n
    %
    Hpm = 1;
    Hp  = [1 0];
    for i=1:n-1
        tmp = Hp;
        Hp = [Hp 0] - [0 0 i*Hpm];
        Hpm = tmp;
    end

    %
    % Single dimensional weights and points
    %
    xi = roots(Hp)';
    W = factorial(n)./(n^2*polyval(Hpm,xi).^2);
    

%%
% Gaussian integration exactly, using linearization
% and sigma-points.
%
    dt = 0.01;
    t  = 5;
    T = (0:dt:t);

    m_ga = x0;
    P_ga = eps;
    
    m_li = x0;
    P_li = 0;

    m_si = x0;
    P_si = 0;

    dx = 0.01;
    grid = -20:dx:20;
    
    for k=1:length(T)
        gd_ga = 1/sqrt(2*pi*P_ga)*exp(-1/(2*P_ga)*(grid - m_ga).^2);
        gd_ga = gd_ga ./ sum(gd_ga) / dx;

        dm_ga = sum(tanh(grid) .* gd_ga .* dx);
        dP_ga = sum(2*(1-tanh(grid).^2).*gd_ga*dx)*P_ga + 1;
        m_ga = m_ga + dm_ga * dt;
        P_ga = P_ga + dP_ga * dt;

        dm_li = tanh(m_li);
        dP_li = 1 + 2*(1 - tanh(m_li)^2)*P_li;
        m_li = m_li + dm_li * dt;
        P_li = P_li + dP_li * dt;
        
        dm_si = sum(W .* tanh(m_si + sqrt(P_si) * xi));
        dP_si = 1 + 2 * sum(W .* sqrt(P_si) .* xi .* tanh(m_si + sqrt(P_si) * xi));
        m_si = m_si + dm_si * dt;
        P_si = P_si + dP_si * dt;
        
    end

    gd_ex = trans_dens(grid,t);
    gd_ex = gd_ex ./ sum(gd_ex) / dx;
    m_ex = sum(grid .* gd_ex .* dx);
    P_ex = sum(grid.^2 .* gd_ex .* dx) - m_ex^2;
            
    m_ta = x0 + tanh(x0) * t;
    P_ta = x0^2 + 2 * x0 * tanh(x0) * t + t + t^2 - m_ta^2;

    gd_ta = 1/sqrt(2*pi*P_ta)*exp(-1/(2*P_ta)*(grid - m_ta).^2);
    gd_ta = gd_ta ./ sum(gd_ta) / dx;
    
    gd_li = 1/sqrt(2*pi*P_li)*exp(-1/(2*P_li)*(grid - m_li).^2);
    gd_li = gd_li ./ sum(gd_li) / dx;
    
    gd_si = 1/sqrt(2*pi*P_si)*exp(-1/(2*P_si)*(grid - m_si).^2);
    gd_si = gd_si ./ sum(gd_si) / dx;
    
    xx = -15:0.1:15;
    figure(1); clf;
    plot(xx,trans_dens(xx,t),grid,gd_ta,grid,gd_ga,'--',grid,gd_si,'-.',grid,gd_li,':');
    legend('Exact density','Taylor','Grid Gaussian fit','Gauss-Hermite','Linearization','Location','northwest');

    
%%
% Test sigma-point methods for simulation
%

    dt = 1;
    t  = 5;
    T = (0:dt:t);
    eu_steps = 10;

    % Lock random seed
    if exist('rng') % Octave doesn't have rng
      rng(1);
    else
      randn('state',2);
    end
    
    XX3 = zeros(1,10000);
    for n=1:size(XX1,2)
        x = x0;
        for k=1:length(T)-1
            m_si = x;
            P_si = 0;

            ddt = dt/eu_steps;
            for i=1:eu_steps
                dm_si = sum(W .* tanh(m_si + sqrt(P_si) * xi));
                dP_si = 1 + 2 * sum(W .* sqrt(P_si) .* xi .* tanh(m_si + sqrt(P_si) * xi));
                m_si = m_si + dm_si * ddt;
                P_si = P_si + dP_si * ddt;
            end

            x = m_si + sqrt(P_si) * randn;
        end
        XX3(n) = x;
    end

    figure(1); clf;
    nh = hist(XX3,bins);
    hold off;
    bar(bins,nh/sum(nh)/(bins(2)-bins(1)));
    hold on;
    plot(xx,trans_dens(xx,t));

%%
% Test the simulation with local linearizations
%

    dt = 1;
    t  = 5;
    T = (0:dt:t);

    % Lock random seed
    if exist('rng') % Octave doesn't have rng
      rng(1);
    else
      randn('state',2);
    end
    
    XX4 = zeros(1,10000);

    for n=1:size(XX1,2)
        x1 = x0;
        x2 = x0;

        for k=1:length(T)
            x  = x1;
            f  = tanh(x);
            df = 1 - tanh(x)^2;
            d2f = 2*tanh(x)*(tanh(x)^2 - 1);
            F = 1/dt*log(1 + f/(x*df)*(exp(df*dt) - 1));

            m1 = exp(F*dt)*x;
            P1 = 1/(2*F)*(exp(2*F*dt)-1);

            x  = x2;
            f  = tanh(x);
            df = 1 - tanh(x)^2;
            d2f = 2*tanh(x)*(tanh(x)^2 - 1);
            
            G = df;
            a = d2f/2;
            b = f - x*df;

            %
            % dm/dt = G m + a t + b
            % dP/dt = 2 G P + 1
            %
            % m(t) = exp(G*t)*x0 + int_0^t exp(G (t-s)) [a s + b] ds
            %      = (b*(exp(G*t) - 1))/G - a*(t/G - (exp(G*t) - 1)/G^2)
            % P(t) = 1/(2G) [exp(2 G dt) - 1]

            m2 = exp(G*dt)*x + (b*(exp(G*dt) - 1))/G - a*(dt/G - (exp(G*dt) - 1)/G^2);
            P2 = 1/(2*G)*(exp(2*G*dt)-1);
            
            x1 = m1 + sqrt(P1) * randn;
            x2 = m2 + sqrt(P2) * randn;            
        end
        XX4(n) = x1;
        XX5(n) = x2;
    end
    
    figure(1); clf
    subplot(1,2,1);
    nh = hist(XX4,bins);
    hold off;
    bar(bins,nh/sum(nh)/(bins(2)-bins(1)));
    hold on;
    plot(xx,trans_dens(xx,t));
    
    subplot(1,2,2);
    nh = hist(XX5,bins);
    hold off;
    bar(bins,nh/sum(nh)/(bins(2)-bins(1)));
    hold on;
    plot(xx,trans_dens(xx,t));
    
%%
% Plot the actual figures
%   

    % Gaussian approximations
    figure(1); clf; hold on
    
      % Exact density
      fill(xx,trans_dens(xx,t),1,'FaceColor',[.7 .7 .7],'EdgeColor',[.7 .7 .7])
    
      % The lines
      h = plot(grid,gd_ta,'-',grid,gd_si,'--',grid,gd_li,'-.');
      set(h(1),'Color',[0 0 0]);
      set(h(1),'LineWidth',1);
      set(h(2),'Color',[0 0 0]);
      set(h(2),'LineWidth',1);
      set(h(3),'Color',[0 0 0]);
      set(h(3),'LineWidth',1);

      legend('Exact density','Exact Gaussian','Gauss--Hermite','Linearization','Location','northwest');
      xlabel('$x$');
      xlim([-15 15]);
      ylim([0 0.2]);
    
    % Approximations of the histogram

    % Sigma-point
    figure(2); clf; hold on

      % Exact density
      fill(xx,trans_dens(xx,t),1,'FaceColor',[.7 .7 .7],'EdgeColor',[.7 .7 .7])

      % Approximation
      nh = hist(XX3,bins);
      stairs(bins-(bins(2)-bins(1))/2,nh/sum(nh)/(bins(2)-bins(1)),'-k');
      
      xlabel('$x$');
      ylabel('$p(x)$');
      xlim([-15 15]);
      ylim([0 0.15]);
    
    % Ozaki
    figure(3); clf; hold on

      % Exact density
      fill(xx,trans_dens(xx,t),1,'FaceColor',[.7 .7 .7],'EdgeColor',[.7 .7 .7])

      % Approximation
      nh = hist(XX4,bins);
      stairs(bins-(bins(2)-bins(1))/2,nh/sum(nh)/(bins(2)-bins(1)),'-k');

      xlabel('$x$');
      ylabel('$p(x)$');
      xlim([-15 15]);
      ylim([0 0.15]);

    % Shoji
    figure(4); clf; hold on

      % Exact density
      fill(xx,trans_dens(xx,t),1,'FaceColor',[.7 .7 .7],'EdgeColor',[.7 .7 .7])

      % Approximation    
      nh = hist(XX5,bins);
      stairs(bins-(bins(2)-bins(1))/2,nh/sum(nh)/(bins(2)-bins(1)),'-k');
      
      xlabel('$x$');
      ylabel('$p(x)$');
      xlim([-15 15]);
      ylim([0 0.15]);

    % Ito-Taylor
    figure(5); clf; hold on

      % Exact density
      fill(xx,trans_dens(xx,t),1,'FaceColor',[.7 .7 .7],'EdgeColor',[.7 .7 .7])

      % Approximation
      nh = hist(XX2,bins);
      stairs(bins-(bins(2)-bins(1))/2,nh/sum(nh)/(bins(2)-bins(1)),'-k');
      
      xlabel('$x$');
      ylabel('$p(x)$');
      xlim([-15 15]);
      ylim([0 0.15]);

