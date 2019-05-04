%% Example 9.23: Pathwise series expansion of Benes SDE
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
    dt = 1;
    x0 = 0.5;
    t  = 5;
    
    trans_dens = @(xx,t) 1./sqrt(2*pi*t).*cosh(xx)./cosh(x0).*exp(-0.5*t).*exp(-1./(2*t).*(xx-x0).^2);

    xx = -15:0.1:15;
    plot(xx,trans_dens(xx,t));

%%
% Expansion approximation
%
  % Lock random seed
  if exist('rng') % Octave doesn't have rng
    rng(1);
  else
    randn('state',2);
  end

    NN = [1:20 50 100]; % Series components
    nmc = 10000;
    xsamp = zeros(nmc,length(NN));
    
    trajx = cell(1,length(NN));
    trajt = cell(1,length(NN));
    
    opt = odeset;
    if exist('OCTAVE_VERSION', 'builtin') ~= 0
        % Octave integrates in quite long steps by default
        opt.MaxStep = 0.05;
    end

    for i=1:nmc
        z0 = randn(100,1);
        for j=1:length(NN)
            N = NN(j);
            z = z0(1:N);
            T = t;
            ode_f = @(t,x) tanh(x) ...
                + z' * cos((2*(1:N)'-1)*pi*t/T/2) * sqrt(2/T);

            Tspan = [0 T];
            [T,X] = ode45(ode_f,Tspan,x0,opt);
%            fprintf('Number of steps = %d.\n',length(T));

            xsamp(i,j) = X(end);

            if isempty(trajx{j})
                trajx{j} = X;
                trajt{j} = T;
            end
        end
        if rem(i,100) == 0
            fprintf('%d/%d\n',i,nmc);
        end
    end
    
%%
% Visualization
%
    figure(1); clf; hold on
    
      % Exact density
      fill(xx,trans_dens(xx,t),1,'FaceColor',[.7 .7 .7],'EdgeColor',[.7 .7 .7])

      % Approximation    
      bins = linspace(-15,15,100);
      i = find(NN == 10);
      nh = hist(xsamp(:,i),bins);
      stairs(bins-(bins(2)-bins(1))/2,nh/sum(nh)/(bins(2)-bins(1)),'-k');
    
      xlabel('$x$');
      ylabel('$p(x)$');
      xlim([-15 15]);
      ylim([0 0.15]);

    figure(2); clf; hold on

      % Exact density
      fill(xx,trans_dens(xx,t),1,'FaceColor',[.7 .7 .7],'EdgeColor',[.7 .7 .7])

      % Approximation     
      bins = linspace(-15,15,100);
      i = find(NN == 100);
      nh = hist(xsamp(:,i),bins);
      stairs(bins-(bins(2)-bins(1))/2,nh/sum(nh)/(bins(2)-bins(1)),'-k');
    
      xlabel('$x$');
      ylabel('$p(x)$');
      xlim([-15 15]);
      ylim([0 0.15]);
    
    figure(3); clf;
      i = find(NN == 10);
      j = find(NN == 100);
      h = plot(trajt{i},trajx{i},trajt{j},trajx{j});
      set(h(1),'Color',[0.8 0.8 0.8]);
      set(h(1),'LineWidth',1.5);
      set(h(2),'Color',[0.0 0.0 0.0]);
      set(h(2),'LineWidth',.5);
      xlabel('Time, $t$');
      ylabel('$x(t)$');
     
    figure(4); clf;
      m_ta = x0 + tanh(x0) * t;
      P_ta = x0^2 + 2 * x0 * tanh(x0) * t + t + t^2 - m_ta^2;
      ml = mean(xsamp);
      vl = std(xsamp);
      h = semilogx(NN,abs(ml-m_ta*ones(size(NN))),'-',...
                 NN,abs(vl-sqrt(P_ta)*ones(size(NN))),'--');
      set(h(1),'Color',[0.0 0.0 0.0]);
      set(h(1),'LineWidth',1);
      set(h(2),'Color',[0.0 0.0 0.0]);
      set(h(2),'LineWidth',1);
      xlabel('Number of terms, $N$');
      ylabel('Absolute error');
      legend('Error in mean','Error in STD');
      xlim([1 100]);
    
    