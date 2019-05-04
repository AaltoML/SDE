%% Example 8.15: Duffing van der Pol oscillator
%
% Copyright: 
%   2018 - Simo Särkkä and Arno Solin
%
% License:
%   This software is provided under the MIT License. See the accompanying 
%   LICENSE file for details.

%% Duffing van der Pol

  % Time-span
  tspan = 0:2^-5:20;

  % Parameters
  alpha = 1;

  % Define arrow (for visalization)
  arrow1 = [-1 1 0 -1; -.5 -.5 2 -.5]';
    
  % The model
  f = @(x,t) [x(2,:); x(1,:).*(alpha - x(1,:).^2)-x(2,:)];
  L = @(x,t) [zeros(1,size(x,2)); x(1,:)];

  
%% ODE  
  
  figure(1); clf; hold on
  
  for j=1:10
    %x = rk4(f,tspan,[-2-.2*j; 0]);
    x = rk4simple(f,tspan,[-2-.2*j; 0]);
    %[~,x] = ode45(@(t,x) f(x,t),tspan,[-2-.2*j; 0]); x = x';
    
    % Plot trajectory
    plot(x(1,:),x(2,:),'-k','LineWidth',.25)
    
    % Plot direction
    uv = f(x,[]); ind = 10;
    newquiver(x(1,ind),x(2,ind),uv(2,ind),-uv(1,ind), ...
      'X',arrow1,'scale',.04*[1 14.8/8.8])
    
  end
  
  % Set axis limits
  axis([-4.4 4.4 -4.8 10]), axis square
  %set(gca,'XTick',-4:2:4,'YTick',-4:2:8)
  xlabel('$x_1$'); ylabel('$x_2$')
  box on

%% SDE  
    
  figure(2); clf; hold on
  figure(3); clf; hold on
  
  for j=1:10
      
    figure(2);  
      
    % Lock random seed
    if exist('rng') % Octave doesn't have rng
      rng(3,'twister')  
    else
      randn('state',2);
      rand('state',2);
    end
      
    % Use the strong order 1.0 method
    x = srkS10scalarnoise(f,L,tspan,[-2-.2*j; 0],.5^2);
    
    % Plot trajectory
    plot(x(1,:),x(2,:),'-k','LineWidth',.25)
    
    % Plot direction
    uv = f(x,[]); ind = 10;
    newquiver(x(1,ind),x(2,ind),uv(2,ind),-uv(1,ind), ...
      'X',arrow1,'scale',.04*[1 14.8/8.8])
    
    figure(3);
    plot(tspan,x(1,:),'-k')
    plot(tspan,x(2,:),'-','Color',[.7 .7 .7])
  
  end
  
  % Set axis limits
  figure(2)
  axis([-4.4 4.4 -4.8 10]), axis square
  %set(gca,'XTick',-4:2:4,'YTick',-4:2:8)
  xlabel('$x_1$'); ylabel('$x_2$')
  box on
  
  % Set axis limits
  figure(3)
  xlim([0 20])
  xlabel('Time, $t$'); ylabel('$x$')
  legend('$x_1(t)$','$x_2(t)$')
  box on
  
  
%% Weak approximation

  % Time-span
  tspan = 0:2^-4:20;

  % Reset random seed
  % Lock random seed
  if exist('rng') % Octave doesn't have rng
    rng(1,'twister')  
  else
    randn('state',2);
    rand('state',2);
  end

  % Initial point
  x0 = [-3;0];

  % Number of smaples
  x = zeros(2,10000);
  
  % Iterate
  for j=1:size(x,2)
          
    % Weak SRK scheme
    foo = srkW20(f,L,tspan,x0,.5^2,true);    
    
    % Store
    x(:,j) = foo(:,end);
    
    % Report
    if rem(j,100)==0,
      figure(4); clf
        hist(x(1,1:j),ceil(sqrt(j)))
        drawnow
      j
    end
    
  end
  
  
%% Histogram
  
  % Bins for histogram
  t = linspace(min(x(1,:)),max(x(1,:)),64);
  
  figure(5); clf

    % Show solution w2.0
    n = histc(x(1,:),t);
    stairs(t-(t(2)-t(1))/2,n/size(x,2),'-k')

    % Label
    xlabel('$x_1$')
    
    % Ticks
    box off
    xlim([-2.2 2.2])
    %set(gca,'XTick',0:.2:1.2)
    
  figure(6); clf

    % Show solution w2.0
    plot(x(1,:),x(2,:),'.k')

    % Label
    xlabel('$x_1$')
    ylabel('$x_2$')
    
    % Ticks
    box on
    