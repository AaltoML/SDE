%% Example 2.9: Numerical solution of ODEs
%
% Copyright: 
%   2018 - Simo Särkkä and Arno Solin
%
% License:
%   This software is provided under the MIT License. See the accompanying 
%   LICENSE file for details.

%% Simulate model

  % Specify method (add your method here)
  method_name = {'Euler','Heun','RK4'};
  method_handle = {@euler,@heun,@rk4simple};
  lineStyles = {'-','--','-.'};
  
  % Set parameters
  t1 = 10;
  dt = .1;
  g = 1;
  v = 2;
  q = 0.02;
  x0 = [1;0];

  % Dynamics
  F = [0 1; -v^2 -g];
  f = @(x,t) F*x;
  
  % Exact solution
  t_exact = linspace(0,t1,500);
  x_exact = zeros(size(F,1),numel(t_exact));
  for j=1:numel(t_exact)
    x_exact(:,j) = expm(F*t_exact(j))*x0;
  end
    
  % Allocate space for approximate results
  t = cell(length(method_name),1);
  x = cell(length(method_name),1);
  
  % Run each method
  for j=1:length(method_name)
    [x{j},t{j}] = feval(method_handle{j},f,[0:dt:t1]',x0);
  end  
  
  % Show result
  figure(1); clf; hold on
  
  % Plot exact
  plot(t_exact,x_exact(1,:),'-','LineWidth',2,'Color',[.5 .5 .5])
  
  % Plot numerical results
  for j=1:length(method_name)
    plot(t{j},x{j}(1,:),lineStyles{j},'Color','k')
  end

  % Pimp up the plot
  legend(['Exact' method_name])
  box on
  xlabel('Time, $t$')
  ylabel('$x_1$')
  ylim([-.8 1])
  


%% Error plots

  % Step sizes
  dt = logspace(-3,-1,8);
  
  % Errors
  err = zeros(length(method_name),numel(dt));
  
  % For each step length
  for i=1:numel(dt)
   
    % Time steps
    t = 0:dt(i):t1;
      
    % Exact result
    x_exact = zeros(size(F,1),numel(t));
    for j=1:numel(t)
      x_exact(:,j) = expm(F*t(j))*x0;
    end

    % Run each method
    for j=1:length(method_name)
      
      % Evaluate
      x = feval(method_handle{j},f,t,x0);
      
      % Calcualte absolute error
      err(j,i) = mean(abs(x(:)-x_exact(:)));
      
    end
  end

  figure(2); clf; hold on
    for j=1:length(method_name)
      plot(dt,err(j,:),lineStyles{j},'Color','k')
    end
    plot(dt,err,'x','Color','k')
    set(gca,'XScale','log','YScale','log')
    box on
    xlabel('$\Delta t$')
    ylabel('Absolute error, $|\hat{\vx} - \vx|$')
    set(gca,'YTick',10.^[-12 -8 -4 0], ...
            'YTickLabel',{'$10^{-12}$','$10^{-8}$','$10^{-4}$','$10^{0}$'})
    