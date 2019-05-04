%% Example 7.12: Doob's h-transform
%
% Copyright: 
%   2018 - Simo Särkkä and Arno Solin
%
% License:
%   This software is provided under the MIT License. See the accompanying 
%   LICENSE file for details.

%% For the OU model

  % Lock random seed
  if exist('rng') % Octave doesn't have rng
      rng(0,'twister');
  else
      randn('state',2);
  end

  % Parameters
  lambda = 1;
  q = 1;
  x0 = 0;
  T = 1;
  xT = 5;

  % Define the derived functions
  a = @(t) exp(-lambda*(T-t));
  sigma2 = @(t) q/(2*lambda) * (1-exp(-2*lambda*(T-t)));

  % Specify the model to simulate from
  f = @(x,t) -lambda*x + q*a(t)/sigma2(t)*(xT - a(t)*x);
  L = @(x,t) 1;
  
  % Simualte trajectories
  t = linspace(0,T,100);
  n = 100;
  
  % Allocate space
  x = zeros(n,numel(t));
  
  % Simulate using Euler-Maruyama
  for i=1:n
    x(i,:) = eulermaruyama(f,L,t,x0,q);
  end
  
  figure(1); clf; hold on
  
    % Plot sample trajectories  
    plot(t,x','-','Color',[.5 .5 .5],'LineWidth',0.1)
  
    % Plot dashed line showing x(T)
    plot([0 T],xT*[1 1],'--k')
    
    box on
    xlabel('Time, $t$'); ylabel('$x(t)$')
    set(gca,'Layer','Top')
    set(gca,'YTick',[0:2:4 xT],'YTickLabel',{'$0$','$2$','$4$','$x_T$'})
    
    
    