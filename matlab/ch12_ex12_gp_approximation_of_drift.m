%% Example 12.12: GP approximation of drift functions (double-well)
%
% Copyright: 
%   2018 - Simo Särkkä and Arno Solin
%
% License:
%   This software is provided under the MIT License. See the accompanying 
%   LICENSE file for details.

%%

  % Lock random seed
  if exist('rng') % Octave doesn't have rng
    rng(2,'twister')
  else
    randn('state',2);
  end

  % The double-well model
  f = @(x,t) 4*(x-x.^3);
  L = @(x,t) 1;
  q = 1;
  x0 = 0;
  
  dt = 0.01;
  t = 0:dt:(20000-1)*dt;
  
  % Simulate
  [x,tspan] = eulermaruyama(f,L,t,x0,q);
  
  % Observed
  ty = t(1:20:end);
  y = x(1:20:end);
  
  %ty = ty(1:1000);
  %y = y(1:1000);
  Dt = ty(2)-ty(1);
  
  % Show
  figure(1); clf; hold on
    plot(t,x,'-','Color',[.5 .5 .5])
    plot(ty,y,'-k')

  % Plot the true f
  figure(2); clf
  subplot(211)
    z = linspace(-2.5,2.5,100);
    plot(z,f(z,0))
    xlim([-2.5 2.5]); ylim([-5 5])
  subplot(212)
    hist(y,z)
    xlim([-2.5 2.5]);
    
  % Define GP model
  sigma_lin2 = 10^2;
  sigma_se2 = 20^2;
  ell = 1;
  k = @(x,y) sigma_lin2*x(:)*y(:)' + sigma_se2*exp(-(x(:)-y(:)').^2/2/ell^2);
  
  % GP prediction
  d = (y(2:end)-y(1:end-1))'/Dt;
  gpf_mean = @(x) k(x,y(1:end-1)) * ((k(y(1:end-1),y(1:end-1))+q/Dt*eye(numel(y)-1))\d);
  gpf_cov = @(x) k(x,x) - k(x,y(1:end-1))*((k(y(1:end-1),y(1:end-1))+q/Dt*eye(numel(y)-1))\k(y(1:end-1),x));
    
  figure(3); clf; hold on
    z = linspace(-2.5,2.5,100);
 
    % Predict
    m = gpf_mean(z);
    v = diag(gpf_cov(z));
    
    % Plot mean
    h1 = plot(z,m,'-k');
        
    % Plot uncertainty
    h2 = fill([z'; flipud(z')],[m+1.96*sqrt(v); flipud(m-1.96*sqrt(v))],1);
    set(h2,'EdgeColor',[.7 .7 .7],'FaceColor',[.7 .7 .7])
    
    % Plot mean
    plot(z,m,'-k');
    
    % Plot ground-truth
    h3 = plot(z,f(z,0),'--k');
    
    % Limits
    box on
    xlim([-2.5 2.5]); ylim([-5 5])
    xlabel('Input, $x$')
    ylabel('$f(x)$')
    
    % Ticks
    set(gca,'Layer','top')
    set(gca,'XTick',-2:2,'YTick',-4:2:4)
    
    % Legend
    legend([h1 h2 h3],'Mean','95\% quantiles','Ground-truth drift')
    
  figure(4); clf
    plot(ty,y,'-k')
    
    % Limits
    box on
    ylim([-2.5 2.5]);
    xlabel('Time, $t$')
    ylabel('$x(t)$')
    
    % Ticks
    set(gca,'Layer','top')
    set(gca,'YTick',-2:2)
    
    