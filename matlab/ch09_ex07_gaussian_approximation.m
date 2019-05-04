%% Example 9.7: Gaussian approximation of a nonlinear trigonometric SDE
%
% Copyright: 
%   2018 - Simo Särkkä and Arno Solin
%
% License:
%   This software is provided under the MIT License. See the accompanying 
%   LICENSE file for details.

%% Gaussian approximation

  % Lock random seed
  if exist('rng') % Octave doesn't have rng
    rng(1,'twister') 
  else
    randn('state',2);
  end

  normpdf = @(x,m,s) 1./(s .* sqrt(2*pi)) .* exp(-0.5*(x-m).^2 ./ s.^2);

  % Parameters
  tspan = 0:2^-6:10;

  % The model
  f = @(x,t) -(1/10)^2*sin(x).*cos(x).^3;
  L = @(x,t) 1/10*cos(x).^2;

  % Exact solution
  solfun = @(wt,x0) atan(1/10*wt +tan(x0));

  % Intial
  x0 = 1;

  % Sigma-points (1d cubature rule)
  xi = [-1 1];
  wi = [1/2 1/2];
 
  % The state is (m, P)
  fun = @(z,t) [sum(wi.*f(z(1)+sqrt(z(2))*xi,t));
                2*sum(wi.*f(z(1)+sqrt(z(2))*xi,t).*xi*sqrt(z(2,:)))+...
                sum(wi.*L(z(1)+sqrt(z(2))*xi,t).^2)];

  % Integrate using RK4
  [mP] = rk4simple(fun,tspan,[x0;0]);

  % Samples from the exact solution
  x = solfun(sqrt(tspan(end))*randn(1,200000),x0);  
  
  % Bins
  nbins = 64;
  t = linspace(min(x),max(x)+.1,nbins);
  
  % Histogram based ground truth
  n = histc(x,t);
  
  figure(1); clf; hold on
  
    % Show exact
    fill(t,n/numel(x),1,'FaceColor',[.7 .7 .7],'EdgeColor',[.7 .7 .7])
    
    % Exact Gaussian approximation
    ti = linspace(min(x),max(x)+.1,2*nbins);
    npdf = normpdf(ti,mean(x),std(x));
    plot(ti,2*npdf/sum(npdf),'-','Color',[0 0 0],'LineWidth',1)
    
    % Approximation based solution
    npdf = normpdf(ti,mP(1,end),sqrt(mP(2,end)));
    plot(ti,2*npdf/sum(npdf),'--','Color',[.5 .5 .5],'LineWidth',1)
    
    % Label
    xlabel('$x$')
    
    % Limits
    xlim([min(t) max(t)])
    
    % Show legend
    legend('Exact','Exact Gaussian', ...
        'Approximate Gaussian')

    % Set limit to match the latter plots
    xlim([.35 1.4])
    ylim([0 .12])
    