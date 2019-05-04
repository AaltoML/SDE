%% Example 8.6: Simulating from a trigonometric nonlinear SDE
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
    randn('state',1);
    rand('state',1);
  end

  % Parameters
  tspan = 0:1:10;

  % The model
  f = @(x,t) -(1/10)^2*sin(x).*cos(x).^3;
  L = @(x,t) 1/10*cos(x).^2;
  
  % The derivatives
  df  = @(x,t) -1/100*cos(x).^2.*(2*cos(2*x)-1);
  ddf = @(x,t) 1/100*(sin(2*x) + 2*sin(4*x));
  dL  = @(x,t) -1/5*sin(x).*cos(x);
  ddL = @(x,t) -1/5*cos(2*x);
  
  % Exact solution
  solfun = @(wt,x0) atan(1/10*wt +tan(x0));

  % Intial
  x0 = 1;

  
%% Sample  

  x = zeros(1,10000);
  xem   = x;
  xw20  = x;
  xw20g = x;
  
  % Lock random seed
  if exist('rng') % Octave doesn't have rng
    rng(1,'twister') 
  else
    randn('state',1);
    rand('state',1);
  end
    
  for j=1:size(x,2)
    
    % Use Euler-Maruyama
    foo = eulermaruyama_weak(f,L,tspan,x0,1);
    xem(:,j) = foo(:,end);
    
    % Report
    if rem(j,100)==0, j, end
  
  end
  
  % Lock random seed
  if exist('rng') % Octave doesn't have rng
    rng(1,'twister') 
  else
    randn('state',1);
    rand('state',1);
  end
  
  for j=1:size(x,2)
    
    % Weak order 2.0
    foo = w20scalar({f,df,ddf},{L,dL,ddL},tspan,x0,1,false);
    xw20(:,j) = foo(:,end);
    
    % Report
    if rem(j,100)==0, j, end
  
  end
  
  % Lock random seed
  if exist('rng') % Octave doesn't have rng
    rng(1,'twister') 
  else
    randn('state',1);
    rand('state',1);
  end
  
  for j=1:size(x,2)

    % Weak order 2.0 (Gaussian increments)
    foo = w20scalar({f,df,ddf},{L,dL,ddL},tspan,x0,1,true);
    xw20g(:,j) = foo(:,end);
    
    % Store
    x(:,j) = foo(:,end);
    
    % Report
    if rem(j,100)==0, j, end
  
  end
   
  
%% Visualize
  
  % Lock random seed
  if exist('rng') % Octave doesn't have rng
    rng(1,'twister') 
  else
    randn('state',1);
    rand('state',1);
  end

  % Samples from the exact solution
  xe = solfun(sqrt(tspan(end))*randn(1,200000),x0);  
  
  % Bins
  nbins = 64;
  t = linspace(min(xe),max(xe)+.1,nbins);
  
  figure(2); clf; hold on

    % Show solution
    n = histc(xe,t);
    fill(t,n/numel(xe),1,'FaceColor',[.7 .7 .7],'EdgeColor',[.7 .7 .7])
  
    % Show solution w2.0
    n = histc(xw20,t);
    stairs(t-(t(2)-t(1))/2,n/numel(xw20),'-k')
    
    % Limits
    xlim([.35 1.25])
    lims = ylim;
    
    % Label
    xlabel('$x$')
    
    % Ticks
    set(gca,'XTick',0:.2:1.2)
    

  figure(3); clf; hold on

    % Show solution
    n = histc(xe,t);
    fill(t,n/numel(xe),1,'FaceColor',[.7 .7 .7],'EdgeColor',[.7 .7 .7])
  
    % Show solution w2.0 (Gaussian increments)
    n = histc(xw20g,t);
    stairs(t-(t(2)-t(1))/2,n/numel(xw20g),'-k')
    
    % Limits
    %xlim([min(t) max(t)])
    xlim([.35 1.25])
    ylim(lims)
    
    % Show legend
    legend('Exact','Weak order $2.0$')
    
    % Label
    xlabel('$x$')
    
    % Ticks
    set(gca,'XTick',0:.2:1.2)
    
