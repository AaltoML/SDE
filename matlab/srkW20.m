function [x,tspan] = srkW20(f,L,tspan,x0,Q,gaussian)
%% SRKW20 - Numerical SDE solver: Stochastic RK, weak 2.0
%
% Syntax:
%   [x,tspan] = srkW20(f,L,tspan,x0,Q,gaussian)
%
% In:
%   f        - Drift function, f(x,t)
%   L        - Diffusion function, L(x,t)
%   tspan    - Time steps to simulate, [t0,...,tend]
%   x0       - Initial condition
%   Q        - Spectral density (default: standard Brownian motion)
%   gaussian - Use Gaussian increments (default: false)
%
% Out:
%   x      - Solved values
%   tspan  - Time steps
%   
% Description:
%   Integrates the system of stochatic differential equations
%     dx = f(x,t) dt + L(x,t) dbeta,  for x(0) = x0
%   over the time interval defined in tspan.
%
% Copyright: 
%   2018 - Simo Särkkä and Arno Solin
%
% License:
%   This software is provided under the MIT License. See the accompanying 
%   LICENSE file for details.

%%

  % Number of steps
  steps = numel(tspan);
  
  % Allocate space
  x = zeros(size(x0,1),steps);

  % Initial state
  x(:,1) = x0;
  
  % Dimensions
  n = size(x0,1);
  m = size(L(x0,tspan(1)),2);
  
  % CHeck inputs
  if nargin < 5 || isempty(Q)
    Q = eye(m);  
  end
  cQ = chol(Q,'lower');
  
  % Pre-calculate random numbers
  if nargin < 6 || isempty(gaussian) || ~gaussian
    foo = rand(m,steps);
    B = sqrt(3)*(foo < 1/6) - sqrt(3)*(foo > 5/6);
    foo = rand(m,steps);
    Z = -1 + 2*(foo > 1/2);
  else
    B = randn(m,steps);
    foo = rand(m,steps);
    Z = -1 + 2*(foo > 1/2);
  end
    
  % Iterate
  for k=2:steps

    % Time discretization
    dt = tspan(k)-tspan(k-1);

    % Increments
    db  = cQ*sqrt(dt)*B(:,k);
    dbb = (db*db')/2 + ...
          sqrt(dt)/2*tril(sqrt(dt)*ones(m,1)*Z(:,k)',-1) - ...
          sqrt(dt)/2*triu(sqrt(dt)*Z(:,k)*ones(1,m),1) - ...
          dt*eye(m)/2;
        
    % Evaluate only once
    fx = f(x(:,k-1),tspan(k-1));
    Lx = L(x(:,k-1),tspan(k-1));
    
    % Supporting values
    xx  = x(:,k-1) + fx*dt;
    tx0 = xx + Lx*db;
    Lxdbb = Lx*dbb;
    tx2 = bsxfun(@plus,xx,Lxdbb);
    tx3 = bsxfun(@plus,xx,-Lxdbb);
    
    % Iterate to obtain the bar x values
    bx2 = tx2; bx3 = tx3;
    for n=1:m
      ind = setdiff(1:m,n);
      bx2(:,n) = xx + sum(Lx(:,ind)*dbb(ind,n),2)/sqrt(dt); %sum(Lxdbb(:,ind),2)/sqrt(dt);
      bx3(:,n) = xx - sum(Lx(:,ind)*dbb(ind,n),2)/sqrt(dt); %sum(Lxdbb(:,ind),2)/sqrt(dt);
    end
    
    % Evaluate
    Ltx2 = L(tx2,tspan(k));
    Ltx3 = L(tx3,tspan(k));
    Lbx2 = L(bx2,tspan(k));
    Lbx3 = L(bx3,tspan(k));
        
    % Make the step
    x(:,k) = x(:,k-1) + ...
      dt/2*(fx + f(tx0,tspan(k))) + ...
      (  Lx/2 + Ltx2/4 + Ltx3/4)*db + ...
      (Ltx2/2 - Ltx3/2         )*diag(dbb)/sqrt(dt) + ...
      ( -Lx/2 + Lbx2/4 + Lbx3/4)*db + ...
      (Lbx2/2 - Lbx3/2         )*ones(m,1)*sqrt(dt);
      
  end
  
