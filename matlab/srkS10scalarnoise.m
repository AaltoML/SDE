function [x,tspan] = srkS10scalarnoise(f,L,tspan,x0,Q)
%% SRKS10SCALARNOISE - Numerical SDE solver: Stochastic RK, strong 1.0
%
% Syntax:
%   [x,tspan] = srkS10scalarnoise(f,L,tspan,x0,Q)
%
% In:
%   f      - Drift function, f(x,t)
%   L      - Diffusion function, L(x,t)
%   tspan  - Time steps to simulate, [t0,...,tend]
%   x0     - Initial condition
%   Q      - Spectral density (default: standard Brownian motion)
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

  % Check if Q given
  if nargin<5 || isempty(Q), Q = eye(size(L(x0,tspan(1)),2)); end 

  % NB: Only for scalar beta
  if size(Q,1)>1, error('NB: Only for scalar beta.'), end

  % Number of steps
  steps = numel(tspan);
  
  % Allocate space
  x = zeros(size(x0,1),steps);

  % Initial state
  x(:,1) = x0;
  
  % Pre-calculate random numbers
  R = randn(1,steps);
  
  % Iterate
  for k=2:steps

    % Time discretization
    dt = tspan(k)-tspan(k-1);

    % Increment
    db  = sqrt(dt*Q)*R(1,k);
    dbb = 1/2*(db^2 - Q*dt);
        
    % Evaluate only once
    fx = f(x(:,k-1),tspan(k-1));
    Lx = L(x(:,k-1),tspan(k-1));
    
    % Supporting values
    x2  = x(:,k-1) + fx*dt;
    tx2 = x2 + Lx*dbb/sqrt(dt);
    tx3 = x2 - Lx*dbb/sqrt(dt);
    
    % Evaluate the remaining values
    fx2 = f(x2,tspan(k-1)+dt);
    Lx2 = L(tx2,tspan(k-1)+dt);
    Lx3 = L(tx3,tspan(k-1)+dt);
    
    % Step
    x(:,k) = x(:,k-1) + ...
        (fx+fx2)*dt/2 + ...
        Lx*db + ...
        sqrt(dt)/2*(Lx2 - Lx3);
    
  end

