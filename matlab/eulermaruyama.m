function [x,tspan] = eulermaruyama(f,L,tspan,x0,Q)
%% EULERMARUYAMA - Numerical SDE solver: The Euler-Maruyama method
%
% Syntax:
%   [x,tspan] = eulermaruyama(f,L,tspan,x0,Q)
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

  % Cholesky factor of Q
  cQ = chol(Q,'lower');
  
  % Number of steps
  steps = numel(tspan);
  
  % Allocate space
  x = zeros(size(x0,1),steps);

  % Initial state
  x(:,1) = x0;
  
  % Iterate
  for k=2:steps

    % Time discretization
    dt = tspan(k)-tspan(k-1);

    % Increment
    db = sqrt(dt)*cQ*randn(size(Q,1),1);

    % Step
    x(:,k) = x(:,k-1) + ...
        f(x(:,k-1),tspan(k-1))*dt + ...
        L(x(:,k-1),tspan(k-1))*db;
    
  end


