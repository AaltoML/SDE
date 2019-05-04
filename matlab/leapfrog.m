function [x,tspan] = leapfrog(f,s,eta,tspan,x0,Q)
%% LEAPFROG - Numerical SDE solver: The simple leapfrog Verlet method
%
% Syntax:
%   [x,tspan] = leapfrog(f,L,tspan,x0,Q)
%
% In:
%   f      - Function f(x), see below 
%   s      - Function s(x), see below 
%   nu     - Scalar constant, see below
%   tspan  - Time steps to simulate, [t0,...,tend]
%   x0     - Initial condition
%   Q      - Spectral density (default: standard Brownian motion)
%
% Out:
%   x      - Solved values
%   tspan  - Time steps
%   
% Description:
%   Integrates the second-order stochatic differential equation
%     \ddot{x} = f(x) - \nu s^2(x) \dot{x} + s(x) w(t),  for x(0) = x0
%   over the time interval defined in tspan.
%
% References:
%   [1] Burrage, Lenane, and Lythe (2007). NUMERICAL METHODS FOR 
%       SECOND-ORDER STOCHASTIC DIFFERENTIAL EQUATIONS. SIAM J. 
%       SCI. COMPUT.
%
% Copyright: 
%   2018 - Simo Särkkä and Arno Solin
%
% License:
%   This software is provided under the MIT License. See the accompanying 
%   LICENSE file for details.

%%

  % Check if Q given
  if nargin<6 || isempty(Q), Q = 1; end 

  % Square-root of Q
  cQ = sqrt(Q);
  
  % Number of steps
  steps = numel(tspan);
  
  % Allocate space
  x = zeros(size(x0,1),steps);

  % Initial state
  x(:,1) = x0;
  
  % Iterate
  for k=1:steps-1

    % Time discretization
    dt = tspan(k+1)-tspan(k);
    
    % Gaussian increment
    db = sqrt(dt)*cQ*randn(1);
      
    % Supporting half-step value
    hx = x(1,k) + 1/2*x(2,k)*dt;
      
    % Velocity
    x(2,k+1) = x(2,k) - eta*s(hx)^2*x(2,k)*dt + f(hx)*dt + s(hx)*db;
      
    % Position
    x(1,k+1) = hx + 1/2*x(2,k+1)*dt;
    
  end
  
  