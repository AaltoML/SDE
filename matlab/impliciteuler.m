function [x,tspan] = impliciteuler(f,tspan,x0)
%% IMPLICITEULER - Numerical ODE solver: Implicit Euler method
%
% Syntax:
%   [x,tspan] = euler(f,tspan,x0)
%
% In:
%   f      - Function handle, f(x,t)
%   tspan  - Time steps to simulate, [t0,...,tend]
%   x0     - Initial condition
%
% Out:
%   x      - Solved values
%   tspan  - Time steps
%   
% Description:
%   Integrates the system of differential equations
%     x' = f(x,t),  for x(0) = x0
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
  
  % Iterate
  for k=2:steps

    % Time discretization
    dt = tspan(k)-tspan(k-1);

    % Step solve the algebraic problem
    x(:,k) = fsolve(@(z) x(:,k-1) + f(z,tspan(k))*dt - z, ...
        x(:,k-1), optimset('display','none'));
    
  end
  
  