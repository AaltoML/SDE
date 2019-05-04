function [x,tspan] = rk4simple(f,tspan,x0)
%% RK4SIMPLE - Numerical ODE solver: The fourth-order Runge-Kutta method
%
% Syntax:
%   [x,tspan] = rk4simple(f,tspan,x0)
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
%   2016-2018 - Simo Särkkä and Arno Solin
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

    % Stages
    dx1 = f(x(:,k-1),tspan(k-1))*dt;
    dx2 = f(x(:,k-1)+dx1/2,tspan(k-1)+dt/2)*dt;
    dx3 = f(x(:,k-1)+dx2/2,tspan(k-1)+dt/2)*dt;
    dx4 = f(x(:,k-1)+dx3,tspan(k-1)+dt)*dt;

    % Step
    x(:,k) = x(:,k-1) + 1/6*(dx1 + 2*dx2 + 2*dx3 + dx4);
    
  end
  
  