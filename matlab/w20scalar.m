function [x,tspan] = w20scalar(f,L,tspan,x0,Q,gaussian)
%% W20scalar - Weak order 2.0 Ito-Taylor
%
% Syntax:
%   [x,tspan] = w20scalar(f,L,tspan,x0,Q,gaussian)
%
% In:
%   f        - Drift function and derivatives, {f, df, ddf}
%   L        - Diffusion function and derivatives, {L, dL, ddL}
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

  % Check if Q given
  if nargin<5 || isempty(Q), Q = eye(size(L(x0,tspan(1)),2)); end 

  % Check the sort of increments to use
  if nargin<6 || isempty(gaussian), gaussian = false; end
  
  % Extract derivatives
  ddf = f{3};
  df = f{2};
  f = f{1};
  ddL = L{3};
  dL = L{2};
  L = L{1};
  
  % Cholesky factor of Q
  cQ = chol(Q,'lower');
  
  % Number of steps
  steps = numel(tspan);
  
  % Allocate space
  x = zeros(size(x0,1),steps);

  % Initial state
  x(:,1) = x0;
  
  % Pre-calculate random numbers
  RR = rand(size(Q,1),steps);
  R = (RR<1/6)*sqrt(3) - (RR>5/6)*sqrt(3);
  
  % Iterate
  for k=2:steps

    % Time discretization
    dt = tspan(k)-tspan(k-1);

    % Increment
    if gaussian
      db = sqrt(dt)*cQ*randn;
      dz = 1/2*db*dt;    
    else
      db = sqrt(dt)*cQ*R(:,k);
      dz = 1/2*dt^(3/2)*R(:,k);
    end
    
    % Evaluate only once
    fx = f(x(:,k-1),[]);
    dfx = df(x(:,k-1),[]);
    ddfx = ddf(x(:,k-1),[]);
    Lx = L(x(:,k-1),[]);
    dLx = dL(x(:,k-1),[]);
    ddLx = ddL(x(:,k-1),[]);
    
    % Step
    x(:,k) = x(:,k-1) + ...
        fx*dt + ...
        Lx*db + ...
        1/2*Lx*dLx*(db^2 - dt) + ...
        Lx*df(x(:,k-1),[])*dz + ...
        1/2*(fx*dfx + 1/2*Lx^2*ddfx)*dt^2 + ...
        (fx*dLx - 1/2*Lx^2*ddLx)*(db*dt - dz);
    
  end


