function [A,Q] = lti_disc(F,L,Q,dt)
% LTI_DISC Equivalent discrete-time solution of an LTI SDE
%
% Syntax:
%   [A,Q] = lti_disc(F,L,Qc,dt)
%
% In:
%   F  - NxN Feedback matrix
%   L  - NxL Noise effect matrix        (optional, default identity)
%   Qc - LxL Diagonal Spectral Density  (optional, default zeros)
%   dt - Time step                      (optional, default 1)
%
% Out:
%   A - Transition matrix
%   Q - Discrete process covariance matrix
%
% Description:
%   Equivalent discrete-time solution of an LTI SDE of form
%
%     dx/dt = F x + L w,
%
%   where w(t) is a white-noise process with spectral density Qc.
%   Results in the following model:
%
%     x[k] = A x[k-1] + q, q ~ N(0,Q).
%
%   Can be used for integrating the model exactly over time steps, 
%   which are multiples of dt.
%
% Copyright: 
%   2019 - Simo Särkkä and Arno Solin
%
% License:
%   This software is provided under the MIT License. See the accompanying 
%   LICENSE file for details.

  % Check number of arguments and defaults
  if nargin < 1
    error('Too few arguments');
  end
  if nargin < 2
    L = [];
  end
  if nargin < 3
    Q = [];
  end
  if nargin < 4
    dt = [];
  end
  if isempty(L)
    L = eye(size(F,1));
  end
  if isempty(Q)
    Q = zeros(size(F,1),size(F,1));
  end
  if isempty(dt)
    dt = 1;
  end

  % Closed-form integration of the transition matrix
  A = expm(F*dt);

  % Closed-form integration of the covariance matrix
  n   = size(F,1);
  Phi = [F L*Q*L'; zeros(n,n) -F'];
  AB  = expm(Phi*dt)*[zeros(n,n);eye(n)];
  Q   = AB(1:n,:)*A'; % A' = inv(AB((n+1):(2*n),:));