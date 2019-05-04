function x = rejectionsample(q1,epsilon)
%% REJECTIONSAMPLE - Rejection sampling with Gaussian proposals
%
% Syntax:
%   x = rejectionsample(q1,epsilon)
%
% In:
%   q1      - Function handle
%   epsilon - Bounding parameter
%
% Out:
%   x       - Sample
%   
% Description:
%   Rejection sampling with Gaussian proposals.
%
% Copyright: 
%   2018 - Simo Särkkä and Arno Solin
%
% License:
%   This software is provided under the MIT License. See the accompanying 
%   LICENSE file for details.

%% Set defaults

  % Scale constant
  if nargin<2 || isempty(epsilon)
    epsilon = 1;  
  end

  % Proposal distribution (standard Normal)
  % q2 = @(x) normpdf(x,0,1);
  q2 = @(x) exp(-0.5 * x.^2) / sqrt(2*pi);
  
  
%% Draw a sample

  while true
  
    % Sample x from q2
    x = randn(1);
  
    % Sample a uniform random variable
    u = rand(1);
  
    % Accept or reject x
    if u < epsilon*q1(x)/q2(x)
      return
    end
    
  end



