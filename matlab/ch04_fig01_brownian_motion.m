%% Figure 4.1: Brownian motion
%
% Copyright: 
%   2018 - Simo Särkkä and Arno Solin
%
% License:
%   This software is provided under the MIT License. See the accompanying 
%   LICENSE file for details.

%% A sketch of Brownian motion

  if exist('rng') % Octave doesn't have rng
      rng(4,'twister');
  else
      randn('state',2);
  end

  figure(1); clf
    T = 0:0.001:1;
    X = randn(size(T));
    plot([0 T],[0 cumsum(X)*sqrt((T(2)-T(1)))]);
    xlabel('Time, $t$'); ylabel('$\beta(t)$')
