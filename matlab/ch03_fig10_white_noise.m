%% Figure 3.10: White noise
%
% Copyright: 
%   2018 - Simo Särkkä and Arno Solin
%
% License:
%   This software is provided under the MIT License. See the accompanying 
%   LICENSE file for details.

%% A sketch of white noise

  if exist('rng') % Octave doesn't have rng
    rng(0,'twister');
  else
    randn('state',0);
  end

  figure(1); clf
    T = 0:0.001:1;
    X = randn(size(T));
    plot(T,X);
    xlabel('Time, $t$'); ylabel('$w(t)$')

    
