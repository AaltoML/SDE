%% Example 7.3: Karhunen–Loeve series of Brownian motion
%
% Copyright: 
%   2018 - Simo Särkkä and Arno Solin
%
% License:
%   This software is provided under the MIT License. See the accompanying 
%   LICENSE file for details.

%%
% Series expansion
%
    T = 10;
    NN = [10 100 1000];

    wd = [2.5 1.5 0.5];
    cl = [0.8 0.6 0.0];
    
    legs = {};
    
  figure(3); clf
    for i=1:length(NN)
        N = NN(i);
        if exist('rng') % Octave doesn't have rng
            rng(7,'twister')
        else
            randn('state',2)
        end
        n = 1:N;
        t = 0:0.01:T;
        lam = (2*T./(2*n-1)*pi).^2;
        [tt,nn] = meshgrid(t,n);
        phi = sqrt(2/T)*sin((2*nn-1)*pi.*tt./(2*T));
        z = sqrt(lam) .* randn(size(lam));
        h = plot(t,z*phi);
        set(h,'LineWidth',wd(i));
        set(h,'Color',cl(i)*[1 1 1]);
        hold on;
        
        legs{end+1} = sprintf('$N = %d$', NN(i));
    end
    grid on;
    
    xlabel('Time, $t$')
    ylabel('$\beta(t)$')
  
    legend(legs(:));
    