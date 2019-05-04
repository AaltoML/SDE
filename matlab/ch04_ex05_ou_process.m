%% Example 4.5: Solution of the Ornstein–Uhlenbeck process
%
% Copyright: 
%   2018 - Simo Särkkä and Arno Solin
%
% License:
%   This software is provided under the MIT License. See the accompanying 
%   LICENSE file for details.

%% Simulate trajectories from the OU process

  if exist('rng') % Octave doesn't have rng
      rng(10,'twister')
  else
      randn('state',1)
  end
  
  lambda = 0.5;
  q = 1;
  dt = 0.01;
  T = (0:dt:1);
  x0 = 4;
  M = exp(-lambda*T)*x0;
  P = q/(2*lambda)*(1 - exp(-2*lambda*T));

  XX = zeros(50,length(T));
  for n=1:size(XX,1)
    x = x0;
    for k=1:length(T)
      XX(n,k) = x;
      x = x - lambda * x * dt + sqrt(dt)*randn;
    end
  end
  
  
  figure(1); clf; hold on
  
    % Shade the 95% quantiles
    fill([T fliplr(T)],[M+1.96*sqrt(P) fliplr(M-1.96*sqrt(P))],1, ...
      'FaceColor',[.9 .9 .9],'EdgeColor',[.9 .9 .9])
  
    % Plot realizations
    h1 = plot(T,XX,'-','Color',[.5 .5 .5],'LineWidth',0.5);
    
    % Plot mean and quantiles
    h2 = plot(T,M,'k-','LineWidth',1);
    h34 = plot(T,M+1.96*sqrt(P),'--k', ...
               T,M-1.96*sqrt(P),'--k','LineWidth',0.7);
    
    %set(h(1:3),'Linewidth',2)
    
    legend([h2(1) h34(1) h1(1)],'Mean','95\% quantiles','Realizations')
    xlabel('Time, $t$'); ylabel('$x(t)$')
    
    ylim([0 5])
    box on
    set(gca,'Layer','Top')
    