%% Figure 3.2: Two views of Brownian motion
%
% Copyright: 
%   2018 - Simo Särkkä and Arno Solin
%
% License:
%   This software is provided under the MIT License. See the accompanying 
%   LICENSE file for details.

%%

  if exist('rng') % Octave doesn't have rng
    rng(10,'twister')
  else
    randn('state',3);
  end
  
  t = linspace(0,10,1000);
  dt = t(2)-t(1);
  x = [0 cumsum(sqrt(dt)*randn(1,numel(t)-1))];
  
  % For Octave
  norminv = @(x,m,s) -sqrt(2)*erfcinv(2*x) .* s + m;
  normpdf = @(x,m,s) 1./(s .* sqrt(2*pi)) .* exp(-0.5*(x-m).^2 ./ s.^2);
  
  % Upper and lower bounds
  ub = norminv(0.975,zeros(size(t)),sqrt(t));
  lb = norminv(0.025,zeros(size(t)),sqrt(t));
  
  figure(1); clf
    plot(t,x,'-k', ...
         t,0*t,':k',...
         t,ub,'--k', ...
         t,lb,'-.k','LineWidth',0.5)
    
  % Legend
  legend('Sample path of $\beta(t)$','Mean', ...
    'Upper 95\% quantile','Lower 95\% quantile')  
    
  xlabel('Time, $t$')
  ylabel('$\beta(t)$')

  
%% Mesh  

  t = linspace(0,10,21); t(1) = .1;
  x = linspace(-10,10,21);
  [T,X] = meshgrid(t,x);
    
  P = normpdf(X(:),0.*T(:),sqrt(T(:)));
  P = reshape(P,size(T));
  
  figure(2); clf
  
    % Show surface
    if exist('surf2patch') % Matlab    
        p=surf2patch(T,X,P);
        h=patch(p,'EdgeColor','k','FaceColor',[.7 .7 .7],'FaceAlpha',.8,'LineWidth',.25);

        xlabel('Time, $t$')
        ylabel('$\beta(t)$')
        zlabel('$p(\beta(t))$') 

        camproj('perspective')
    else % Octave
        surfl(T,X,P);
        colormap gray;
    end
    
    axis vis3d, grid on, box on
    xlim([0 10])
    ylim([-10 10])
    zlim([0 1])
    
    view([25 25])
     