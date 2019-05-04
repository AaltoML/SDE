%% Example 9.18: Hermite expansion of Benes SDE
%
% Copyright: 
%   2018 - Simo Särkkä and Arno Solin
%
% License:
%   This software is provided under the MIT License. See the accompanying 
%   LICENSE file for details.

%% Plot the Hermite expansion example

    x0 = 1/2;
    trans_dens = @(xx,t) 1./sqrt(2*pi*t).*cosh(xx)./cosh(x0).*exp(-0.5*t).*exp(-1./(2*t).*(xx-x0).^2);

    %
    % % Symbolic solution (requires the Symbolic Math Toolbox)
    % syms x;
    % f = matlabFunction(tanh(x));
    % df = matlabFunction(diff(tanh(x),1));
    % d2f = matlabFunction(diff(tanh(x),2));
    % d3f = matlabFunction(diff(tanh(x),3));
    % d4f = matlabFunction(diff(tanh(x),4));
    % d5f = matlabFunction(diff(tanh(x),5));
    %
    
    f = @(x) tanh(x);
    df = @(x) 1 - tanh(x)^2;
    d2f = @(x) 2*tanh(x)*(tanh(x)^2 - 1);
    d3f = @(x) -2*(tanh(x)^2 - 1)^2 - 4*tanh(x)^2*(tanh(x)^2 - 1);
    d4f = @(x) 16*tanh(x)*(tanh(x)^2 - 1)^2 + 8*tanh(x)^3*(tanh(x)^2 - 1);
    d5f = @(x) -16*(tanh(x)^2 - 1)^3 - 16*tanh(x)^4*(tanh(x)^2 - 1) - 88*tanh(x)^2*(tanh(x)^2 - 1)^2;
    
    xx = -15:0.1:15;

    figure(1); clf; hold on
    
      % Solve at t=2
      t = 2;
      herm_x = herm_exp(f(x0),df(x0),d2f(x0),d3f(x0),d4f(x0),d5f(x0),t,xx,x0);

      % Exact density
      fill(xx,trans_dens(xx,t),1,'FaceColor',[.7 .7 .7],'EdgeColor',[.7 .7 .7])
      
      % Series approximation
      plot(xx,herm_x,'-k','LineWidth',1)

      ylim([-0.21 0.25]);      
      legend('Exact density','Series approximation','Location','south');
      xlabel('$x$');
      ylabel('$p(x)$');

    figure(2); clf; hold on

      % Solve at t=5    
      t = 5;
      herm_x = herm_exp(f(x0),df(x0),d2f(x0),d3f(x0),d4f(x0),d5f(x0),t,xx,x0);

      % Exact density
      fill(xx,trans_dens(xx,t),1,'FaceColor',[.7 .7 .7],'EdgeColor',[.7 .7 .7])
      
      % Series approximation
      plot(xx,herm_x,'-k','LineWidth',1)
      
      ylim([-0.21 0.25]);
      xlabel('$x$');
      ylabel('$p(x)$');
    