%% Example 8.11: Comparison of ODE solvers
%
% Copyright: 
%   2018 - Simo Särkkä and Arno Solin
%
% License:
%   This software is provided under the MIT License. See the accompanying 
%   LICENSE file for details.

%% Plot phase portrait

  % Define the system in terms of the differential equation
  f = @(x) [x(1,:)-x(2,:)-x(1,:).^3; 
            x(1,:)+x(2,:)-x(2,:).^3];
  
  % Define arrow
  arrow1 = [-1 1 0 -1; -.5 -.5 2 -.5]';
  
  % Define arrow without head
  arrow2 = .4*[-1 -1 1 1 -1;-.1 .1 .1 -.1 -.1]';
  
  % Phase portrait figure
  figure(1); clf; hold on
  
    % Discrete grid
    x = linspace(-2,2,15);
    y = linspace(-2,2,15);
    [X Y] = meshgrid(x,y);
    
    % Form vector field
    UV = f([X(:) Y(:)]');
    U = reshape(UV(1,:),size(X));
    V = reshape(UV(2,:),size(Y));
    
    % Plot field
    X(X==0 & Y==0) = nan;
    newquiver(X,Y,U,V,'X',arrow2, ...
        'color',[.7 .7 .7],'EdgeColor',[.7 .7 .7])
    
    % Change axis
    axis equal; axis([-2.2 2.2 -2.2 2.2]);
    set(gca,'XTick',[-2:4:2],'YTick',[-2:4:2])
        
    % Show some trajectories (start from outside)
    x0 = 3*[cos(linspace(0,2*pi,25)); sin(linspace(0,2*pi,25))];
    
    tin = 0:2^-4:10;
    for i=1:6
      try
        xout = euler(@(x,t) f(x),tin,x0(:,i)); xout = xout';
        plot(xout(:,1),xout(:,2),'-k','LineWidth',.5,'Color',[.5 .5 .5])
        uv = f(xout'); ind = 4;
        newquiver(xout(ind,1),xout(ind,2),uv(2,ind),-uv(1,ind), ...
          'X',arrow1,'scale',.04)
      end
    end
    for i=7:12
      try
        xout = heun(@(x,t) f(x),tin,x0(:,i)); xout = xout';
        plot(xout(:,1),xout(:,2),'-k','LineWidth',.5)
        uv = f(xout'); ind = 4;
        newquiver(xout(ind,1),xout(ind,2),uv(2,ind),-uv(1,ind), ...
          'X',arrow1,'scale',.04)
      end
    end
    for i=13:18
      try
        xout = impliciteuler(@(x,t) f(x),tin,x0(:,i)); xout = xout';
        plot(xout(:,1),xout(:,2),'-k','LineWidth',.5,'Color',[.5 .5 .5])
        uv = f(xout'); ind = 4;
        newquiver(xout(ind,1),xout(ind,2),uv(2,ind),-uv(1,ind), ...
          'X',arrow1,'scale',.04)
      end      
    end
    for i=19:24
      try
        xout = rk4simple(@(x,t) f(x),tin,x0(:,i)); xout = xout';
        plot(xout(:,1),xout(:,2),'-k','LineWidth',.5)
        uv = f(xout'); ind = 4;
        newquiver(xout(ind,1),xout(ind,2),uv(2,ind),-uv(1,ind), ...
          'X',arrow1,'scale',.04)
      end      
    end
    
    % Labels
    text( 1.9, 1.9,'Forward Euler','HorizontalAlignment','right')
    text(-1.9, 1.9,'Heun','HorizontalAlignment','left')
    text(-1.9,-1.9,'Backward Euler','HorizontalAlignment','left')
    text( 1.9,-1.9,'RK4','HorizontalAlignment','right')
    
    % Show some trajectories (start from inside)
    x0 = .02*[cos(linspace(0,2*pi,5)); sin(linspace(0,2*pi,5))];
    for i=1:size(x0,2)
      try
        xout = rk4simple(@(x,t) f(x),linspace(0,10,48),x0(:,i)); xout = xout';
        plot(xout(:,1),xout(:,2),'-k','LineWidth',.5)
        uv = f(xout'); ind = 18;
        newquiver(xout(ind,1),xout(ind,2),uv(2,ind),-uv(1,ind), ...
          'X',arrow1,'scale',.04)
      end
    end
    
    % Axis labels
    xlabel('$x_1$')
    ylabel('$x_2$')
    
     % Plot fixed point
    plot(0,0,'ok','MarkerFaceColor','w','MarkerSize',3)
    
  % Set figure options
  set(gcf,'Color','w')
  
 