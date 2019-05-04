%% Example 
%
% Copyright: 
%   2018 - Simo Särkkä and Arno Solin
%
% License:
%   This software is provided under the MIT License. See the accompanying 
%   LICENSE file for details.

%% Leapfrog / Verlet integration

dts = logspace(-2,0,10);
err = nan(10,numel(dts));

for j=1:numel(dts)
    fprintf('Running step %d/%d\n',j,numel(dts));

% Leapfrog / Verlet integration

  % Lock seed
  if exist('rng') % Octave doesn't have rng
    rng(j,'twister')
  else
    randn('state',j);
  end

  % Parameters
  n = 25000;
  g = 1; % g>0
  eta = 1;
  q = 1;
  x0 = 1;
  v0 = 0;

  % Time discretization
  dt = dts(j); %0.1;
  t = 0:dt:10;
  
  % Specify the model
  f = @(x) -g*x;
  s = @(x) 1;
   
  % Allocate space and set initial conditions
  x = zeros(n,numel(t)); x(:,1) = x0;
  v = zeros(n,numel(t)); v(:,1) = v0;
  
  % For each trajectory
  for i=1:n

    % Simulate one trajectory from the leapfrog method
    z = leapfrog(f,s,eta,t,[x0; v0],q);
    
    x(i,:) = z(1,:);
    v(i,:) = z(2,:);
    
  end
  
  figure(1); clf
  subplot(211)
    plot(t,x(1:5,:), ...
         t,mean(x), ...
         t,mean(x)+1.96*std(x),'--', ...
         t,mean(x)-1.96*std(x),'--')
    xlabel('Time, t'), ylabel('Position, x(t)')     
  subplot(212)
    plot(t,v(1:5,:), ...
         t,mean(v), ...
         t,mean(v)+1.96*std(v),'--', ...
         t,mean(v)-1.96*std(v),'--')
    xlabel('Time, t'), ylabel('Velocity, v(t)')     
  
    
% Euler-Maruyama

  % Same random seed
  % Lock seed
  if exist('rng') % Octave doesn't have rng
    rng(j,'twister')
  else
    randn('state',j);
  end

  % Specify the model
  f = @(x,t) [0 1; -g -eta]*x;
  L = @(x,t) [0; 1];
  
  % Allocate space and set initial conditions
  x_em = zeros(n,numel(t)); x_em(:,1) = x0;
  v_em = zeros(n,numel(t)); v_em(:,1) = v0;
  
  % For each trajectory
  for i=1:n

    % Simulate one trajectory from the leapfrog method
    z = eulermaruyama(f,L,t,[x0; v0],q);
    
    x_em(i,:) = z(1,:);
    v_em(i,:) = z(2,:);
    
  end
  
  figure(1); clf
  subplot(211)
    plot(t,x_em(1:5,:), ...
         t,mean(x_em), ...
         t,mean(x_em)+1.96*std(x_em),'--', ...
         t,mean(x_em)-1.96*std(x_em),'--')
    xlabel('Time, t'), ylabel('Position, x(t)')
  subplot(212)
    plot(t,v(1:5,:), ...
         t,mean(v_em), ...
         t,mean(v_em)+1.96*std(v_em),'--', ...
         t,mean(v_em)-1.96*std(v_em),'--')
    xlabel('Time, t'), ylabel('Velocity, v(t)')
     
  
% Closed-form solution

  % This is a LTI SDE of the form
  F = [0 1; -g -eta];
  L = [0; 1];
  Qc = q;

  % Discretize
  [A,Q] = lti_disc(F,L,Qc,dt);
  
  % Evaluate
  M = zeros(2,numel(t)); M(:,1) = [x0; v0];
  P = zeros(2,2,numel(t));
  
  % Loop through
  for k=1:numel(t)-1
    M(:,k+1) = A*M(:,k);
    P(:,:,k+1) = A*P(:,:,k)*A'+Q;
  end
  
  figure(2); clf
  subplot(211)
    plot(t,M(1,:), ...
         t,M(1,:)'+1.96*sqrt(squeeze(P(1,1,:))),'--', ...
         t,M(1,:)'-1.96*sqrt(squeeze(P(1,1,:))),'--')
  subplot(212)
    plot(t,M(2,:), ...
         t,M(2,:)'+1.96*sqrt(squeeze(P(2,2,:))),'--', ...
         t,M(2,:)'-1.96*sqrt(squeeze(P(2,2,:))),'--')

     
% Calculate errors
  
  err(1,j) = P(1,1,end)-std(x(:,end))^2;
  err(2,j) = P(2,2,end)-std(v(:,end))^2;
  
  err(3,j) = P(1,1,end)-std(x_em(:,end))^2;
  err(4,j) = P(2,2,end)-std(v_em(:,end))^2;

  err(5,j) = M(1,end)-mean(x(:,end));
  err(6,j) = M(2,end)-mean(v(:,end));
  
  err(7,j) = M(1,end)-mean(x_em(:,end));
  err(8,j) = M(2,end)-mean(v_em(:,end));
  
  err(9,j)  = mean(abs(squeeze(P(1,1,:))'-std(x).^2));
  err(10,j) = mean(abs(squeeze(P(1,1,:))'-std(x_em).^2));
  
  
end


%% Compose results

  figure(1); clf
    loglog(dts,abs(err(9,:)),'-k', ...
           dts,abs(err(10,:)),'--k')
    xlabel('Time step length, $\Delta t$')
    ylabel('Mean absolute error in $\sigma_x^2$')
    legend('Leapfrog Verlet', 'Euler--Maruyama')
    %ylim([0.4 10])
    xlim([min(dts) max(dts)])
    set(gca,'XTick',[0.01 0.1 1],'XTickLabel',[0.01 0.1 1])
    set(gca,'YTick',[1e-3 1e-2 1e-1 1e0 1e1], ...
            'YTickLabel',{'','$10^{-2}$','$10^{-1}$','$10^{0}$','$10^1$'})
