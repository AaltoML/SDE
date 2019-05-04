%% Examples 12.6 and 12.11: Batch and sequential solution to GP regression
%
% Copyright: 
%   2018 - Simo Särkkä and Arno Solin
%
% License:
%   This software is provided under the MIT License. See the accompanying 
%   LICENSE file for details.

%% Batch GP regression

  % Lock random seed
  if exist('rng') % Octave doesn't have rng
    rng(0,'twister')
  else
    randn('state',0);
  end

  % Define parameters
  sigma2 = 0.1^2;
  magnSigma2 = 1^2;
  ell = .1;
  
  % Define covariance function: Matern 3/2
  C = @(x,y) magnSigma2*(1+sqrt(3)*abs(x(:)-y(:)')/ell).* ...
      exp(-sqrt(3)*abs(x(:)-y(:)')/ell);
  
  % Number of inputs
  n = 10;
  to = sort(rand(n,1));
  
  % Simulate data
  Koo = C(to,to);
  cK = chol(Koo);
  yo = cK'*randn(size(cK,1),1) + sqrt(sigma2)*randn(n,1);
  
  % Evaluate test points
  t = linspace(0,1,500)';
  Kto = C(t,to);
  Ktt = C(t,t);
  
  % Solve the GP regression problem
  L = chol(Koo+sigma2*eye(n),'lower');
  m_gp = Kto*((Koo+sigma2*eye(n))\yo);
  cov_gp = Ktt - Kto*((Koo+sigma2*eye(n))\Kto');
  var_gp = diag(cov_gp);
  
  % Samples form the posterior GP
  ns = 5;
  cK = chol(cov_gp);
  xs = m_gp + cK'*randn(size(cK,1),ns);
  
  % Plot
  figure(1); clf; hold on
  
    % Plot mean
    h1 = plot(t,m_gp,'-k');
    
    % Plot uncetainty
    h2 = fill([t; flipud(t)],[m_gp+1.96*sqrt(var_gp); flipud(m_gp-1.96*sqrt(var_gp))],1);
    set(h2,'EdgeColor',[.7 .7 .7],'FaceColor',[.7 .7 .7])
    
    % Plot sample trajectories
    h4 = plot(t,xs,'-','Color',[.5 .5 .5],'LineWidth',.5);
    
    % Plot mean
    plot(t,m_gp,'-k');
    
    % Plot observations
    h3 = plot(to,yo,'+k');
    
    % Axis options
    set(gca,'layer','top')
    box on
    ylim([-3 4])
    xlabel('Input, $\xi$'), ylabel('Output, $x(\xi)$')
    legend([h1, h2, h3, h4(1)],'Mean','95\% quantiles','Observations','Samples')
    
    
%% GP regression by Kalman filtering and RTS smoothing

  % Define the LTI SDE model: Matern 3/2
  lambda = sqrt(3)/ell;
  F = [0 1; -lambda^2 -2*lambda];
  L = [0; 1];
  Q = 4*lambda^3*magnSigma2;
  H = [1 0];
  Pinf = magnSigma2*diag([1, lambda^2]);
  
  % Combine test and observation times
  [tt,ind] = sort([to; t]);
  
  % Time steps with observations
  doUpdate = (ind<=n);
    
  % Set initial mean and covariance
  m = zeros(size(F,1),1);
  P = Pinf;
  
  % Allocate space for results
  MF = zeros(size(m,1),numel(tt));
  PF = zeros(size(m,1),size(m,1),numel(tt));
  MP = zeros(size(m,1),numel(tt));
  PP = zeros(size(m,1),size(m,1),numel(tt));
  GS = zeros(size(m,1),size(m,1),numel(tt));
  
  % Run the Kalman filter
  for k=1:numel(tt)
      
    % Time step
    if k==1
      dt = 0;
    else
      dt = tt(k)-tt(k-1);  
    end
      
    % Solve the LTI SDE for this step
    A = expm(F*dt);
    Q = Pinf - A*Pinf*A';
    
    % Kalman prediction
    mp = A*m;
    Pp = A*P*A' + Q;
    
    % Pre-calculate smoother gain
    Gs = (P*A')/Pp;

    % Do Kalman update if there is an observation on this step
    if doUpdate(k)
      v = yo(ind(k)) - H*m;
      S = H*P*H' + sigma2;
      K = P*H'/S;
      m = m + K*v;
      P = P - K*H*P;
    else
      m = mp;
      P = Pp;
    end
    
    % Store results
    MF(:,k) = m;
    PF(:,:,k) = P;
    MP(:,k) = mp;
    PP(:,:,k) = Pp;
    GS(:,:,k) = Gs;
    
  end
  
  % Run the RTS smoother
  MS = MF;
  PS = PF;
  ms = MS(:,end);
  Ps = PS(:,:,end);
  for j=numel(tt)-1:-1:1
    ms = MF(:,j) + GS(:,:,j+1)*(ms-MP(:,j+1));
    Ps = PF(:,:,j) + GS(:,:,j+1)*(Ps - PP(:,:,j+1))*GS(:,:,j+1)';
    MS(:,j) = ms;
    PS(:,:,j) = Ps;    
  end
    
  % Extract values
  m_f = (H*MF)';
  var_f = arrayfun(@(k) H*PF(:,:,k)*H',1:numel(tt))';
  m_s = (H*MS)';
  var_s = arrayfun(@(k) H*PS(:,:,k)*H',1:numel(tt))';
  
  % Only the test inputs
  ii = ind>n;
  m_f   = m_f(ii);
  var_f = var_f(ii);
  m_s   = m_s(ii);
  var_s = var_s(ii);
  

  % Plot
  figure(2); clf; hold on
  
    % Plot uncetainty
    h2 = fill([t; flipud(t)], ...
              [m_f+1.96*sqrt(var_f); flipud(m_f-1.96*sqrt(var_f))],1);
    set(h2,'EdgeColor',[.7 .7 .7],'FaceColor',[.7 .7 .7])
    
    plot(t,m_gp+1.96*sqrt(var_gp),'--k', ...
         t,m_gp-1.96*sqrt(var_gp),'--k');
        
    % Plot mean
    plot(t,m_f,'-k')
    
    % Plot observations
    plot(to,yo,'+k')
    
    % Axis options
    set(gca,'layer','top')
    box on    
    xlabel('Input, $t$'), ylabel('Output, $x(t)$')
    ylim([-3 4])
    
  % Plot
  figure(3); clf; hold on
  
    % Plot mean
    h1 = plot(t,m_s,'-k');
    
    % Plot uncetainty
    h2 = fill([t; flipud(t)], ...
              [m_s+1.96*sqrt(abs(var_s)); flipud(m_s-1.96*sqrt(abs(var_s)))],1);
    set(h2,'EdgeColor',[.7 .7 .7],'FaceColor',[.7 .7 .7])
    
    % Plot mean
    plot(t,m_s,'-k')
    
    % Plot observations
    h3 = plot(to,yo,'+k');
 
    h4 = plot(t,m_gp+1.96*sqrt(var_gp),'--k', ...
              t,m_gp-1.96*sqrt(var_gp),'--k');
    
    % Axis options
    set(gca,'layer','top')
    box on    
    xlabel('Input, $t$'), ylabel('Output, $x(t)$')
    ylim([-3 4])
    legend([h1, h2, h3, h4(1)],'Mean','95\% quantiles','Observations','Batch solution')
  
  
  