%% Examples 11.5 and 11.9: Parameter estimation in Ornstein-Uhlenbeck
%
% Copyright: 
%   2018 - Simo Särkkä and Arno Solin
%
% License:
%   This software is provided under the MIT License. See the accompanying 
%   LICENSE file for details.

%%
% Simulate data
%
    % Lock random seed
    if exist('rng') % Octave doesn't have rng
      rng(1);
    else
      randn('state',1);
    end
    
    dt = 0.1;
    lambda = 0.5;
    q = 1;
    steps = 100;

    X = zeros(1,steps);
    T = zeros(1,steps);
    
    x = 0;
    t = 0;
    A = exp(-lambda*dt);
    Q = (q - q*exp(-2*dt*lambda))/(2*lambda);
    for k=1:steps
        X(k) = x;
        T(k) = t;
        x = A*x + sqrt(Q) * randn;
        t = t + dt;
    end
    
    figure(1); clf
    plot(T,X);
    
%%
% Closed-form maximum likelihood
%
    lambda_est = -1/dt * log(sum(X(1:end-1) .* X(2:end)) / ...
        sum(X(1:end-1) .* X(1:end-1)))
    q_est = 1/(length(X)-1) * (2*lambda_est / (1 - exp(-2*lambda_est*dt))) * ...
        sum((X(2:end) - exp(-lambda_est*dt) * X(1:end-1)).^2)
    
%%
% Optimized maximum likelihood
%
    a_f = @(lambda,q) exp(-lambda * dt);
    s_f = @(lambda,q) q/(2*lambda) * (1 - exp(-2*lambda*dt));
    
    ell = @(p) 0.5 * sum(log(2*pi*s_f(p(1),p(2))) + 1/s_f(p(1),p(2)) * (X(2:end) - a_f(p(1),p(2)) * X(1:end-1)).^2);
    
    p = fminsearch(ell,[1 1]);
    lambda_est_2 = p(1)
    q_est2 = p(2)

%%
% Optimized maximum likelihood from Euler-Maruyama
%
    a_f = @(lambda,q) (1 - lambda * dt);
    s_f = @(lambda,q) q*dt;
    
    ell_em = @(p) 0.5 * sum(log(2*pi*s_f(p(1),p(2))) + 1/s_f(p(1),p(2)) * (X(2:end) - a_f(p(1),p(2)) * X(1:end-1)).^2);
    
    p = fminsearch(ell_em,[1 1]);
    lambda_est_em = p(1)
    q_est_em = p(2)
    
%%
% Do some Metropolis-Hastings on the exact model.
%
    % Note that we do a change of variables such that
    % theta = (log q, log lambda), which
    % which leads to a correction term 1/q/lambda
    
    % Lock random seed
    if exist('rng') % Octave doesn't have rng
      rng(1);
    else
      randn('state',2);
    end
    
    nmcmc = 10000;
    samp = [];

    theta = [0 0];

    en = 0;
    accepted = 0;
    
    for j=1:nmcmc

        % Draw candidate
        new_theta = theta + 0.5*randn(1,2);
        
        % Evaluate energy
        new_en = ell(exp(new_theta)) - sum(new_theta);
        
        % Accept or reject
        if j == 1
            en = new_en;
        end
        a = min(1,exp(en - new_en));
        u = rand;

        if u <= a
            en = new_en;
            theta = new_theta;
            accepted = accepted + 1;
            %fprintf('Accepted: %f %f\n',theta);
        else
            %fprintf('Rejected\n');
        end
        samp = [samp; theta];

        if rem(j,1000) == 0
            fprintf('Acceptance rate %f\n',accepted/j);
        end
    end

    figure(1); clf;
    subplot(1,2,1);
    
      plot(exp(samp(:,1)),exp(samp(:,2)),'.');
      xlabel('Parameter $\lambda$');
      ylabel('Parameter $q_c$')

    subplot(1,2,2);
      [N,XE,YE] = ndhist(exp(samp(:,1:2)),20);
      pcolor(XE,YE,N);
      colormap(flipud(gray))
      shading flat;
      xlabel('Parameter $\lambda$');
      ylabel('Parameter $q_c$')

      
%%
% Do some Metropolis-Hastings on the E-M approximation.
%
    % Lock random seed
    if exist('rng') % Octave doesn't have rng
      rng(1);
    else
      randn('state',2);
    end
    
    nmcmc = 10000;
    samp_em = [];

    theta = [0 0];

    en = 0;
    accepted = 0;
    
    for j=1:nmcmc

        % Draw candidate
        new_theta = theta + 0.5*randn(1,2);
        
        % Evaluate energy
        new_en = ell_em(exp(new_theta)) - sum(new_theta);
        
        % Accept or reject
        if j == 1
            en = new_en;
        end
        a = min(1,exp(en - new_en));
        u = rand;

        if u <= a
            en = new_en;
            theta = new_theta;
            accepted = accepted + 1;
            %fprintf('Accepted: %f %f\n',theta);
        else
            %fprintf('Rejected\n');
        end
        samp_em = [samp_em; theta];

        if rem(j,1000) == 0
            fprintf('Acceptance rate %f\n',accepted/j);
        end
    end
    
    figure(1); clf;
    subplot(1,2,1);
      plot(exp(samp_em(:,1)),exp(samp_em(:,2)),'.');
      xlabel('Parameter $\lambda$');
      ylabel('Parameter $q_c$')

    subplot(1,2,2);
      [N,XE,YE] = ndhist(exp(samp_em(:,1:2)),20);
      pcolor(XE,YE,N);
      colormap(flipud(gray))
      shading flat;
      xlabel('Parameter $\lambda$');
      ylabel('Parameter $q_c$')


%%
% Fokker-Planck approximation
%
    L  = 15;
    h = 0.05;
    
    x_grid = (-L:h:L)';
    n = length(x_grid);

    f = -lambda * x_grid;
    
    Fa = q/2/h^2+f/h;
    Fb = -q/h^2-f/h;  
    Fc = repmat(q/2/h^2,n,1);
    
    F = spdiags([Fa Fb Fc],-1:1,n,n);

    fprintf('Computing the matrix exponential...');
    
    T_pde = expm(full(F*dt));
    pcolor(T_pde)
    shading flat

    fprintf('done.\n');
    
    t=0;
    pp = zeros(size(x_grid));
    ind = find(x_grid == 0);
    ind = find(x_grid == 0);
    pp(ind) = 1/h;
    for k=1:steps
         pp = T_pde * pp;
         pp_norm = pp ./ sum(pp) / h;
         t = t + dt;         
         %plot(x_grid,pp_norm)
         %pause(0.1);
         %drawnow;
    end
    
    figure(1); clf;
    plot(x_grid,pp_norm)
    
  %%
  % Plot the whole likelihood surfaces and MCMC results
  %
    [ll,qq] = meshgrid(0.01:0.01:2,0.01:0.01:2);
    lh    = zeros(size(ll));
    lh_em = zeros(size(ll));
    for i=1:size(ll,1)
        for j=1:size(ll,2)
            lh(i,j)    = ell([ll(i,j) qq(i,j)]);
            lh_em(i,j) = ell_em([ll(i,j) qq(i,j)]);
        end
    end
    
    figure(1); clf
    subplot(2,2,1); 
    imagesc(ll(1,:),qq(:,1),exp(-lh));
    hold on;
    
    plot(lambda_est,q_est,'*')
    xlabel('$\lambda$');
    ylabel('$q$');
    title('Exact likelihood');
    ax = axis; axis xy, box on
    colormap(flipud(gray))
    set(gca,'Layer','Top')
    set(gca,'XTick',0:.5:2.5,'YTick',0:.5:2.5)
    
    subplot(2,2,2); 
    n = size(samp,1);
    for i=1:1000:n
       ind = i:min(i+999,n);
       plot(exp(samp(ind,1)),exp(samp(ind,2)),'.','MarkerSize',3,'Color',[.5 .5 .5]);
       hold on;
    end
    xlabel('$\lambda$');
    ylabel('$q$');
    title('Exact MCMC');
    grid on; box on
    axis(ax);
    set(gca,'XTick',0:.5:2.5,'YTick',0:.5:2.5)

    subplot(2,2,3); 
    imagesc(ll(1,:),qq(:,1),exp(-lh_em));
    hold on;
    
    plot(lambda_est_em,q_est_em,'*')
    xlabel('$\lambda$');
    ylabel('$q$');
    title('Pseudo-likelihood');
    ax = axis; axis xy, box on
    colormap(flipud(gray))
    set(gca,'Layer','Top')
    set(gca,'XTick',0:.5:2.5,'YTick',0:.5:2.5)
        
    subplot(2,2,4); 
    n = size(samp_em,1);
    for i=1:1000:n
      ind = i:min(i+999,n);
      plot(exp(samp_em(ind,1)),exp(samp_em(ind,2)),'.','MarkerSize',3,'Color',[.5 .5 .5]);
      hold on;
    end
    xlabel('$\lambda$');
    ylabel('$q$');
    title('Pseudo-MCMC');
    grid on; box on
    axis(ax);
    set(gca,'XTick',0:.5:2.5,'YTick',0:.5:2.5)
   
%%
% Plot final figures
%

  % Exact likelihood
  figure(1); clf; hold on;
    
    imagesc(ll(1,:),qq(:,1),exp(-lh));
    
    h = plot(lambda_est,q_est,'*');
    set(h,'Color',[1 1 1]);
    xlabel('$\lambda$');
    ylabel('$q$');
    ax = axis; axis xy, box on
    colormap(flipud(gray))
    set(gca,'Layer','Top')
    set(gca,'XTick',0:.5:2.5,'YTick',0:.5:2.5)
    
  % Exact MCMC
  figure(2); clf;
    n = size(samp,1);
    for i=1:1000:n
       ind = i:min(i+999,n);
       plot(exp(samp(ind,1)),exp(samp(ind,2)),'.','MarkerSize',3,'Color',[.5 .5 .5]);
       hold on;
    end
    xlabel('$\lambda$');
    ylabel('$q$');
    grid on; box on
    axis(ax);
    set(gca,'XTick',0:.5:2.5,'YTick',0:.5:2.5)
    
  % Pseudo-likelihood
  figure(3); clf;
    imagesc(ll(1,:),qq(:,1),exp(-lh_em));
    hold on;
    h = plot(lambda_est_em,q_est_em,'*');
    set(h,'Color',[1 1 1]);
    xlabel('$\lambda$');
    ylabel('$q$');
    ax = axis; axis xy, box on
    colormap(flipud(gray))
    set(gca,'Layer','Top')
    set(gca,'XTick',0:.5:2.5,'YTick',0:.5:2.5)

  % Pseudo-MCMC
  figure(4); clf;
    n = size(samp_em,1);
    for i=1:1000:n
      ind = i:min(i+999,n);
      plot(exp(samp_em(ind,1)),exp(samp_em(ind,2)),'.','MarkerSize',3,'Color',[.5 .5 .5]);
      hold on;
    end
    xlabel('$\lambda$');
    ylabel('$q$');
    grid on; box on
    axis(ax);
    set(gca,'XTick',0:.5:2.5,'YTick',0:.5:2.5)
    
  % Trajectory
  figure(5); clf;
    h = plot(T,X);
    set(h,'Color',[0 0 0]);
    set(h,'LineWidth',1);
    xlabel('Time, $t$');
    ylabel('$x(t)$');
    