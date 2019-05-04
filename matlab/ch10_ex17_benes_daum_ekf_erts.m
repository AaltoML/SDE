%% Examples 10.17, 10.26, and 10.38: Benes-Daum and EKF/ERTS examples
%
% Copyright: 
%   2018 - Simo Särkkä and Arno Solin
%
% License:
%   This software is provided under the MIT License. See the accompanying 
%   LICENSE file for details.


%%
% Benes and Benes-Daum filter tests
%
    % Check that we have the correct normalization constants and
    % moments. We should have (for x = x(t) with t >= t_k)
    %
    % p(x | Y_k) = 
    % 1/sqrt{2 pi P} exp(-(x-m)^2/(2P)) exp(-1/2 P) cosh(x) / cosh(m)
    %
    % and
    %
    %    E[x] = m + tanh(m) P
    %  E[x^2] = m^2 + P^2 + 2 m tanh(m) P + P
    %    V[x] = P + [1 - tanh(m)^2] P^2
    
    dx = 0.001;
    xx = -15:dx:15;
    dens = @(x,m,P) 1/sqrt(2*pi*P)*exp(-P/2)/cosh(m) * exp(-(x-m).^2/(2*P)) .* cosh(x);

    m = 1;
    P = 5;
    pp = dens(xx,m,P);
    figure(1); clf; plot(xx,pp);
    
    sum(pp .* dx)
    [sum(xx .* pp * dx)        m+tanh(m)*P]
    [sum(xx.^2 .* pp * dx)  m^2+P^2+2*m*tanh(m)*P+P]
    [sum((xx-m-tanh(m)*P).^2 .* pp * dx) P+(1-tanh(m)^2)*P^2]
    
%%
% Generate data
%
  % Lock random seed
  if exist('rng') % Octave doesn't have rng
    rng(1) 
  else
    randn('state',9);
  end
    
    steps = 500;
    meas_mod = 20;
    dt = 0.01;
    R = 0.1;

    X = [];
    Y = [];
    T = [];
    x0 = 0;
    x = x0;
    t = 0;
    for k=1:steps
        x = x + tanh(x) * dt + sqrt(dt)*randn;
        if rem(k,meas_mod) == 0
            y = x + sqrt(R)*randn;
        else
            y = NaN;
        end
        t = t + dt;
        X(k) = x;
        Y(k) = y;
        T(k) = t;
    end
    
    figure(1); clf
    plot(T,X,T,Y,'o');
    
%%
% Benes-Daum filter
%
    EB = zeros(1,length(Y));
    VB = zeros(1,length(Y));
    MMB = zeros(1,length(Y));
    PPB = zeros(1,length(Y));

    mb0 = 0.1;
    Pb0 = 1;
    
    mb = mb0;
    Pb = Pb0;
    
    m0 = mb + Pb * tanh(mb);
    P0 = Pb + (1-tanh(mb)^2) * Pb^2;

    for k=1:length(Y)
        %
        % The Benes-Daum filter
        %
        Pb = Pb + dt;
        if ~isnan(Y(k))
            mb = mb + Pb / (R + Pb)*(Y(k) - mb);
            Pb = Pb - Pb^2/(R + Pb);
        end
        MMB(:,k) = mb;
        PPB(:,k) = Pb;
        EB(:,k) = mb + Pb * tanh(mb);
        VB(:,k) = Pb + (1-tanh(mb)^2) * Pb^2;
    end
    
    figure(1); clf;
    subplot(2,1,1);
      plot(T,MMB,T,EB);
    subplot(2,1,2);
      plot(T,PPB,T,VB);
    
    rmse_b = sqrt(mean((EB - X).^2))
    
%%
% Determine the quantiles of Benes solution
%
    QQ1 = zeros(size(MMB));
    QQ2 = zeros(size(MMB));
    for k=1:length(MMB)
        pp = dens(xx,MMB(:,k),PPB(:,k));
        cu = cumsum(pp ./ sum(pp));
        ind1 = find(cu >= 2.5/100,1);
        ind2 = find(cu >= 97.5/100,1);
        QQ1(k) = xx(ind1);
        QQ2(k) = xx(ind2);
    end
    
    figure(1); clf;
    plot(T,EB,T,QQ1,'--',T,QQ2,'--');
    
%%
% EKF solution to the same problem
%
    me = m0; 
    Pe = P0;
    MMe = zeros(1,length(Y));
    PPe = zeros(1,length(Y));
    ekf_steps = 5;

    for k=1:length(Y)
        dt2 = dt/ekf_steps;
        for i=1:ekf_steps
            Pe = Pe + dt2 * (2 * (1 - tanh(me)^2) * Pe + 1);
            me = me + dt2 * tanh(me);
        end
        if ~isnan(Y(k))
            me = me + Pe / (R + Pe) * (Y(k) - me);
            Pe = Pe - Pe^2/(R + Pe);
        end
        MME(:,k) = me;
        PPE(:,k) = Pe;
    end

    rmse_e = sqrt(mean((MME - X).^2))
    
    figure(1); clf;
    subplot(2,2,1);
    plot(T,EB,T,QQ1,'--',T,QQ2,'--');
    subplot(2,2,2);
    plot(T,MME,T,MME-1.96*sqrt(PPE),'--',T,MME+1.96*sqrt(PPE),'--');
    
    subplot(2,2,3);
    plot(T,EB,T,MME)
    
    subplot(2,2,4);
    plot(T,VB,T,PPE)
    
%%
% ERTS Smoother
%
    MMS = MME;
    PPS = PPE;
    ms = MME(end);
    Ps = PPE(end);
    for k=length(Y)-1:-1:1
        m = MME(k);
        P = PPE(k);
        dms = dt * tanh(m) + dt * ((1 - tanh(m)^2) * P + 1)/P * (ms - m);
        dPs = dt * 2 * ((1 - tanh(m)^2) * P + 1)/P * Ps - dt;
        ms = ms - dms;
        Ps = Ps - dPs;
        MMS(k) = ms;
        PPS(k) = Ps;
    end

    rmse_s = sqrt(mean((MMS - X).^2))
    
    figure(1); clf;
    plot(T,MMS,T,MMS-1.96*sqrt(PPS),'--',T,MMS+1.96*sqrt(PPS),'--');
    
%%
% Plot the final figures
%
    
    %
    % Benes-Daum
    %
    figure(1); clf; hold on
  
    % Plot mean
    h0 = plot(T,EB,'-k');
    set(h0,'LineWidth',1);
    set(h0,'Color',[0.5 0.5 0.5]);

    % Plot uncetainty
    h1 = fill([T'; flipud(T')], [QQ1'; flipud(QQ2')],1);
    set(h1,'EdgeColor',[.7 .7 .7],'FaceColor',[.7 .7 .7])
    
    % Plot signal
    h2 = plot(T,X,'-k');
    set(h2,'LineWidth',1);
    
    % Plot mean
    h3 = plot(T,EB,'-k');
    set(h3,'LineWidth',1);
    set(h3,'Color',[0.5 0.5 0.5]);
    
    % Plot observations
    h4 = plot(T,Y,'+k');
    
    box on   
    xlabel('Time, $t$'), ylabel('Signal, $x(t)$')
    
    legend([h0, h1, h2, h4],'BD mean','BD 95\% quantiles','Signal','Observations')   

    %
    % Benes-Daum vs EKF
    %
    figure(2); clf; hold on
    
    % Plot BD mean
    h1 = plot(T,EB,'-k');
    set(h1,'LineWidth',1);
    set(h1,'Color',[0.5 0.5 0.5]);

    % Plot EKF mean
    h2 = plot(T,MME,'--k');
    set(h2,'LineWidth',1);
    %set(h2,'Color',[0.0 0.0 0.0]);
    
    % Plot BD uncertainty
    h3 = fill([T'; flipud(T')], [QQ1'; flipud(QQ2')],1);
    set(h3,'EdgeColor',[.7 .7 .7],'FaceColor',[.7 .7 .7])
    
    % Plot BD mean
    h4 = plot(T,EB,'-k');
    set(h4,'LineWidth',1);
    set(h4,'Color',[0.5 0.5 0.5]);
    
    % Plot EKF mean
    h5 = plot(T,MME,'--k');
    set(h5,'LineWidth',1);
    set(h5,'Color',[0.0 0.0 0.0]);
    
    % Plot EKF uncertainty
    h6 = plot(T,MME-1.96*sqrt(PPE),'-.k',T,MME+1.96*sqrt(PPE),'-.k');
    set(h6,'LineWidth',.5);
    
    box on    
    xlabel('Time, $t$'), ylabel('Signal, $x(t)$')
    legend([h1, h2, h3, h6(1)],'BD mean','EKF mean','BD 95\% quantiles','EKF 95\% quantiles')   

    %
    % ERTS
    %
    figure(3); clf; hold on
  
    % Plot mean
    h0 = plot(T,MMS,'-k');
    set(h0,'LineWidth',1);
    set(h0,'Color',[0.5 0.5 0.5]);

    QQ1S = MMS-1.96*sqrt(PPS);
    QQ2S = MMS+1.96*sqrt(PPS);
    
    % Plot uncetainty
    h1 = fill([T'; flipud(T')], [QQ1S'; flipud(QQ2S')],1);
    set(h1,'EdgeColor',[.7 .7 .7],'FaceColor',[.7 .7 .7])
    
    % Plot signal
    h2 = plot(T,X,'-k');
    set(h2,'LineWidth',1);
    
    % Plot mean
    h3 = plot(T,MMS,'-k');
    set(h3,'LineWidth',1);
    set(h3,'Color',[0.5 0.5 0.5]);
    
    % Plot observations
    h4 = plot(T,Y,'+k');
 
    box on
    xlabel('Time, $t$'), ylabel('Signal, $x(t)$')
    
    legend([h0, h1, h2, h4],'ERTS mean','ERTS 95\% quantiles','Signal','Observations')   

