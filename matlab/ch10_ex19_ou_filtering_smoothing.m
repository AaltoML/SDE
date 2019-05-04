%% Examples 10.19, 10.21, 10.29, 10.33: Kalman filtering/smoothing of OU
%
% Copyright: 
%   2018 - Simo Särkkä and Arno Solin
%
% License:
%   This software is provided under the MIT License. See the accompanying 
%   LICENSE file for details.

%%
% This is a double check of the discretization
%
    dt = 0.1;
    lambda = 0.3;
    q = pi;
    
    a1 = exp(-lambda*dt);
    S1 = (q/(2*lambda))*(1 - exp(-2*lambda*dt));
    
    [a2,S2] = lti_disc(-lambda,1,q,dt)

%%
% Simulate data
%
    % Lock random seed
    if exist('rng') % Octave doesn't have rng
      rng(103,'twister');
    else
      randn('state',0);
    end
    
    dt = 0.01;
    lambda = 0.5;
    q = 1;
    R = 0.1;
    steps = 300;

    X = zeros(1,steps);
    Y = zeros(1,steps);
    T = zeros(1,steps);
    
    m0 = 0;
    P0 = 1;
    
    x = m0 + sqrt(P0) * randn;
    t = 0;
    A = exp(-lambda*dt);
    Q = (q - q*exp(-2*dt*lambda))/(2*lambda);
    for k=1:steps
        x = A*x + sqrt(Q) * randn;
        if rem(k,10) == 0
            y = x + sqrt(R) * randn;
        else
            y = NaN;
        end
        t = t + dt;
        X(k) = x;
        Y(k) = y;
        T(k) = t;
    end
    
    figure(1); clf
    plot(T,X,T,Y,'+');
    
%%
% Kalman filter (implement in 3 ways for debugging)
%

    MM = zeros(1,length(Y));
    PP = zeros(1,length(Y));

    % This is the debug
    m = m0;
    P = P0;
    for k=1:length(Y)
        m = A*m;
        P = A*P*A' + Q;
        if ~isnan(Y(k))
            S = P + R;
            K = P / S;
            m = m + K * (Y(k) - m);
            P = P - K * S * K';
        end
    end
    m1 = m;
    P1 = P;

    % This is corresponding to Example 10.19
    m = m0;
    P = P0;
    eu_steps = 100;
    for k=1:length(Y)
        dt2 = dt / eu_steps;
        for j=1:eu_steps
            m = m + (-lambda * m) * dt2;
            P = P + (-2 * lambda * P + q) * dt2;
        end
        
        if ~isnan(Y(k))
            m = m + P / (P + R) * (Y(k) - m);
            P = P - P^2 / (P + R);
        end
        MM(k) = m;
        PP(k) = P;
    end
    m2 = m;
    P2 = P;
    
    % This is the actual discretized one (Example 10.21)
    m = m0;
    P = P0;
    for k=1:length(Y)
          m = exp(-lambda * dt) * m;
          P = exp(-2 * lambda * dt) * P + (q / (2 * lambda)) * (1 - exp(-2 * lambda * dt));
          if ~isnan(Y(k))
              m = m + P / (P + R) * (Y(k) - m);
              P = P - P^2 / (P + R);
          end
          MM(k) = m;
          PP(k) = P;
    end

    [m m1 m2]
    [P P1 P2]
    
    figure(1);
    plot(T,MM,T,Y,'+',T,X,'--',T,MM-1.96*sqrt(PP),'-.',T,MM+1.96*sqrt(PP),'-.');
    
%%
% RTS smoother, again implement in 3 different ways
%
 
    % This is a reference solution
    kf_m = MM;
    kf_P = reshape(PP,[1 1 length(PP)]);
    ms = kf_m(:,end);
    Ps = kf_P(:,:,end);
    rts_m = zeros(size(m,1),size(Y,2));
    rts_P = zeros(size(P,1),size(P,2),size(Y,2));
    rts_m(:,end) = ms;
    rts_P(:,:,end) = Ps;
    for k=size(kf_m,2)-1:-1:1
        mp = A*kf_m(:,k);
        Pp = A*kf_P(:,:,k)*A'+Q;
        Ck = kf_P(:,:,k)*A'/Pp;
        ms = kf_m(:,k) + Ck*(ms - mp);
        Ps = kf_P(:,:,k) + Ck*(Ps - Pp)*Ck';
        rts_m(:,k) = ms;
        rts_P(:,:,k) = Ps;
    end
    ms1 = ms;
    Ps1 = Ps;
    
    % This is corresponding to Example 10.29
    MMS2 = MM;
    PPS2 = PP;
    ms = MM(end);
    Ps = PP(end);
    for k=length(Y)-1:-1:1
        MMf = zeros(1,eu_steps);
        PPf = zeros(1,eu_steps);
        dt2 = dt / eu_steps;
        m = MM(k);
        P = PP(k);
        for j=1:eu_steps
            MMf(j) = m;
            PPf(j) = P;
            m = m + (-lambda * m) * dt2;
            P = P + (-2 * lambda * P + q) * dt2;
        end

        for j=eu_steps:-1:1
            dms = -lambda * ms + (q / PPf(j)) * (ms - MMf(j));
            dPs = 2 * (-lambda + q / PPf(j)) * Ps - q;
            ms = ms - dms * dt2;
            Ps = Ps - dPs * dt2;
        end

        MMS2(k) = ms;
        PPS2(k) = Ps;
    end
    ms2 = ms;
    Ps2 = Ps;
    
    %
    % This is the actual one from Example 10.33
    %
    MMS = MM;
    PPS = PP;
    
    ms = MM(end);
    Ps = PP(end);
    for k=length(Y)-1:-1:1
        mp = exp(-lambda * dt) * MM(k);
        Pp = exp(-2 * lambda * dt) * PP(k) + (q / (2 * lambda)) * (1 - exp(-2 * lambda * dt));
        ms = MM(k) + (PP(k) * exp(-lambda * dt) / Pp) * (ms - mp);
        Ps = PP(k) + (PP(k) * exp(-lambda * dt) / Pp)^2 * (Ps - Pp);

        MMS(k) = ms;
        PPS(k) = Ps;
    end
    
    [ms ms1 ms2]
    
    figure(1); clf;
    subplot(2,1,1);
    plot(T,MMS,'r-',T,rts_m,'b--',T,MMS2,'g:')
    
    subplot(2,1,2);
    plot(T,PPS,'r-',T,squeeze(rts_P),'b--',T,PPS2,'g:')
    
%%
% Plot the final figures
%
    figure(1); clf; hold on

    % Plot KF mean
    h1 = plot(T,MM,'-k');
    set(h1,'LineWidth',1);
    set(h1,'Color',[0.5 0.5 0.5]);
    
    % Plot KF uncertainty
    QQ1 = MM - 1.96 * sqrt(PP); 
    QQ2 = MM + 1.96 * sqrt(PP); 
    h2 = fill([T'; flipud(T')], [QQ1'; flipud(QQ2')],1);
    set(h2,'EdgeColor',[.7 .7 .7],'FaceColor',[.7 .7 .7])
    
    % Plot signal
    h3 = plot(T,X,'-k');
    set(h3,'LineWidth',1);
    
    % Plot KF mean
    h4 = plot(T,MM,'-k');
    set(h4,'LineWidth',1);
    set(h4,'Color',[0.5 0.5 0.5]);
    h5 = plot(T,Y,'+k');
    box on
    xlabel('Time, $t$'), ylabel('Signal, $x(t)$')
    
    legend([h1, h2, h3, h5],'KF mean','KF 95\% quantiles','Signal','Observations')   
    
    figure(2); clf; hold on

    % Plot RTS mean
    h1 = plot(T,MMS,'-k');
    set(h1,'LineWidth',1);
    set(h1,'Color',[0.5 0.5 0.5]);
    
    % Plot RTS uncertainty
    QQ1 = MMS - 1.96 * sqrt(PPS); 
    QQ2 = MMS + 1.96 * sqrt(PPS); 
    h2 = fill([T'; flipud(T')], [QQ1'; flipud(QQ2')],1);
    set(h2,'EdgeColor',[.7 .7 .7],'FaceColor',[.7 .7 .7])
    
    % Plot signal
    h3 = plot(T,X,'-k');
    set(h3,'LineWidth',1);
    
    % Plot RTS mean
    h4 = plot(T,MMS,'-k');
    set(h4,'LineWidth',1);
    set(h4,'Color',[0.5 0.5 0.5]);
    
    h5 = plot(T,Y,'+k');

    box on
    xlabel('Time, $t$'), ylabel('Signal, $x(t)$')
    
    legend([h1, h2, h3, h5],'RTS mean','RTS 95\% quantiles','Signal','Observations')   
