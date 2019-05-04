%% Example 3.10: Stochastic spring model
%
% Copyright: 
%   2018 - Simo Särkkä and Arno Solin
%
% License:
%   This software is provided under the MIT License. See the accompanying 
%   LICENSE file for details.

%%
% Simulate model
%
if exist('rng') % Octave doesn't have rng
    rng(1,'twister');
else
    randn('state',1);
end

t1 = 10;
dt = 0.01;
steps = round(t1/dt);
g = 1;
v = 2;
q = 0.02;
x0 = [1;0];
nmc = 10;

F = [0 1; -v^2 -g];
L = [0;1];
[A,Q] = lti_disc(F,L,q,dt);

M = [];
T = [];
m = x0;
t = 0;
for k=1:steps
    m = A*m;
    t = t + dt;
    
    M = [M m];
    T = [T t];
end

X1 = [];
X2 = [];
for n=1:nmc
    X = [];
    x = x0;
    for k=1:steps
        x = A*x + chol(Q,'lower')*randn(size(Q,1),1);
        X = [X x];
    end
    X1 = [X1;X(1,:)];
    X2 = [X2;X(2,:)];
end

figure(1); clf; hold on
  plot(T,M(1,:),'-k','linewidth',1)
  plot(T,X1','-','color',[.5 .5 .5])
  plot(T,M(1,:),'-k','linewidth',1)
  legend('Deterministic solution','Realizations of the SDE');
  xlabel('Time, $t$'); ylabel('Displacement, $x(t)$')
  box on

%%
% Mean and covariance trajectories
%

m = x0;
P = zeros(2,2);
MM = [];
VV = [];
for k=1:steps
    m = A*m;
    P = A*P*A' + Q;
    
    MM = [MM m];
    VV = [VV diag(P)];
end

figure(2); clf; hold on
  
  % Draw mean solution
  plot(T,MM(1,:),'-k','linewidth',1);

  % Draw 95% shaded quantiles
  ind = [1:10:steps steps];
  hp = fill([T(ind)'; flipud(T(ind)')], ...
            [MM(1,ind)'-1.96*sqrt(VV(1,ind))'; flipud(MM(1,ind)'+1.96*sqrt(VV(1,ind))')],1);
  set(hp,'EdgeColor',[.9 .9 .9],'FaceColor',[.9 .9 .9])

  % Draw realization
  plot(T,X1','-','color',[0.5 0.5 0.5])
    
  % Draw mean on top of everything
  plot(T,MM(1,:),'-k','linewidth',1);
  
  legend('Mean solution','95\% quantiles','Realizations of the SDE');
  xlabel('Time, $t$'); ylabel('Displacement, $x(t)$')
  box on
  set(gca,'Layer','Top')
