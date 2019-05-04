%% Example 8.24: Exact simulation of sine diffusion (EA1)
%
% Copyright: 
%   2018 - Simo Särkkä and Arno Solin
%
% License:
%   This software is provided under the MIT License. See the accompanying 
%   LICENSE file for details.

%%

  % Parameters
  T = 1.5;
  n = 100000;

  % Define the SDE model
  f = @(x,t) sin(x);
  L = @(x,t) 1;
  
  % Define d/dx f(x)
  df = @(x) cos(x);
 
  % Define psi(x) = int_0^x sin(u) du
  psi = @(x) 1-cos(x);
  
  % Bounds: k1 < .5*(f(x)^2 + df(x)) < k2
  k1 = -0.5;
  k2 =  0.6250;
  M = (k2-k1);

  % Initialize h(x) = exp(psi(x)?x2/2T)/c
  c = integral(@(u) exp(psi(u)-u.^2/(2*T)),-inf,inf);
  h = @(u) 1/c*exp(psi(u)-u.^2/(2*T));

  % Initialize: phi(x) =  .5*f(x) + .5*df(x)?k1
  phi = @(x) 1/2*(f(x,0).^2+df(x))-k1;
  
  % Choose an epsilon for rejection sampling
  epsilon = 1/70;
  
%%  

  % Lock seed
  if exist('rng') % Octave doesn't have rng
    rng(0,'twister')
  else
    randn('state',0);
    rand('state',0);
  end


% Allocate space for samples
xea = zeros(n,1);

for i=1:n

  % Proposal skeleton
  while (true)

    % Step #1
    lambda = T*M;
    
    if exist('poissrnd')
        tau = poissrnd(lambda);
    else
        % Knuth's algorithm for generating tau = Poisson(lambda)
        ell = exp(-lambda);
        tau = 0;
        p = 1;
        while p > ell
            tau = tau + 1;
            u = rand;
            p = p * u;
        end
        tau = tau - 1;
    end
    
    % Step #2
    t = T*rand(1,tau);
    x = M*rand(1,tau);
  
    % Step #3
    beta = zeros(1,2+tau);
    beta(1) = 0;
    beta(2) = rejectionsample(h,epsilon);
    tb = [0 T t];
  
    % The rest of the points by Brownian bridges
    for j=3:numel(beta)

      % Lower
      ind = find(tb(1:j-1)<tb(j)); 
      [t1,i1] = max(tb(ind));
      b1 = beta(ind); b1 = b1(i1);
     
      % Upper
      ind = find(tb(1:j-1)>tb(j)); 
      [t2,i2] = min(tb(ind));
      b2 = beta(ind); b2 = b2(i2);
     
      % The Brownian bridge
      mu = b1+(b2-b1)/(t2-t1)*(tb(j)-t1); 
      sigma2 = (t2-tb(j))*(tb(j)-t1)/(t2-t1);     
      beta(j) = mu + sqrt(sigma2)*randn(1);
     
    end

    % Step #4
    N = sum(x<phi(beta(3:end)));
  
    % Step #5: If $N = 0$, go to step 6 and otherwise go to step 1.
    if N==0, break; end
  
  end
  
  % Save sample
  xea(i) = beta(2);
  
  if rem(i,floor(n/100))==0, fprintf('%i/%i\n',i,n), end    
    
end
  
%% Visualize

figure(1); clf; hold on
  [N,edges]=hist(xea,100);
  stairs(edges-(edges(2)-edges(1))/2,N/sum(N)/(edges(2)-edges(1)),'-k');
  
  % Format plot
  xlim([-5 5])
  xlabel('$x$')
  

  