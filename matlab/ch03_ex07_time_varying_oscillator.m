%% Example 3.7: Stochastic resonators with time-varying frequency
%
% Copyright: 
%   2018 - Simo Särkkä and Arno Solin
%
% License:
%   This software is provided under the MIT License. See the accompanying 
%   LICENSE file for details.

%%

  % Lock random seed
  if exist('rng') % Octave doesn't have rng
      rng(0,'twister')
  else
      randn('state',1);
  end

  % Time
  t = linspace(0,10,512);

  % Define the sinc function
  sinc = @(x) sin(x)./x;

  % Set time-varying frequency trajectory
  f = 1+.5*abs(t-5);
      
  % Parameters  
  H = [1 0];
  
  figure(1); clf; hold on
  
  for j=1:4
  
      % Allocate space
      x = zeros(2,numel(t));    
      z = x;

      % For harmonics
      for i=1:2

          % Diffusion
          Qc = 10^-i;

          % Initial state
          x(:,1) = [0; 1];

          % Loop through time points
          for k=2:numel(t)    

            % Time delta
            dt = t(k)-t(k-1);  

            % Define F and L
            F = [0 2*pi*i*f(k); -2*pi*i*f(k) 0];  
            L = [0; 1];

            % Solve A and Q
            [A,Q] = lti_disc(F,L,Qc,dt);

            % The next state
            x(:,k) = A*x(:,k-1) + chol(Q,'lower')*randn(size(Q,1),1);

          end

          z = z + H*x;

      end

      plot(t,5*(j-1)+z,'-k','LineWidth',.5)
  
  end
  
  xlabel('Time, $t$')
  xlim([0 10.5])
  ylim([-3 18])
  
  figure(2); clf
    plot(t,f)
    xlabel('Time, $t$')
    ylabel('Frequency [Hz]')

    
    