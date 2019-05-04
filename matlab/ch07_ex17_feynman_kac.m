%% Example 7.17: Solution of an elliptic PDE using SDE simulation
%
% Copyright: 
%   2018 - Simo Särkkä and Arno Solin
%
% License:
%   This software is provided under the MIT License. See the accompanying 
%   LICENSE file for details.

%%
% The equation to be considered is
%
%   f_1(x,y) du/dx + f_2 du/dy
%   + q/2 [d^2u(x,y)/dx^2 + d^2u(x,y)/dy^2] - r u(x,y) = 0
%   u(x,0) = psi_x1(x)
%   u(x,L) = psi_x2(x)
%   u(0,y) = psi_y1(y)
%   u(L,y) = psi_y2(y)
% 
% Finite differences approxmation
%
% The discretizations are
%
%  d^2u(x,y)/dx^2 + d^2u(x,y)/dy^2 =~ (u(x+h,y) + u(x,y+h) - 4 u(x,y) + u(x-h,y) + u(x,y-h))/h^2
%  du/dx =~ (u(x+h,y) - u(x-h,y)) / (2h)
%  du/dy =~ (u(x,y+h) - u(x,y-h)) / (2h)
%
% Giving
%
%  f_1(x,y) (u(x+h,y) - u(x-h,y)) / (2h)
%  + f_2(x,y) (u(x,y+h) - u(x,y-h)) / (2h)
%  + q/2 (u(x+h,y) + u(x,y+h) - 4 u(x,y) + u(x-h,y) + u(x,y-h))/h^2
%  -r u(x,y) = 0
%
% i.e.
%
%  [q/2/h^2 + f_1(x,y)/(2h)] u(x+h,y) 
%  + [q/2/h^2 + f_2(x,y)/(2h)] u(x,y+h)
%  - [2q/h^2 + r] u(x,y)
%  + [q/2/h^2 - f_1(x,y)/(2h)] u(x-h,y) 
%  + [q/2/h^2 - f_2(x,y)/(2h)] u(x,y-h) = 0
%
% with discretizations 0:h:(n+1)h the boundary conditions translate to
%
% when x=0:
%  [q/2/h^2 + f_1(x,y)/(2h)] u(x+h,y) 
%  + [q/2/h^2 + f_2(x,y)/(2h)] u(x,y+h)
%  - [2q/h^2 + r] u(x,y)
%  + [q/2/h^2 - f_2(x,y)/(2h)] u(x,y-h)
%  = -[q/2/h^2 - f_1(x,y)/(2h)] psi_y1(y) 
%
% when x=L:
%  [q/2/h^2 + f_2(x,y)/(2h)] u(x,y+h)
%  - [2q/h^2 + r] u(x,y)
%  + [q/2/h^2 - f_1(x,y)/(2h)] u(x-h,y) 
%  + [q/2/h^2 - f_2(x,y)/(2h)] u(x,y-h)
%  = -[q/2/h^2 + f_1(x,y)/(2h)] psi_y2(y) 
%
% when y=0:
%  [q/2/h^2 + f_1(x,y)/(2h)] u(x+h,y) 
%  + [q/2/h^2 + f_2(x,y)/(2h)] u(x,y+h)
%  - [2q/h^2 + r] u(x,y)
%  + [q/2/h^2 - f_1(x,y)/(2h)] u(x-h,y) 
%  = -[q/2/h^2 - f_2(x,y)/(2h)] psi_x1(x)
%
% when y=L:
%  [q/2/h^2 + f_1(x,y)/(2h)] u(x+h,y) 
%  - [2q/h^2 + r] u(x,y)
%  + [q/2/h^2 - f_1(x,y)/(2h)] u(x-h,y) 
%  + [q/2/h^2 - f_2(x,y)/(2h)] u(x,y-h)
%  = -[q/2/h^2 + f_2(x,y)/(2h)] psi_x2(x)
%
% which can be added to the system as additional equations.
% If we let
%
%   U = [u(h,h) u(2h,h) ... u(nh,h) u(h,2h) u(2h,2h) ... u(nh,nh)]^T
%
% then we have for j,i=1,...,n
%
%    u(ih,jh) = U((j-1)*n+i)
% 

  % The model
  lam1 = 0.1;
  lam2 = 0.1;
  L = 1;    
  f1 = @(x,y) -lam1*x;
  f2 = @(x,y) -lam2*y;
  psi_x1 = @(x) sin(2*pi*x);
  psi_x2 = @(x) 1-x;
  psi_y1 = @(y) y;
  psi_y2 = @(y) sin(2*pi*y);
  q = 0.1;
  r = 1;
  
  % Discretization
  n = 50;
  h = L/(n+1);
  
  xgrid = h:h:n*h;
  ygrid = h:h:n*h;
  [YY,XX] = meshgrid(xgrid,ygrid);
  X = XX(:);
  Y = YY(:);

  % Evaluate boundary and f(x,y)
  FF1 = f1(XX,YY);
  FF2 = f2(XX,YY);
  F1 = FF1(:);
  F2 = FF2(:);
  Px1 = psi_x1(xgrid);
  Px2 = psi_x2(xgrid);
  Py1 = psi_y1(ygrid);
  Py2 = psi_y2(ygrid);
  
%  [q/2/h^2 + f_1(x,y)/(2h)] u(x+h,y) 
%  + [q/2/h^2 + f_2(x,y)/(2h)] u(x,y+h)
%  - [2q/h^2 + r] u(x,y)
%  + [q/2/h^2 - f_1(x,y)/(2h)] u(x-h,y) 
%  + [q/2/h^2 - f_2(x,y)/(2h)] u(x,y-h)

  F = spdiags((-2*q/h^2 - r)*ones(n*n,1),0,n*n,n*n);
  b = spalloc(n*n,1,4*n);
  for i=1:n
    for j=1:n
      if j~=1
          % u(x,y-h)
          F((j-1)*n+i,(j-2)*n+i) = q/2/h^2 - F2((j-1)*n+i)/(2*h);
      else
          % boundary y=0
          b((j-1)*n+i) = b((j-1)*n+i) - (q/2/h^2 - F2((j-1)*n+i)/(2*h)) * Px1(i);
      end
      if j~=n
          % u(x,y+h)
          F((j-1)*n+i,(j-0)*n+i) = q/2/h^2 + F2((j-1)*n+i)/(2*h);
      else
          % boundary y=L
          b((j-1)*n+i) = b((j-1)*n+i) - (q/2/h^2 + F2((j-1)*n+i)/(2*h)) * Px2(i);
      end
      if i~=1
          % u(x-h,y)
          F((j-1)*n+i,(j-1)*n+i-1) = q/2/h^2 - F1((j-1)*n+i)/(2*h);
      else
          % boundary x=0
          b((j-1)*n+i) = b((j-1)*n+i) - (q/2/h^2 - F1((j-1)*n+i)/(2*h)) * Py1(j);
      end
      if i~=n
          % u(x+h,y)
          F((j-1)*n+i,(j-1)*n+i+1) = q/2/h^2 + F1((j-1)*n+i)/(2*h);
      else
          % boundary x=L
          b((j-1)*n+i) = b((j-1)*n+i) - (q/2/h^2 + F1((j-1)*n+i)/(2*h)) * Py2(j);
      end
    end
  end
  
%  clf;
%  subplot(2,1,1);
%  spy(F);
  
    U  = F\b;
    UU = reshape(full(U),size(XX));
    
    figure(1); clf;
    surf(XX,YY,UU);
    shading interp;
    if exist('camproj') % Octave doesn't have this
        camproj perspective;
    end
    colormap winter;
    colorbar;
    
    xlabel('x');
    ylabel('y');
    
%%
% Feynman-Kac approximation
%
    % Lock random seed
    if exist('rng') % Octave doesn't have rng
      rng(1);
    else
      randn('state',2);
    end
    
    dt = 0.01;
    a1 = exp(-lam1*dt);
    a2 = exp(-lam2*dt);
    sq1 = sqrt(q/(2*lam1)*(1 - exp(-2*lam1*dt)));
    sq2 = sqrt(q/(2*lam2)*(1 - exp(-2*lam2*dt)));
    
    UFK = zeros(size(XX));
    NN  = zeros(size(XX));
    
    squ_x = [0 L L 0];
    squ_y = [0 0 L L];
    
    stored_SX1 = [];
    stored_SX2 = [];
    stored_whos_in = [];
    
    nmc = 100;
    for mc=1:nmc
        SX1 = XX(:);
        SX2 = YY(:);

        count = 0;
        whos_in = true(1,length(SX1));
        num_in  = sum(whos_in);
        while num_in > 0
            SX1(whos_in) = a1 * SX1(whos_in) + sq1 * randn(size(SX1(whos_in)));
            SX2(whos_in) = a2 * SX2(whos_in) + sq2 * randn(size(SX2(whos_in)));

            ind = find(whos_in(:) & (SX1(:) < 0));
            whos_in(ind) = false;
            UFK(ind) = UFK(ind) + psi_y1(SX2(ind));
            NN(ind) = NN(ind) + 1;

            ind = find(whos_in(:) & (SX1(:) > L));
            whos_in(ind) = false;
            UFK(ind) = UFK(ind) + psi_y2(SX2(ind));
            NN(ind) = NN(ind) + 1;

            ind = find(whos_in(:) & (SX2(:) < 0));
            whos_in(ind) = false;
            UFK(ind) = UFK(ind) + psi_x1(SX1(ind));
            NN(ind) = NN(ind) + 1;

            ind = find(whos_in(:) & (SX2(:) > L));
            whos_in(ind) = false;
            UFK(ind) = UFK(ind) + psi_x2(SX1(ind));
            NN(ind) = NN(ind) + 1;

            num_in = sum(whos_in);

            count = count + 1;
            
            if mc == 1 && count == 200
                stored_whos_in = whos_in;
                stored_SX1 = SX1;
                stored_SX2 = SX2;
                clf;
                plot(SX1(:),SX2(:),'.',squ_x,squ_y,'k');
                pause(1);
            end

            if rem(count,100) == 0

                subplot(2,2,1);
                pcolor(XX,YY,UU);
                colorbar;

                subplot(2,2,2);
                pcolor(XX,YY,UFK ./ NN);
                colorbar;

                subplot(2,2,3);
                plot(SX1(:),SX2(:),'.',squ_x,squ_y,'k');

                subplot(2,2,4);
                pcolor(XX,YY,abs(UU - UFK ./ NN));
                title(sprintf('MC iteration %d/%d',mc,nmc));
                colorbar;

                drawnow;
                fprintf('count = %d, num_in = %d\n',count,num_in);
            end

        end
    end
    
    % Consider caching the results for later use
    % save fkac1
    
%%
% Plot the results
%
    figure(1); clf;
    imagesc(XX(:,1),YY(1,:),UU);
    caxis([-1 1])
    shading flat
    colormap gray
    colorbar
    axis xy
    
    figure(2); clf;
    imagesc(XX(:,1),YY(1,:),UFK ./ NN);
    caxis([-1 1])    
    shading flat
    colormap gray
    colorbar;
    axis xy

    figure(3); clf;
    plot(stored_SX1(stored_whos_in),stored_SX2(stored_whos_in),'kx',...
         stored_SX1(~stored_whos_in),stored_SX2(~stored_whos_in),'ko',squ_x,squ_y,'k');
    axis xy
     
    figure(4); clf;
    imagesc(XX(:,1),YY(1,:),abs(UU - UFK ./ NN));
    shading flat
    colormap gray
    colorbar
    axis xy
    
    nmc
    
    err = UU - UFK ./ NN;
    rmse = mean(err(:).^2)