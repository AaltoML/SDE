function newquiver(X,Y,U,V,varargin)
% NEWQUIVER - Plot velocity vectors as arrows
% 
% Syntax:
%   NEWQUIVER(X,Y,U,V,options)
%
% In:
%   X       - Grid of x coordinates
%   Y       - Grid of y coordinates
%   U       - Vector direction (x component)
%   V       - Vector direction (y component)
%   options - name-value pairs (see below)
%
% Out:
%   (none)
%
% Description:
%   NEWQUIVER(X,Y,U,V) plots velocity vectors as arrows with components 
%   (u,v) at the points (x,y).  The matrices X,Y,U,V must all be the same 
%   size and contain corresponding position and velocity components 
%   (X and Y can also be vectors to specify a uniform grid). 
%   NEWQUIVER automatically scales the arrows to fit within the grid.
%
%   Additional option-value pairs can be given.
%   Accepted options/values are:
%    * scale     : A two-element vector for scaling the shape 
%    * color     : Fill color (rgbkcy...) (default 'k': black)
%    * EdgeColor : Edge color of shape  (default 'k': black)
%    * LineWidth : Line width (default 0.5)
%    * X         : Arrow shape as a two-column vector (default: arrow)
%
% See also:
%    quiver
%
% Copyright: 
%   2012-2018 - Simo Särkkä and Arno Solin
%
% License:
%   This software is provided under the MIT License. See the accompanying 
%   LICENSE file for details.

%% Define arrow shape

  % Define arrow shape
  arrow.X = [8.7185878,  4.0337352 % Upper right corner
            -2.2072895,  0         % Arrow tip
             8.7185878, -4.0337352 % Lower right corner
             6.9831476, -1.6157441 % Inner lower corner
             6.9831476, -0.8       % Inner lower shaft start
             28,        -0.8       % Lower shaft end
             28,         0.8       % Upper shaft end
             6.9831476,  0.8       % Inner upper shaft start
             6.9831476,  1.6157441 % Inner upper corner
             8.7185878,  4.0337352]; % Upper right corner
  
   % Flip arrow left to right
   arrow.X(:,1) = -arrow.X(:,1);
   
   % Normalize length to one and center
   arrow.X(:,1) = arrow.X(:,1)-min(arrow.X(:,1));
   arrow.X = arrow.X/max(arrow.X(:,1));
   arrow.X(:,1) = arrow.X(:,1)-.5;


%% Set default options

  % Default options structure
  defaultopt.scale       = .5*[min(X(X(:)>min(X(:))))-min(X(:)) ...
                               min(Y(Y(:)>min(Y(:))))-min(Y(:))]; 
  defaultopt.color       = 'k';
  defaultopt.EdgeColor   = 'k';
  defaultopt.LineWidth   = .5;
  defaultopt.X           = arrow.X;
 
  
%% Get options

  % If there are options submitted
  if length(varargin) > 1
  
    % Options
    options = cell2struct(varargin(2:2:end),varargin(1:2:end),2);
  
    % Merge the default and submitted options
    fnames = fieldnames(defaultopt);
    for i=1:length(fnames)
      %val = getfield(options,fnames{i});
      if ~isfield(options,fnames{i})
        options=setfield(options,fnames{i},getfield(defaultopt,fnames{i}));
      end
    end
  else
    options = defaultopt;  
  end
  
  % Make sure there is a scaling factor in both directions
  if numel(options.scale)~=2, options.scale = [1 1]*options.scale(1); end
  
   
%% Get arrow directions 
   
  % Calculate arrow direction
  TH = atan2(V,U);

  hold on
  
  for i=1:numel(X)
   
    % Roatate
    Z = [cos(TH(i)) -sin(TH(i)); sin(TH(i)) cos(TH(i))]*options.X';
   
    % Scale
    Z = bsxfun(@times,Z,options.scale(:));
    
    % Move
    Z = bsxfun(@plus,Z,[X(i); Y(i)]);
   
    % Plot
    patch(Z(1,:),Z(2,:),1,'FaceColor',options.color, ...
        'EdgeColor',options.EdgeColor, ...
        'LineWidth',options.LineWidth,'LineStyle','-')
   
  end
