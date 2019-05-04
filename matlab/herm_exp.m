function [ herm_x, ser, c_list, H_list ] = herm_exp(f, df, d2f, d3f, d4f, d5f, dt, x, x0)
%% HERM_EXP - Hermite expansion
%
% Syntax:
%   [herm_x, ser, c_list, H_list] = herm_exp(f, df, d2f, d3f, d4f, d5f, dt, x, x0)
%   
% Description:
%   Supporting file for 'ch09_ex18_hermite_expansion.m'.
%
% Copyright: 
%   2018 - Simo Särkkä and Arno Solin
%
% License:
%   This software is provided under the MIT License. See the accompanying 
%   LICENSE file for details.

%%

    c_list = {};
    c_list{end+1} = ones(size(x));
    c_list{end+1} = -f.*dt^(1/2) - (2.*f.*df + d2f) .* dt^(3/2) / 4 ...
     - (4.*f.*df.^2 + 4.*f.^2.*d2f + 6.*df.*d2f + 4.*f.*d3f + d4f) .* dt^(5/2) / 24;
    c_list{end+1} = (f.^2 + df).*dt/2 + (6.*f.^2.*df + 4.*df.^2 + 7.*f.*d2f + 2.*d3f).*dt^2/12 ...
      + (28 .* f.^2 .* df.^2 + 28 .* f.^2 .* d3f + 16 .* df.^3 + 16 .* f.^3 .* d2f + 88.*f.*df.*d2f ...
        + 21.*d2f.^2 + 32.*df.*d3f + 16.*f.*d4f + 3.*d5f) .* dt^3 / 96;
    c_list{end+1} = -(f.^3 + 3.*f.*df + d2f).*dt^(3/2)/6 - (12.*f.^3.*df + 28.*f.*df.^2 ...
      + 22.*f.^2.*d2f + 24.*df.*d2f + 14.*f.*d3f + 3.*d4f) .* dt^(5/2) / 48;
    c_list{end+1} = (f.^4 + 6.*f.^2.*df + 3.*df.^2 + 4.*f.*d2f + d3f).*dt^2/24 ...
      + (20.*f.^4.*df + 50.*f.^3.*d2f + 100.*f.^2.*df.^2 + 50.*f.^2.*d3f + 23.*f.*d4f ...
        + 180.*f.*df.*d2f + 40.*df.^3 + 34.*d2f.^2 + 52.*df.*d3f + 4.*d5f).*dt^3/240;
    c_list{end+1} = -(f.^5 + 10.*f.^3.*df + 15.*f.*df.^2 + 10.*f.^2.*d2f ...
      + 10.*df.*d2f + 5.*f.*d3f + d4f).*dt^(5/2)/120;
    c_list{end+1} = (f.^6 + 15.*f.^4.*df + 15.*df.^3 + 20.*f.^3.*d2f + 15.*df.*d3f + 45.*f.^2.*df.^2 ...
      + 10.*d2f.^2 + 15.*f.^2.*d3f + 60.*f.*df.*d2f + 6.*f.*d4f + d5f).*dt^3/720;

    H_list = {};
    H_list{end+1} = ones(size(x));
    H_list{end+1} = -((x-x0)/sqrt(dt));
    H_list{end+1} = -1 + ((x-x0)/sqrt(dt)).^2;
    H_list{end+1} = 3.*((x-x0)/sqrt(dt)) - ((x-x0)/sqrt(dt)).^3;
    H_list{end+1} = 3 - 6.*((x-x0)/sqrt(dt)).^2 + ((x-x0)/sqrt(dt)).^4;
    H_list{end+1} = -15.*((x-x0)/sqrt(dt)) + 10.*((x-x0)/sqrt(dt)).^3 - ((x-x0)/sqrt(dt)).^5;
    H_list{end+1} = -15 + 45.*((x-x0)/sqrt(dt)).^2 - 15.*((x-x0)/sqrt(dt)).^4 + ((x-x0)/sqrt(dt)).^6;

    ser = zeros(size(x));
    for i=1:length(c_list)
        ser = ser + c_list{i} .* H_list{i};
    end

    herm_x = exp(-(x-x0).^2/(2*dt)) / sqrt(2*pi*dt) .* ser;
    
end

