function beta = volumeCompressibility(S)
% computes the volume compressibility of an elasticity tensor
%
% Syntax
%   beta = volumeCompressibility(S)
%
% Input
%  C - elastic stiffness @tensor
%
% Output
%  beta - volume compressibility
%
% Description
%
% $$\beta(x) = S_{iikk}$$
%


% compute tensor product
beta = double(EinsteinSum(S,[-1 -1 -2 -2]));
