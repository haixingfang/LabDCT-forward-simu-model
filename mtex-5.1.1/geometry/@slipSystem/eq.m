function b = eq(sS1,sS2,varargin)
% ? sS1 == sS2
%
% Input
%  v1, v2 - @vector3d
%
% Output
%  b - boolean
%
% Options
%  antipodal - include antipodal symmetry
%

if sS1.CS ~= sS2.CS || ...
    (length(sS1)>1 && length(sS2)>1 && any(size(sS2)~=size(sS1)))
  b = false;
  return
end

b = sS1.b == sS2.b & sS1.n == sS2.n;
