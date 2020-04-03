function v = rotate(v,q)
% rotate vector3d by rotation or orientation
%
% Syntax
%   v = rotate(v,20*degree) % rotation about the z-axis
%   rot = rotation('Euler',10*degree,20*degree,30*degree)
%   v = rotate(v,rot)
%
% Description
%  Either |v| or |rot| are single elements or both have the same size. The
%  ouptut |v| will have the same size as the biger of both input arrays.
%
% Input
%  v - @vector3d
%  q - @quaternion
%
% Output
%  r = q * v;

if isnumeric(q), q = axis2quat(zvector,q);end

a = q.a; b = q.b; c = q.c; d = q.d;
[x,y,z] = double(v);

n = b.^2 + c.^2 + d.^2;
s = 2*(x.*b + y.*c + z.*d);

a_2 = 2*a;
a_n  = a.^2 - n;
  
v.x = a_2.*(c.* z - y.*d) + s.*b + a_n.*x;
v.y = a_2.*(d.* x - z.*b) + s.*c + a_n.*y;
v.z = a_2.*(b.* y - x.*c) + s.*d + a_n.*z;
  
if isa(q,'rotation')
  if numel(q.i)>1
    i = logical(q.i);
    v.x(i) = -v.x(i);
    v.y(i) = -v.y(i);
    v.z(i) = -v.z(i);
  elseif q.i
    v.x = -v.x;
    v.y = -v.y;
    v.z = -v.z;
  end
end

% if output has symmetry set it to Miller
if isa(q,'orientation') && isa(q.SS,'crystalSymmetry')
   
  v = Miller(v,q.SS);
  
else % convert to vector3d

  v = vector3d(v);
  
end
