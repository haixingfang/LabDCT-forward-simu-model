function f = dot(sAF1,sAF2,varargin)
%
% syntax
%  f = dot(sAF,v)
%  f = dot(sAF1,sAF2)
%  f = dot(sAF,sVF)
%
% Input
%  sAF,sAF1,sAF2 - @S2AxisFieldHarmonic
%  sVF  - @S2VectorFieldHarmonic
%  v    - @vector3d
%
% Output
%  f - @S2FunHarmonic
%

% first argument should be S2AxisField
if isa(sAF1,'vector3d')
  f = dot(sAF2,sAF1);
  return
end

% the simple case
if isa(sAF2,'vector3d')
  
  [x,y,z] = double(sAF2);
  ghat = [x(:).*x(:),2*x(:).*y(:),y(:).*y(:),2*x(:).*z(:),2*y(:).*z(:),z(:).*z(:)];
 
  f = sAF1.sF;
  f.fhat = sum(f.fhat .* repmat(ghat,size(f.fhat,1),1),2);
    
  f = sqrt(f,varargin{:});
  
  return
end

f = S2FunHarmonic.quadrature(@(v) local_dot(v),varargin{:});

function d = local_dot(v)
  
  f1 = sAF1.sF.eval(v);
  f2 = sAF2.sF.eval(v);

  % add weights for the dot product
  f1(:,[2 4 5]) = 2*f1(:,[2 4 5]);
  
  d = sqrt(sum(f1 .* f2,2));
  
end

end
