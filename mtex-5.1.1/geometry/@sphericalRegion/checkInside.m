function inside = checkInside(sR,v,varargin)
% check for points to be inside the spherical region
  
% for antipodal symmetry check v and -v
if (sR.antipodal || check_option(varargin,'antipodal')) && ...
    ~check_option(varargin,'removeAntipodal')
  sR.antipodal = false;
  inside = checkInside(sR,v) | checkInside(sR,-v);
  return
end

if check_option(varargin,'removeAntipodal')
  v.antipodal = false;
end

% verify all conditions are satisfies
inside = true(size(v));
v = normalize(vector3d(v));
for i = 1:length(sR.N)
  inside = inside & (dot(subSet(sR.N,i),v) >= sR.alpha(i)-1e-4);
end
 
end
