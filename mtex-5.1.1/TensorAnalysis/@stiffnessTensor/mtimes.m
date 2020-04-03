function sigma = mtimes(C,eps)
% tensor product C * eps

if isa(eps,'stressTensor')

  % TODO: check symmetry
  sigma = stressTensor(...
    EinsteinSum(C,[-1 -2 1 2],eps,[-1 -2]));
  
else
  
  sigma = mtimes@tensor(C,eps);

end
