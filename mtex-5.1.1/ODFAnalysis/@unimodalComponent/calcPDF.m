function pdf = calcPDF(component,h,r,varargin)
% compute the pole density function for a given unimodal component

if nargin > 2 && min(length(h),length(r)) > 0 && ...
    (check_option(varargin,'old') || max(length(h),length(r)) < 10)
  
  pdf = component.psi.RK_symmetrised(...
    component.center,h,r,component.weights,...
    component.CS,component.SS,varargin{:});
  return
  
end

antipodal = extract_option(varargin,'antipodal');
if length(h) == 1 % pole figure

  sh = symmetrise(h,'unique');
  pdf = S2FunHarmonicSym.quadrature(component.center*sh,...
    repmat(component.weights(:),1,length(sh)),component.SS,antipodal{:});

  % convolve with kernel function
  pdf = 4 * pi * conv(pdf,component.psi)./ length(sh);

  % maybe we should evaluate
  if nargin > 2 && isa(r,'vector3d'), pdf = pdf.eval(r); end
  
else % inverse pole figure

  sr = symmetrise(r,component.SS,'unique');
  pdf = S2FunHarmonicSym.quadrature(inv(component.center)*sr,...
    repmat(component.weights(:),1,length(sr)),component.CS,antipodal{:});
  
  % convolve with kernel function
  pdf = 4 * pi * conv(pdf,component.psi) ./ length(sr);

  % maybe we should evaluate
  if isa(h,'vector3d'), pdf = pdf.eval(h); end
  
end

if ~isnumeric(pdf)
  pdf.antipodal = component.antipodal;
end
