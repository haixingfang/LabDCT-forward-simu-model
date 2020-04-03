function map = findByLocation( ebsd, xy )
% select EBSD data by spatial coordinates
%
% Input
%  ebsd - @EBSD
%  xy - list of [x(:) y(:)] coordinates, respectively [x(:) y(:) z(:)]
%
% Output
%  ebsd - @EBSD subset
%
% Example
%   mtexdata fo
%   plotx2east
%   plot(ebsd)
%   p = [10000 5000] %ginput(1)
%   g = findByLocation(ebsd,p)
%
% See also
% EBSD/findByLocation grain2d/findByOrientation

if all(isfield(ebsd.prop,{'x','y','z'}))
  x_D = [ebsd.prop.x,ebsd.prop.y,ebsd.prop.z];
elseif all(isfield(ebsd.prop,{'x','y'}))
  x_D = [ebsd.prop.x,ebsd.prop.y];
else
  error('mtex:findByLocation','no Spatial Data!');
end

delta = 1.5*mean(sqrt(sum(diff(ebsd.unitCell).^2,2)));

x_Dm = x_D-delta;  x_Dp = x_D+delta;

nd = sparse(length(ebsd),size(xy,1));
dim = size(x_D,2);
for k=1:size(xy,1)
  
  candit = find(all(bsxfun(@le,x_Dm,xy(k,:)) & bsxfun(@ge,x_Dp,xy(k,:)),2));
  dist = sqrt(sum(bsxfun(@minus,x_D(candit,:),xy(k,:)).^2,2));
  [dist, i] = min(dist);
  nd(candit(i),k) = 1;
  
end

map = find(any(nd,2));
