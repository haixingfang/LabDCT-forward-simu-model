function h = plotAxisDistribution(obj,varargin)
% plot axis distribution
%
% Syntax
%
%   plotAxisDistribution(cs)        % random axis distribution
%   plotAxisDistribution(cs1,cs2)   % random misorientation axis distribution
%   plotAxisDistribution(mori)      % axes in crystal coordinates
%   plotAxisDistribution(mori,cs1)  % axes in crystal coordinates
%   plotAxisDistribution(ori1,ori2) % axes in specimen coordinates
%   plotAxisDistribution(odf)
%   
% Input
%  odf - @ODF
%  mori - @misorientation
%  ori1,ori2 - @orientation
%  cs1, cs2 - @crystalSymmetry
%
% Options
%  RESOLUTION - resolution of the plots
%
% Flags
%  antipodal - include [[AxialDirectional.html,antipodal symmetry]]
%  complete  - plot entire (hemi)--sphere
%
% See also
% S2Grid/plot savefigure Plotting Annotations_demo ColorCoding_demo PlotTypes_demo
% SphericalProjection_demo

[mtexFig,isNew] = newMtexFigure(varargin{:});

if isa(obj,'symmetry')
  cs1 = obj;
  if nargin > 1 && isa(varargin{1},'symmetry')
    cs2 = varargin{1};
  else
    cs2 = specimenSymmetry;
  end
elseif isa(obj,'orientation')
  
  obj = axis(obj,varargin{:});
  
else
  cs1 = obj.CS;
  cs2 = obj.SS;
end


if isa(obj,'vector3d')

  h = plot(obj,'symmetrised','FundamentalRegion',varargin{:});
  
else

  if isa(obj,'orientation') && isa(varargin{1},'crystalSymmetry')
    dcs = varargin{1};
  else
    dcs = calcDisjoint(cs1,cs2,varargin{:});
  end
    
  % plotting grid
  sR = fundamentalSector(dcs,varargin{:});
  if isa(obj,'symmetry')
    h = plotS2Grid(sR,'resolution',.5*degree,varargin{:});
  else
    h = plotS2Grid(sR,'resolution',2.5*degree,varargin{:});
  end

  % plot  
  varargin = [delete_option(varargin,'complete'),dcs.plotOptions];
  density = pos(calcAxisDistribution(obj,h,varargin{:}));
  h = smooth(h,density,dcs,varargin{:});

end
  
if isNew % finalize plot
  name = inputname(1);
  set(gcf,'Name',['Axis Distribution of ',name]);
  mtexFig.drawNow('figSize',getMTEXpref('figSize'),varargin{:}); 
end

if nargout == 0, clear h; end

  function dcs = calcDisjoint(cs1,cs2,varargin)
    
    dcs = disjoint(cs1.properGroup,cs2.properGroup);
    if check_option(varargin,'antipodal') || ...
        ((isa(obj,'orientation') || isa(obj,'ODF')) && obj.antipodal)
      dcs = dcs.Laue;
    end
  end

end

function d = pos(d)

if min(d(:)) > -1e-5, d = max(d,0);end

end
