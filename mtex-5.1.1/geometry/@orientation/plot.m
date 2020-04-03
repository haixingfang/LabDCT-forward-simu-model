function varargout = plot(ori,varargin)
% annotate a orientation to an existing plot
%
% Syntax
%   plot(ori)
%   plot(ori,'label','o1')
%
% Input
%  ori - @orientation
%
% Options
%
% See also
% orientation/scatter orientation/plotPDF orientation/plotODF
% orientation/plotIPDF vector3d/text

[mtexFig,isNew] = newMtexFigure(varargin{:});

if isNew || isappdata(mtexFig.gca,'orientationPlot')
  
  [varargout{1:nargout}] = scatter(ori,varargin{:});
  return;

elseif isappdata(mtexFig.parent,'ODFSections')

  oS = getappdata(mtexFig.parent,'ODFSections');
  [varargout{1:nargout}] = oS.plot(ori,varargin{:});
  return
  
end

if check_option(varargin,'add2all')
  allAxes = mtexFig.children;
else
  allAxes = get_option(varargin,'parent',mtexFig.currentAxes);
end
varargin = delete_option(varargin,{'add2all','parent'});

% plotting
for ax = allAxes(:).'
  switch get(ax,'tag')
  
    case 'pdf' % pole figure annotations
      
      [varargout{1:nargout}] = plotPDF(ori,varargin{:},'parent',ax,'noTitle');
    
    case 'ipdf' % inverse pole figure annotations
      
      [varargout{1:nargout}] = plotIPDF(ori,varargin{:},'parent',ax,'noTitle');
  
    case 'odf' % ODF sections plot
    
      [varargout{1:nargout}] = plotSection(ori,varargin{:},'parent',ax);
    
    otherwise
    
      [varargout{1:nargout}] = scatter(ori,varargin{:},'parent',ax);
  end
end
