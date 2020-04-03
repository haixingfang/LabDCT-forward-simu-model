classdef ipfColorKey < orientationColorKey
  % defines an orientation mapping based on a certain inverse pole figure
  %   Detailed explanation goes here
  
  properties
    inversePoleFigureDirection
    dirMap 
  end
    
  methods
    
    function oM = ipfColorKey(varargin)
      oM = oM@orientationColorKey(varargin{:});
      
      if isa(oM.CS2,'crystalSymmetry')
        oM.inversePoleFigureDirection = Miller(0,0,1,oM.CS2);
      else
        oM.inversePoleFigureDirection = zvector;
      end
      
      oM.dirMap = getClass(varargin,'directionColorKey',[]);
      if isempty(oM.dirMap), oM.dirMap = HSVDirectionKey(oM.CS1); end
      
    end
    
    function plot(oM,varargin)
      
    
      [~,caxes] = plot(oM.dirMap,'doNotDraw',varargin{:});
      mtexFig = gcm;
      
      mtexTitle(caxes(1),char(oM.inversePoleFigureDirection,'LaTeX'),varargin{:});
            
      name = oM.CS1.pointGroup;
      if ~isempty(oM.CS1.mineral), name = [oM.CS1.mineral ' (' name ')']; end
        
      set(mtexFig.parent,'name',['IPF key for ' name])      
      set(caxes,'tag','ipdf')
      setappdata(caxes,'CS',oM.CS1);
      setappdata(caxes,'inversePoleFigureDirection',oM.inversePoleFigureDirection);
            
      try
        mtexFig.drawNow('figSize',getMTEXpref('figSize'),varargin{:});
      end

    end
        
    function rgb = orientation2color(oM,ori)
    
      if ~(ori.CS.properSubGroup <= oM.CS1)
        warning('The symmetry of the ipf key and the orientations does not fit.')
      end
      
      % compute crystal directions
      ori.CS = oM.CS1;
      h = inv(ori) .* oM.inversePoleFigureDirection;
      
      % colorize fundamental region
      rgb = oM.Miller2Color(h);
      
    end
    
    
    function rgb = Miller2Color(oM,h)
      
      rgb = oM.dirMap.direction2color(h);
      
    end
    
  end  
end
