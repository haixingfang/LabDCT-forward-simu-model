classdef ipfTSLKey < ipfColorKey
  % orientation mapping as it is used by TSL and HKL software
      
  methods
    function oM =ipfTSLKey(varargin)

      oM = oM@ipfColorKey(varargin{:}); 
      oM.dirMap = TSLDirectionKey(oM.CS1);
      
    end
  end
    
end
