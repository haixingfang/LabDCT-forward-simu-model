classdef ChristoffelTensor < tensor
  
  methods
    function sT = ChristoffelTensor(varargin)

      sT = sT@tensor(varargin{:});
      
    end
  end
  
   
  methods (Static = true)
    function C = load(varargin)
      T = load@tensor(varargin{:});
      C = ChristoffelTensor(T);
    end
  end
end