classdef (InferiorClasses = {?rotation,?quaternion}) orientation < rotation
% orientation - class representing orientations
%
% This MTEX class represents orientations and misorientations. 
%
% orientation('Euler',phi1,Phi,phi2,cs)  defines an orientation in Euler angles
%

properties
  
  CS = crystalSymmetry('1');   % crystal symmetry 
  SS = specimenSymmetry('1');  % specimen symmetry or crystal symmetry
  antipodal = false
  
end

methods

  function o = orientation(varargin)
    % defines an orientation
    %
    % Syntax
    %   ori = orientation(rot,cs,ss)
    %   ori = orientation('Euler',phi1,Phi,phi2,cs,ss)
    %   ori = orientation('Euler',alpha,beta,gamma,'ZYZ',cs,ss)
    %   ori = orientation('Miller',[h k l],[u v w],cs,ss)
    %   ori = orientation(name,cs,ss)
    %   ori = orientation('axis,v,'angle',omega,cs,ss)
    %   ori = orientation('matrix',A,cs)
    %   ori = orientation('map',u1,v1,u2,v2,cs)
    %   ori = orientation('quaternion',a,b,c,d,cs)
    %
    % Input
    %  rot       - @rotation
    %  cs, ss    - @symmetry
    %  u1,u2     - @Miller
    %  v, v1, v2 - @vector3d
    %  name      - named orientation
    %    currently available:
    %
    %    * 'Cube', 'CubeND22', 'CubeND45', 'CubeRD'
    %    * 'Goss', 'invGoss'
    %    * 'Copper', 'Copper2'
    %    * 'SR', 'SR2', 'SR3', 'SR4'
    %    * 'Brass', 'Brass2'
    %    * 'PLage', 'PLage2', 'QLage', 'QLage2', 'QLage3', 'QLage4'
    %
    % Ouptut
    %  ori - @orientation
    %
    % See also
    % quaternion_index orientation_index
    
    % find and remove symmetries
    args  = cellfun(@(s) isa(s,'symmetry'),varargin,'uniformoutput',true);
    sym = varargin(args);
    varargin(args) = [];
    
    % call rotation constructor
    o = o@rotation(varargin{:});

    if nargin == 0, return;end
    
    % set symmetry
    if ~isempty(varargin) && isa(varargin{1},'orientation')
      o.CS = varargin{1}.CS;
      o.SS = varargin{1}.SS;
      o.antipodal = varargin{1}.antipodal;
    elseif ~isempty(varargin) && ischar(varargin{1}) && strcmpi(varargin{1},'map')
      if isa(varargin{2},'Miller'), o.CS = varargin{2}.CS; end
      if isa(varargin{3},'Miller'), o.SS = varargin{3}.CS; end
    else
      try %#ok<TRYNC>
        a = get_option(varargin,'axis');
        o.CS = a.CS;
        o.SS = a.CS;
      end
    end
    if ~isempty(sym), o.CS = sym{1};end
    if length(sym) > 1, o.SS = sym{2};end
    
    % empty constructor -> done
    if isempty(varargin), return; end
    
    % copy constructor
    switch class(varargin{1})
       
      case 'char'
  
        switch lower(varargin{1})
          
          case 'miller'
            
            if ~isa(o.CS,'crystalSymmetry')
              o.CS = varargin{2}.CS;
            end
            o = orientation(Miller2quat(varargin{2:3},o.CS),o.CS,o.SS);       
  
          otherwise
        
            if exist([varargin{1},'Orientation'],'file') 

              % there is a file defining this specific orientation
              o = eval([varargin{1},'Orientation(o.CS,o.SS)']);
            
            end
        end
    end
    o.antipodal = o.antipodal | check_option(varargin,'antipodal');
    if o.antipodal && o.CS ~= o.SS
      warning('antipodal symmetry is only meaningfull for misorientations between the same phase.')
    end
  end
end

methods (Static = true)
  
  function ori = nan(varargin)    
    s = varargin(cellfun(@isnumeric,varargin));
    q = quaternion.nan(s{:});
    ori = orientation(q,varargin{:});
  end
  
  function ori = id(varargin)
    id = find(~cellfun(@isnumeric,varargin),1)-1;
    q = quaternion.id(varargin{1:id});
    ori = orientation(q,varargin{id+1:end});
  end
  
  function ori = rand(varargin)    
    s = varargin(cellfun(@isnumeric,varargin));
    q = quaternion.rand(s{:});
    ori = orientation(q,varargin{:});
  end
  
  function ori = cube(varargin)
    ori = orientation('Euler',0,0,0,varargin{:});
  end
  
  function ori = cubeND22(varargin)
    ori = orientation('Euler',22*degree,0,0,varargin{:});
  end
  
  function ori = cubeND45(varargin)
    ori = orientation('Euler',45*degree,0,0,varargin{:});
  end
  
  function ori = cubeRD(varargin)
    ori = orientation('Euler',0,22*degree,0,varargin{:});
  end
  
  function ori = goss(varargin)
    ori = orientation('Euler',0,45*degree,0,varargin{:});
  end
  
  function ori = copper(varargin)
    ori = orientation('Euler',90*degree,35*degree,45*degree,varargin{:});
  end
  
  function ori = copper2(varargin)
    ori = orientation('Euler',270*degree,30*degree,45*degree,varargin{:});
  end
  
  function ori = SR(varargin)
    ori = orientation('Euler',53*degree,35*degree,63*degree,varargin{:});
  end
  
  function ori = SR2(varargin)
    ori = orientation('Euler',233*degree,35*degree,63*degree,varargin{:});
  end
  
  function ori = SR3(varargin)
    ori = orientation('Euler',307*degree,35*degree,27*degree,varargin{:});
  end
  
  function ori = SR4(varargin)
    ori = orientation('Euler',127*degree,35*degree,27*degree,varargin{:});
  end
  
  function ori = brass(varargin)
    ori = orientation('Euler',35*degree,45*degree,0*degree,varargin{:});
  end
  
  function ori = brass2(varargin)
    ori = orientation('Euler',325*degree,45*degree,0*degree,varargin{:});
  end
  
  function ori = PLage(varargin)
    ori = orientation('Euler',65*degree,45*degree,0*degree,varargin{:});
  end
  
  function ori = PLage2(varargin)
    ori = orientation('Euler',245*degree,45*degree,0*degree,varargin{:});
  end
  
  function ori = QLage(varargin)
    ori = orientation('Euler',65*degree,20*degree,0*degree,varargin{:});
  end
  
  function ori = QLage2(varargin)
    ori = orientation('Euler',245*degree,20*degree,0*degree,varargin{:});
  end
  
  function ori = QLage3(varargin)
    ori = orientation('Euler',115*degree,160*degree,0*degree,varargin{:});
  end
  
  function ori = QLage4(varargin)
    ori = orientation('Euler',295*degree,160*degree,0*degree,varargin{:});
  end
  
  function ori = invGoss(varargin)
    ori = orientation('Euler',90*degree,45*degree,0*degree,varargin{:});
  end

  function mori = Bain(csGamma,csAlpha)
    %
    % Syntax:
    %   mori = Bain(csGamma,csAlpha)
    %
    % Input
    %  csGamma - parent @crystalSymmetry (cubic fcc)
    %  csAlpha - child @crystalSymmetry (cubic bcc)
    %
    
    mori = orientation('map',Miller(1,0,0,csGamma),Miller(1,0,0,csAlpha),...
      Miller(0,1,0,csGamma,'uvw'),Miller(0,1,1,csAlpha,'uvw'));
  end
  
  function mori = KurdjumovSachs(csGamma,csAlpha)
    %
    % Syntax:
    %   mori = KurdjumovSachs(csGamma,csAlpha)
    %
    % Input
    %  csGamma - parent @crystalSymmetry (cubic fcc)
    %  csAlpha - child @crystalSymmetry (cubic bcc)
    %
    
    mori = orientation('map',Miller(1,1,1,csGamma),Miller(0,1,1,csAlpha),...
      Miller(-1,0,1,csGamma,'uvw'),Miller(-1,-1,1,csAlpha,'uvw'));
  end
  
  function mori = NishiyamaWassermann(csGamma,csAlpha)
    %
    % Syntax:
    %   mori = NishiyamaWassermann(csGamma,csAlpha)
    %
    % Input
    %  csGamma - parent @crystalSymmetry (cubic fcc)
    %  csAlpha - child @crystalSymmetry (cubic bcc)
    %
    
    mori = orientation('map',Miller(1,1,1,csGamma),Miller(0,1,1,csAlpha),...
      Miller(1,1,-2,csGamma,'uvw'),Miller(0,-1,1,csAlpha,'uvw'));
  end
  
  function mori = Pitsch(csGamma,csAlpha)
    %
    % Syntax:
    %   mori = Pitch(csGamma,csAlpha)
    %
    % Input
    %  csGamma - parent @crystalSymmetry (cubic fcc)
    %  csAlpha - child @crystalSymmetry (cubic bcc)
    %
    
    mori = orientation('map',Miller(0,1,0,csGamma),Miller(1,0,1,csAlpha),...
      Miller(1,0,1,csGamma,'uvw'),Miller(-1,1,1,csAlpha,'uvw'));
  
    %mori = orientation('map',Miller(1,1,0,csGamma),Miller(1,1,1,csAlpha),...
    %  Miller(0,0,1,csGamma,'uvw'),Miller(-1,1,0,csAlpha,'uvw'));
  
  end
  
  function mori = GreningerTrojano(csGamma,csAlpha)
    %
    % Syntax:
    %   mori = GreningerTrojano(csGamma,csAlpha)
    %
    % Input
    %  csGamma - parent @crystalSymmetry (cubic fcc)
    %  csAlpha - child @crystalSymmetry (cubic bcc)
    %
    % cube cube
    mori = inv(orientation('Euler',2.7*degree,46.6*degree,7.5*degree,csAlpha,csGamma));

    %mori = orientation('map',Miller(1,1,1,csGamma),Miller(1,1,0,csAlpha),...
    %  Miller(5,12,17,csGamma,'uvw'),Miller(17,17,7,csAlpha,'uvw'));
    
  end
  
    
end

end

