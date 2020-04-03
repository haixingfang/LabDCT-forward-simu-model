function [h,ax] = scatter(v,varargin)
%
% Syntax
%   scatter(v)              % plot the directions v
%   scatter(v,data)         % colorize directions according to data
%   scatter(v,'label',text) % plot text below markers
%   scatter(v,'label',text,'textaboveMarker') % plot text above markers
%
% Input
%  v     - @vector3d
%  data  - double
%  rgb   - a list of rgb color values
%
% Options
%  Marker            - 's','o','diamont','p'
%  MarkerFaceColor   - 'r','g','w','k','b'
%  MarkerEdgeColor   - 'r','g','w','k','b'
%  MarkerColor       - shortcut for the above two
%  MarkerSize        - size of the markers in pixel
%  DynamicMarkerSize - scale marker size when plot is resized
%
% Output
%
% See also
% vector3d/text

% initialize spherical plots
opt = delete_option(varargin,...
  {'lineStyle','lineColor','lineWidth','color','edgeColor','MarkerSize','Marker','MarkerFaceColor','MarkerEdgeColor','MarkerColor'});
sP = newSphericalPlot(v,opt{:},'doNotDraw');
varargin = delete_option(varargin,'parent');

h = [];

for i = 1:numel(sP)

  % project data
  [x,y] = project(sP(i).proj,v,varargin{:});

  % check that there is something left to plot
  if all(isnan(x) | isnan(y)), continue; end
    
  % add some nans if lines are plotted
  if check_option(varargin,'edgecolor')
    
    % find large gaps
    d = sqrt(diff(x([1:end,1])).^2 + diff(y([1:end,1])).^2);
    ind = find(d > diff(sP(i).bounds([1,3])) / 5);
    
    % and fill the gaps with nans
    for k = 1:numel(ind)
      x = [x(1:ind(k)+k-1);nan;x(ind(k)+k:end)];
      y = [y(1:ind(k)+k-1);nan;y(ind(k)+k:end)];
    end
  end
  
  % default arguments
  patchArgs = {'parent',sP(i).hgt,...
    'vertices',[x(:) y(:)],...
    'faces',1:numel(x),...
    'facecolor','none',...
    'edgecolor','none',...
    'marker','o',...
    };

  % markerSize
  if ~check_option(varargin,{'scatter_resolution','MarkerSize'})
    res = max(v.resolution,1*degree);
  else
    res = get_option(varargin,'scatter_resolution',1*degree);
  end
  MarkerSize  = get_option(varargin,'MarkerSize',min(8,50*res));
  
  patchArgs = [patchArgs,{'MarkerSize',MarkerSize}]; %#ok<AGROW>

  % dynamic markersize
  if check_option(varargin,'dynamicMarkerSize') || ...
      (~check_option(varargin,'MarkerSize') && length(v)>20)
    patchArgs = [patchArgs {'tag','dynamicMarkerSize','UserData',MarkerSize}]; %#ok<AGROW>
  end
    
  % ------- colorcoding according to the first argument -----------
  if ~isempty(varargin) && isa(varargin{1},'crystalShape')
    
    h(i) = plot(x,y,zUpDown * varargin{1}.diameter,varargin{1},'parent', sP(i).hgt,varargin{2:end});
    sP(i).updateBounds(0.1);
  
  elseif ~isempty(varargin) && isnumeric(varargin{1}) && ~isempty(varargin{1})
      
    % extract colorpatchArgs{3:end}coding
    cdata = varargin{1};
    if numel(cdata) == length(v)
      cdata = reshape(cdata,[],1);
      sP(i).updateMinMax(cdata);
    else
      cdata = reshape(cdata,[],3);
    end

    if numel(MarkerSize) > 1
      
      h(i) = optiondraw(scatter(x(:),y(:),MarkerSize(:),cdata,'filled',...
        'parent',sP(i).hgt),varargin{:}); %#ok<AGROW>

      set(get(get(h(i),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
            
    else % draw patches
    
      h(i) = optiondraw(patch(patchArgs{:},...
        'facevertexcdata',cdata,...
        'markerfacecolor','flat',...
        'markeredgecolor','flat'),varargin{2:end}); %#ok<AGROW>
      set(get(get(h(i),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    end
      
  else % --------- colorcoding according to nextStyle -----------------
      
    % get color
    if check_option(varargin,{'MarkerColor','MarkerFaceColor'})
      mfc = get_option(varargin,'MarkerColor','none');
      mfc = get_option(varargin,'MarkerFaceColor',mfc);
    else % cycle through colors
      [ls,mfc] = nextstyle(sP(i).ax,true,true,~ishold(sP(i).ax)); %#ok<ASGLU>
    end
    mec = get_option(varargin,'MarkerEdgeColor',mfc);
  
    % draw patches
    if numel(MarkerSize) > 1
      
      h(i) = optiondraw(scatter(x(:),y(:),MarkerSize(:),'parent',sP(i).hgt,...
        'MarkerFaceColor',mfc,'MarkerEdgeColor',mec),varargin{:}); %#ok<AGROW>      
    
    else
       
      h(i) = optiondraw(patch(patchArgs{:},...
        'MarkerFaceColor',mfc,...
        'MarkerEdgeColor',mec),varargin{:}); %#ok<AGROW>
      % remove from legend
      set(get(get(h(i),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    
      % since the legend entry for patch object is not nice we draw an
      % invisible scatter dot just for legend
      if check_option(varargin,'DisplayName')
        holdState = get(sP(i).ax,'nextPlot');
        set(sP(i).ax,'nextPlot','add');
        optiondraw(scatter([],[],'parent',sP(i).ax,'MarkerFaceColor',mfc,...
          'MarkerEdgeColor',mec),varargin{:});
        set(sP(i).ax,'nextPlot',holdState);
      end
    end
  end

  % set resize function for dynamic marker sizes
  try
    hax = handle(sP(i).ax);
    hListener(1) = handle.listener(hax, findprop(hax, 'Position'), ...
      'PropertyPostSet', {@localResizeScatterCallback,sP(i).ax});
    % save listener, otherwise  callback may die
    setappdata(hax, 'dynamicMarkerSizeListener', hListener);
  catch    
    if ~isappdata(hax, 'dynamicMarkerSizeListener')
      hListener = addlistener(hax,'Position','PostSet',...
        @(obj,events) localResizeScatterCallback(obj,events,sP(i).ax));
%      localResizeScatterCallback([],[],sP(i).ax);
      setappdata(hax, 'dynamicMarkerSizeListener', hListener);
    end
    %disp('some Error!');
  end

  % plot labels
  if check_option(varargin,{'text','label','labeled'})
    text(v,get_option(varargin,{'text','label'}),'parent',sP(i).ax,'addMarkerSpacing',varargin{:});
  end

  if isappdata(sP(1).parent,'mtexFig')
    mtexFig = getappdata(sP(1).parent,'mtexFig');
    mtexFig.drawNow('figSize',getMTEXpref('figSize'),varargin{:});
  end
end

if nargout == 0
  clear h;
else
  ax = [sP.ax];
end

end


% ---------------------------------------------------------------
function localResizeScatterCallback(h,e,hax)
% get(fig,'position')

hax = handle(hax);

% get markerSize
markerSize = get(findobj(hax,'type','patch'),'MarkerSize');
if isempty(markerSize)
  markerSize = 0;
elseif iscell(markerSize)
  markerSize = [markerSize{:}];
end

markerSize = max(markerSize);

% correct text positions
t = findobj(hax,'Tag','setBelowMarker');
correctTextPostion(t,markerSize,-1);

t = findobj(hax,'Tag','setAboveMarker');
correctTextPostion(t,markerSize,1);

% ------------- scale scatterplots -------------------------------
u = findobj(hax,'Tag','dynamicMarkerSize');

if isempty(u), return;end

p = get(u(1),'parent');
while ~isgraphics(p,'axes'), p = get(p,'parent'); end


unit = get(p,'unit');
set(p,'unit','pixel')
pos = get(p,'position');
l = min([pos(3),pos(4)]);
if l < 0, return; end 

for i = 1:length(u)
  d = get(u(i),'UserData');
  o = get(u(i),'MarkerSize');
  %n = l/350 * d;
  n = l/250 * d;
  if abs((o-n)/o) > 0.05, set(u(i),'MarkerSize',n);end
  
end

set(p,'unit',unit);

end

function correctTextPostion(t,markerSize,direction)
% adjust text positions

for it = 1:length(t)
  
  xy = get(t(it),'UserData');
  if any(isnan(xy)), continue; end
  set(t(it),'unit','data','position',[xy,0]);
  set(t(it),'unit','pixels');
  xy = get(t(it),'position');
  if isappdata(t(it),'extent')
    extend = getappdata(t(it),'extent');
  else
    extend = get(t(it),'extent');
    setappdata(t(it),'extent',extend);
  end
  margin = get(t(it),'margin');
  xy(2) = xy(2) + direction*(extend(4)/2 + margin + markerSize/2 + 5);
    
  set(t(it),'position',xy);
  set(t(it),'unit','data');
end

end
