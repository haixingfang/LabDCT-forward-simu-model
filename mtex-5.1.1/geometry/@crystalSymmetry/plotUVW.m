function plotUVW(cs,varargin)
% plot symmetry
%
% Input
%  s - symmetry
%
% Output
%
% Options
%  antipodal      - include [[AxialDirectional.html,antipodal symmetry]]

mtexFig = newMtexFigure(varargin{:});

% which directions to plot
m = Miller({1,0,0},{0,1,0},{0,0,1},{1,1,0},{0,1,1},{1,0,1},{1,1,1},cs,'uvw');

m = unique(m);
options = [{'symmetrised','labeled','MarkerEdgeColor','k','grid','doNotDraw'},varargin];
  
% plot them
washold = getHoldState(mtexFig.gca);
hold(mtexFig.gca,'all')
for i = 1:length(m)
  m(i).scatter(options{:});
end
hold(mtexFig.gca,washold)


% postprocess figure
setappdata(gcf,'CS',cs);
set(gcf,'tag','ipdf');
mtexFig.drawNow('figSize',getMTEXpref('figSize'),varargin{:});
    
