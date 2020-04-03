% generate grain colors according to grain orientation
% IPF along Z direction
function GrainColor=gen_graincolor(euler_grains)

%%% change the default settings
setMTEXpref('xAxisDirection','east'); % default: 'north'
setMTEXpref('zAxisDirection','outOfPlane'); % same as default
setMTEXpref('bAxisDirection','north'); % default: 'east'
setMTEXpref('aAxisDirection',''); % same as default
setMTEXpref('FontSize',44); % default: 15

cs = crystalSymmetry('cubic');
% ss = specimenSymmetry('orthorhombic');

rot = rotation('Euler',euler_grains(:,1)*degree,euler_grains(:,2)*degree,euler_grains(:,3)*degree);
o = orientation(rot,cs,ss);
ebsd_o=o;
cK = ipfHSVKey(cs);
cK.inversePoleFigureDirection=zvector;
GrainColor = cK.orientation2color(ebsd_o);

