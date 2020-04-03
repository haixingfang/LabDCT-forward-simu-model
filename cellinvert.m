%
% cellinvert returns unit cell in reciprocal space 
%
%     recipr_unit_cell = cellinvert(unit_cell)
%
% INPUT: unit_cell = [a b c alpha beta gamma]
% OUPUT: recipr_unit_cell = [a* b* c* alpha* beta* gamma*]

function recipcell = cellinvert(ucell)

     a   = ucell(1);
     b   = ucell(2);
     c   = ucell(3);
     calp = cos(ucell(4)*pi/180.);
     cbet = cos(ucell(5)*pi/180.);
     cgam = cos(ucell(6)*pi/180.);
     salp = sin(ucell(4)*pi/180.);
     sbet = sin(ucell(5)*pi/180.);
     sgam = sin(ucell(6)*pi/180.);
 
  angular = sqrt(1 - calp^2 - cbet^2 - cgam^2 + 2*calp*cbet*cgam);
  V = a*b*c*angular;                                           %Volume of unit cell
 
 
  %  Calculate reciprocal lattice parameters: 
  astar = b*c*salp/V                          ;
  bstar = a*c*sbet/V                          ;
  cstar = a*b*sgam/V                          ;
  salpstar = V/(a*b*c*sbet*sgam)              ;
  sbetstar = V/(a*b*c*salp*sgam)              ;
  sgamstar = V/(a*b*c*salp*sbet)              ;
  calpstar = (cbet*cgam-calp)/(sbet*sgam)     ;
  cbetstar = (calp*cgam-cbet)/(salp*sgam)     ;
  cgamstar = (calp*cbet-cgam)/(salp*sbet)     ;

  alpstar = acos(calpstar)*180/pi;
  betstar = acos(cbetstar)*180/pi;
  gamstar = acos(cgamstar)*180/pi;
  
  recipcell = [astar bstar cstar alpstar betstar gamstar];