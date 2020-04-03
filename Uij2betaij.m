% Uij2betaij transform the ADP U-matrix into the beta form 
%
% betaij = Uij2betaij(adp,unit_cell)
%
% INPUT:  adp: anisotropic displacement parameter U matrix
%         unit_cell = [a b c alpha beta gamma] 
% OUTPUT: betaij: beta displacement matrix

function betaij = Uij2betaij(adp,cell)

U  = [adp(1) adp(6) adp(5);
      adp(6) adp(2) adp(4); 
      adp(5) adp(4) adp(3)];

% U = 2*pi^2*[U11 U12 U13; U12 U22 U23; U13 U23 U33]
 betaij =zeros(3,3);
 cellstar = cellinvert(cell);
 
for i = 1:3
    for j = 1:3
        betaij(i,j) = 2*pi^2*cellstar(i)*cellstar(j)*U(i,j);
    end
end


