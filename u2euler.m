% Euler angles (phi1, PHI, phi2) from U matrix
% The formalism follows the ID11-3DXRD specs
% Note that there are two solutions
% (phi1, PHI, phi2) AND (phi1 + pi, -PHI, phi2 + pi)
% We pick the one with phi1 in the range [-pi/2 pi/2]
%
% Fails if U(3,2) or U(2,3) = 0 e.g. then U(3,3) = ~1
% If U(3,3) ~ 1 ph1 = ph2 = atan(U(2,1)/U(1,1))/2
% In this case there is only one solution.

function phi_vector = u2euler(U)

ph(1)=acos(U(3,3));
if ph(1) < 0.0001
    ph2(1) = atan(U(2,1)/U(1,1))/2.0
    ph1(1) = ph2(1)
    phi1_deg = ph1(1)*180/pi;
    phi_deg = ph(1)*180/pi;
    phi2_deg = ph2(1)*180/pi;
else
    % there is 2 solutions to each of the equations
    ph1(1)=atan(-U(1,3)/U(2,3));
    ph2(1)=atan(U(3,1)/U(3,2));
    ph(2)=2*pi-ph(1);
    ph1(2)=ph1(1)+pi;
    ph2(2)=ph2(1)+pi;

    %the right combination is found by brute-force, using the restriction on phi1
    minsum = 1000;
    for j=1:2
        for k=1:2
            U2 = euler2u(ph1(2)*180/pi,ph(j)*180/pi,ph2(k)*180/pi);
            Udev = abs(U2-U);
            sum1 = sum(sum(Udev));
            if sum1<minsum
                minsum=sum1;
                mj = j;
                mk = k;
            end
        end
    end

    phi1_deg = ph1(2)*180/pi;
    phi_deg = ph(mj)*180/pi;
    % if ph2(mk) < 0
    %     ph2(mk) =ph2(mk)+2*pi;
    % end
    phi2_deg = ph2(mk)*180/pi;
end
phi_vector = [phi1_deg phi_deg phi2_deg];



