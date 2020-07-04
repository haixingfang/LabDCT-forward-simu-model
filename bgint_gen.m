

% This function generates background intensity according to the mass
% attenuation
% Applies to only Al and Fe at present
% For other elements, they are regarded as equivalent to Fe

if atomparam.atomno==13
    rou=2.70; % density [g/cm^3]
    fileID=fopen('Al_transmission.txt','r');
    Transmission=[];
    while(~feof(fileID))
        textdata=str2num(fgetl(fileID));
        if isempty(textdata)
            continue;
        else
            Transmission=[Transmission; textdata];
        end
    end
    fclose(fileID);
    for i=1:length(Energy)
        for j=1:length(Transmission)-1
            if Transmission(j,1)<=Energy(i)/1e3 && Transmission(j+1,1)>=Energy(i)/1e3
                Eindex=[j j+1];
            end
        end       
        if ~exist('Lx','var')
            Lx=400; % [um]
        end
        IE(i)=I0E(i)*exp(-(Transmission(Eindex(1),2)+(Energy(i)/1e3-(Transmission(Eindex(1),1))) ...
            *(Transmission(Eindex(2),2)-Transmission(Eindex(1),2))/ ...
            (Transmission(Eindex(2),1)-Transmission(Eindex(1),1)))*Lx*1e-4*rou);
    end
end

if atomparam.atomno==26 || atomparam.atomno~=13
    rou=7.87; % density [g/cm^3]
    fileID=fopen('Fe_transmission.txt','r');
    Transmission=[];
    while(~feof(fileID))
        textdata=str2num(fgetl(fileID));
        if isempty(textdata)
            continue;
        else
            Transmission=[Transmission; textdata];
        end
    end
    fclose(fileID);
    for i=1:length(Energy)
        for j=1:length(Transmission)-1
            if Transmission(j,1)<=Energy(i)/1e3 && Transmission(j+1,1)>=Energy(i)/1e3
                Eindex=[j j+1];
            end
        end
        if ~exist('Lx','var')
            Lx=400; % [um]
        end
        IE(i)=I0E(i)*exp(-(Transmission(Eindex(1),2)+(Energy(i)/1e3-(Transmission(Eindex(1),1))) ...
            *(Transmission(Eindex(2),2)-Transmission(Eindex(1),2))/ ...
            (Transmission(Eindex(2),1)-Transmission(Eindex(1),1)))*Lx*1e-4*rou);
    end 
end
bgPeak=sum(IE.*ExpTime); % intensity [photons]

% fit pseudo Voigt function for 2D intensity profile
ygrid=-detysize/2+1:detysize/2;
zgrid=-detzsize/2+1:detzsize/2;
[Ygrid,Zgrid]=meshgrid(ygrid,zgrid);
Ln=1./(1+(Ygrid./(2*Lx/pixelysize/1000)).^2.*(Zgrid./(2*Lx/pixelzsize/1000)).^2);
Gn=exp(-0.5*(Ygrid./(2*Lx/pixelysize/1000)).^2-0.5*(Zgrid./(2*Lx/pixelzsize/1000)).^2);
mun=0.1;
Ibg=bgPeak.*(mun*Ln+(1-mun)*Gn);
% figure;
% imagesc(ygrid,zgrid,Ibg);
% colormap(jet);
bgint=Ibg(fix(length(Ygrid)/4),fix(length(Ygrid)/4))./Ibg(fix(length(Ygrid)/2),fix(length(Ygrid)/2))*(2^16-1);% adopted before Nov 21, 2019
if SampleCylinderFlag==0
    bgint=918.9030; % set as constant, June 22, 2020
end
% bgint=(sum(sum(Ibg(fix(length(Ygrid)/4):fix(length(Ygrid)/3),fix(length(Ygrid)/4):fix(length(Ygrid)/3))))./ ...
%     (fix(length(Ygrid)/3)-fix(length(Ygrid)/4))^2)./ ...
%     Ibg(fix(length(Ygrid)/2),fix(length(Ygrid)/2))*(2^16-1); % modified on Nov 21, 2019
% bgint=(2^12-1)/2;





