function [Transmission rou]=ReadTransData(atomparam_atomno)

if atomparam_atomno==13
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
elseif atomparam_atomno==26
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
else
    sprintf('Warning: you are using transmission data of Al by default. You may want to add transmission data for your own sample!')
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
end

