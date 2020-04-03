%
% Parameters describing experimental setup 
%
% Created by Haixing Fang, March 2019
% hfang@mek.dtu.dk, haixingfang868@gmail.com
%
clear all;
%%% load dipimage, mpt3 and mtex toolbox
load_diplib;
load_mpt3; % see documentation, type 'mptdoc'
load_mtex;

% input crystallography data: input_al, input_fe, input_si, input_ni are available
% other crystallography data can be edited
input_al;
% input_fe;

% Experimental setup

mono_beam=0; % monochoromatic photon beam or not, 0-polychromatic; 1-monochromatic. By default it is polychromatic
% I0 = 1e14; % Beam flux from a synchrotron source 10^14-15 (photons/s/mm2)
I0 = 5e8;    % Beam flux from a lab source 10^8 (photons/s/mm2)
ExpTime=600; % exposure time [s]
if mono_beam==1
    Energy = 59.31;  % photon energy [keV], K alpha1 line of W anode material
end
if mono_beam==0
    [Energy I0E]=Xray_Spectrum_W(I0);
    I0E=abs(I0E);
end
% print('X-ray Spectrum','-dtiff');

% detector of Zeiss Xradia 520
exp_parameters_input;

tilt_x   = 0;         % detector tilt counterclockwise around lab x axis [rad] 
tilt_y   = 0;         % detector tilt counterclockwise around lab y axis [rad] 
tilt_z   = 0;         % detector tilt counterclockwise around lab z axis [rad] 
polfactor = 1;   % Apply polarisation factor
Lorentz = 1;     % Apply Lorentz factor
beampol = 1.0;   % Beam polarisation (1 = totally horizontally polarised)
                 %                   (0 = nonpolarised beam)
% Maximum theta angle for reflection generation
tthetamax = acos(dot([Lsam2det 0 0]./norm([Lsam2det 0 0]), ...
    [Lsam2det 0.5*detysize*pixelysize 0.5*detzsize*pixelzsize]./norm([Lsam2det 0.5*detysize*pixelysize 0.5*detzsize*pixelzsize])));
tthetamax = tthetamax*180/pi; % [deg]


% grain size distribution
grain_file=0; % 0: no grain file; 1: read grain information from external file such as .h file exported by GrainMapper
grain_flag=1; % 0: manually generate; 1: call gen_micro_cylinder to generate grain structure
mesh_flag=1; % 0: no 3D mesh on individual grain; 1: 3D mesh on individual grain
grainsize_select=0; % 0: select all grains; 1: only select grains within the defined range
if mesh_flag==0
    voxel_flag=1; % 0: not to use voxels to mesh the grain; 1: use voxels to mesh the grain
else
    voxel_flag=0;
end

% generagte virtual microstructure
if grain_file==0
    % cylinder structure
    % recommend sample volumes: 100*100*150, 200*200*300 or 400*400*600
    Lx=200; % cylinder diameter [um]
    Lz=300; % cylinder height [um]
    gr_diameter=60; % average grain diameter [um]
    [A_P V_microstructure]=gen_micro_cylinder(Lx,Lz,gr_diameter);
    
%     % gradient structre: from top to bottom, grains are small to large
%     Lx=100; % cylinder diameter [um]
%     Lz=150; % cylinder height [um]
%     gr_diameter=2; % average grain diameter of the top layer [um]
%     [A_P V_microstructure]=gen_micro_cylinder_gradient(Lx,Lz,gr_diameter);    
    
    grains = length(A_P(:,1)); % number of grains
    sampos = A_P(:,2:4)./1000; % coordinates of grain centroids [mm]
    euler_grains=A_P(:,10:12); % Euler angles [deg]: ([0 2*pi], [0 pi], [0 2*pi])
    grainsize=A_P(:,13); % equivalent diameter of grain [um]
    grainvolume=A_P(:,5); % [um^3]
    
    figure('Name','First look at size distribution');
    xbins=min(grainsize):5:max(grainsize);
    [counts,centers] = hist(grainsize,xbins);
    BIN_WIDTH = centers(2)-centers(1);
    BIN_MAX = max(xbins)+1/2*BIN_WIDTH;
    bar(centers,counts./(length(grainsize)*BIN_WIDTH));
    hold on;
%     [parmhat parmci] = lognfit(grainsize);
%     y_probability = BIN_WIDTH*lognpdf(centers,parmhat(1),parmhat(2));
%     y_count = length(grainsize) * y_probability;
%     % hfit = plot(centers,y_probability,'o');
%     % set(hfit,'LineWidth',2);
%     cs = csapi(centers,y_probability);
%     fnplt(cs,2,'r');
%     E=exp(parmhat(1)+parmhat(2)^2/2); % corrected E
%     SD=exp(parmhat(1)+1/2*parmhat(2)^2)*sqrt((exp(parmhat(2)^2)-1)); % corrected SD
    xlabel('d_{grain} (\mum)');
    ylabel('Probability density (\mum^{-1})');
    xlim([min(grainsize)- BIN_WIDTH/2 max(grainsize)+ BIN_WIDTH/2])
    ylim([0 ceil(max(100*counts./(length(grainsize)*BIN_WIDTH)))/100]);
    set(gca,'fontsize',16);
    set(gca,'linewidth',1.5);
    box on;hold off;
    if mesh_flag==1
       MeshGrain=length(A_P(:,1)); % define how many grains you want to mesh for simulations
       for jj=1:MeshGrain
           SingleGrain=Polyhedron(V_microstructure.Set(jj));
           MeshNr(jj)=round(SingleGrain.volume^(1/3)/mean([pixelysize pixelzsize].*1000));
           % Adjust MeshNr for tuning the polyheron size
           % preferably to make sure polyheron size smaller than 12.5 um (at least smaller than 25 um)
           if MeshNr(jj)>12
               MeshNr(jj)=12;
           elseif MeshNr(jj)<2
               MeshNr(jj)=2;
           end
           [Vsub{jj} SubGrain{jj}]=mesh3D(SingleGrain,MeshNr(jj)); % cell number,x,y,z,volume,EqDiameter,EqDiameter of the parent grain
           SubGrain{jj}(:,2:4)=SubGrain{jj}(:,2:4)./1000; % [x,y,z] [mm]
           jj
       end
       GrainColor=gen_graincolor(euler_grains); % get IPF color
       figure;
       hold all;
       for jj=1:MeshGrain
%            ColorGrain=rand(1,3);
           ColorGrain=GrainColor(jj,:);
           V_microstructure.Set(jj).plot('color',ColorGrain,'alpha',0.5,'edgecolor','k','linestyle','none','linewidth',0.75);
%            for kk=1:length(Vsub{jj}.Set) % plot the mesh
%                Vsub{jj}.Set(kk).plot('color',ColorGrain,'alpha',0.15,'edgecolor','k','linewidth',0.25); % Adjust the transparancy
% %                Vsub{jj}.Set(kk).plot('alpha',0.15,'edgecolor','k','linewidth',0.75); % Adjust the transparancy
% %                Vsub{jj}.Set(kk).plot('alpha',0.15); % Adjust the transparancy
%            end
       end
       xlabel('x (\mum)','FontSize',16);
       ylabel('y (\mum)','FontSize',16);
       zlabel('z (\mum)','FontSize',16);
       axis equal;
       axis([-Lx/2 Lx/2 -Lx/2 Lx/2 -Lz/2 Lz/2]);
       set(get(gca,'xlabel'),'rotation',30);
       set(get(gca,'ylabel'),'rotation',-30);
       set(get(gca,'zlabel'),'rotation',90);
       view(3);
       set(gca,'fontsize',16);
       set(gca,'linewidth',1.5);       
       box on;
    end
    if ~isempty(MeshGrain) && MeshGrain>0
        grains = MeshGrain; % number of grains
        sampos = A_P(1:MeshGrain,2:4)./1000; % coordinates of grain centroids [mm]
        euler_grains=A_P(1:MeshGrain,10:12); % Euler angles [deg]: ([0 2*pi], [0 pi], [0 2*pi])
        grainsize=A_P(1:MeshGrain,13); % equivalent diameter of grain [um]
        grainvolume=A_P(1:MeshGrain,5); % [um^3]
    end
end

% read grain information from .h5 file
% .h5 is a grain structure output file from GrainMapper3D
% It can be also customized to write a .h5 file for grain structure obtained from other techniques such as 3DXRD, DAXM etc.
if grain_file==1
    FileFolder=strcat(pwd,'\Examples');
    h5FileName='3Dvolume_8grains_2_5um.h5';
    Prefix_Name=h5FileName(1:end-3);

    DS = readLabDCT(fullfile(FileFolder,h5FileName)); % read h5 file that is exported by LabDCT [X*Y*Z]
    
    % 1: for customized grain selection
    % 0: for selecting all the grains, default
    if grainsize_select==1
        grainsize_range=[40 200]; % [um]
        graincomp_range=[0.65 1]; % completeness
        grainsize_vox=pi/6*grainsize_range.^3*1e-9./(DS.VoxSize(1)*DS.VoxSize(2)*DS.VoxSize(3)); % [voxels]
        eff_index=find(DS.nVox>=grainsize_vox(1) & DS.nVox<=grainsize_vox(2));
        eff_index=eff_index(find(DS.SeedComp(eff_index)>=graincomp_range(1) & DS.SeedComp(eff_index)<=graincomp_range(2)));
        DS.SeedID=DS.SeedID(eff_index);
        DS.SeedComp=DS.SeedComp(eff_index);
        DS.RodVec=DS.RodVec(eff_index,:);
        DS.nVox=DS.nVox(eff_index,:);
        DS.Coord=DS.Coord(eff_index,:);
        Dream3D_FileName=[Prefix_Name '_D' num2str(grainsize_range(1)) '_' num2str(grainsize_range(2)) '.dream3d']; % create a new name for dream3D
        dfile=fullfile(FileFolder,Dream3D_FileName); % file path for creating the dream3D file
        if exist(dfile,'file')==0
            Dream3DWriter(DS,dfile);
        end
        Xdmf_FileName=[Prefix_Name '_D' num2str(grainsize_range(1)) '_' num2str(grainsize_range(2)) '.xdmf']; % create a new name for xdmf file
        LabDCT_XdmfWriter(FileFolder,Dream3D_FileName,Xdmf_FileName);
    else 
        Dream3D_FileName=[Prefix_Name '.dream3d']; % create a new name for dream3D
        dfile=fullfile(FileFolder,Dream3D_FileName); % file path for creating the dream3D file
        if exist(dfile,'file')==0
            Dream3DWriter(DS,dfile);
        end
        Xdmf_FileName=[Prefix_Name '.xdmf']; % create a new name for xdmf file
        if exist(fullfile(FileFolder,Xdmf_FileName),'file')==0
            LabDCT_XdmfWriter(FileFolder,Dream3D_FileName,Xdmf_FileName);
        end
    end
    
    % rotation axis of the DCT scan
    % note:in the lab coordinate system, z-along the beam, y-sample height, x-perpendicular to the beam
    % referenced to the DCT setup at position of 0 degree
    % consider the center shift from the reconstruction
    center_shift(1)= DS.Center(1); % x: along the beam [mm]   (lab-z)
    center_shift(2)= DS.Center(2); % y: perpendicular to the beam [mm] (lab-x)
    center_shift(3)= DS.Center(3); % z: sample height direction, irrelevant [mm] (lab-y)
    
    % In Dipimage, xyz has been denoted as yxz
    if (all(DS.Mask,'all')==0)
        DS.Mask=ones(size(DS.Mask));
    end
    grain_mask=dip_image(DS.Mask,'uint8'); % grain labeled image
    grain_label=dip_image(DS.GIDvol,'int32'); % grain labeled image
    grain_measure = measure(grain_label,grain_mask,{'dimension','DimensionsCube','DimensionsEllipsoid','gravity', ...
                        'MajorAxes','Inertia','P2A','Size','SurfaceArea','Inertia'});
    
    grains = length(DS.SeedID); % number of grains
    sampos(:,1) = (DS.Coord(:,1)-DS.Dimension(1)/2).*DS.VoxSize(1) ...
        +center_shift(1); % coordinates of grain centroids [mm]
    sampos(:,2) = (DS.Coord(:,2)-DS.Dimension(2)/2).*DS.VoxSize(2) ...
        +center_shift(2); % coordinates of grain centroids [mm]
    sampos(:,3) = (DS.Coord(:,3)-DS.Dimension(3)/2).*DS.VoxSize(3)+center_shift(3); % coordinates of grain centroids [mm]
    euler = Euler(rodrigues2quat(vector3d(double(DS.RodVec'))));
    euler_grains=euler.*180/pi;% Euler angles [deg]: ([0 2*pi], [0 pi], [0 2*pi])
    grainvolume=DS.nVox*DS.VoxSize(1)*DS.VoxSize(2)*DS.VoxSize(3)*1e9; % [um^3]
    grainsize=2*(3*grainvolume/(4*pi)).^(1/3); % equivalent diameter of grain [um]
    if mesh_flag==1
        for jj=1:length(DS.SeedID)%length(grain_measure.ID)
            if grainsize_select==0
                grain_selection=dip_image(grain_label==jj,'bin');
            else
                grain_selection=dip_image(grain_label==DS.SeedID(jj),'bin');
            end
            vertices=[];
            
            % replaced by a new method shown below after Feb 19, 2020
            %{
            BoundBox(1,1)=round(grain_measure(DS.SeedID(jj)).Gravity(1)-grain_measure(DS.SeedID(jj)).CartesianBox(1)/2-1);
            BoundBox(1,2)=round(grain_measure(DS.SeedID(jj)).Gravity(1)+grain_measure(DS.SeedID(jj)).CartesianBox(1)/2+1);
            BoundBox(2,1)=round(grain_measure(DS.SeedID(jj)).Gravity(2)-grain_measure(DS.SeedID(jj)).CartesianBox(2)/2-1);
            BoundBox(2,2)=round(grain_measure(DS.SeedID(jj)).Gravity(2)+grain_measure(DS.SeedID(jj)).CartesianBox(2)/2+1);
            BoundBox(3,1)=round(grain_measure(DS.SeedID(jj)).Gravity(3)-grain_measure(DS.SeedID(jj)).CartesianBox(3)/2-1);
            BoundBox(3,2)=round(grain_measure(DS.SeedID(jj)).Gravity(3)+grain_measure(DS.SeedID(jj)).CartesianBox(3)/2+1);
            if BoundBox(1,1)<0
                BoundBox(1,1)=0;
            end
            if BoundBox(1,2)>length(grain_selection(:,1,1))-2
                BoundBox(1,2)=length(grain_selection(:,1,1))-2;
            end
            if BoundBox(2,1)<0
                BoundBox(2,1)=0;
            end
            if BoundBox(2,2)>length(grain_selection(1,:,1))-2
                BoundBox(2,2)=length(grain_selection(1,:,1))-2;
            end
            if BoundBox(3,1)<0
                BoundBox(3,1)=0;
            end
            if BoundBox(3,2)>length(grain_selection(1,1,:))-2
                BoundBox(3,2)=length(grain_selection(1,1,:))-2;
            end
            
%             % using the distance map approach is faster
%             [D,IDX] = bwdist(double(grain_selection));
%             IDX_unique=unique(IDX);
%             for j0=1:length(IDX_unique)
%                 [vertices0(j0,1) vertices0(j0,2) vertices0(j0,3)]=ind2sub(size(D),IDX_unique(j0));% [y x z]
%             end
%             for j0=min(vertices0(:,3)):max(vertices0(:,3))
%                 [Bgb,Lgb] = bwboundaries(double(grain_selection(:,:,j0-1)),'noholes');
%             end
            
            % step size to look for vertices [pixel]
            pixel_step=22;
            stepsize(1)=round(grain_measure.CartesianBox(1,jj)/pixel_step);
            if stepsize(1)<=0
                stepsize(1)=1;
            end
            if stepsize(1)>=3
                stepsize(1)=3;
            end
            stepsize(2)=round(grain_measure.CartesianBox(2,jj)/pixel_step);
            if stepsize(2)<=0
                stepsize(2)=1;
            end
            if stepsize(2)>=3
                stepsize(2)=3;
            end
            stepsize(3)=round(grain_measure.CartesianBox(3,jj)/pixel_step);
            if stepsize(3)<=0
                stepsize(3)=1;
            end
            if stepsize(3)>=3
                stepsize(3)=3;
            end
%             for i=0:stepsize(1):length(grain_selection(:,1,1))-2 % dipimage object starts index from 0
%                 for j=0:stepsize(2):length(grain_selection(1,:,1))-2
%                     for k=0:stepsize(3):length(grain_selection(1,1,:))-2
% %                         if (grain_selection(i,j,k)==0 && grain_selection(i,j,k+1)==1) ...
% %                                 || (grain_selection(i,j,k)==0 && grain_selection(i,j+1,k)==1) ...
% %                                 || (grain_selection(i,j,k)==0 && grain_selection(i+1,j,k+1)==1)
%                         if (grain_selection(i,j,k)==0 && grain_selection(i,j,k+1)==1) ...
%                                 || (grain_selection(i,j,k)==0 && grain_selection(i,j+1,k)==1) ...
%                                 || (grain_selection(i,j,k)==0 && grain_selection(i+1,j,k)==1) ...
%                                 || (grain_selection(i,j,k)==0 && grain_selection(i,j+1,k+1)==1) ...
%                                 || (grain_selection(i,j,k)==0 && grain_selection(i+1,j+1,k)==1) ...
%                                 || (grain_selection(i,j,k)==0 && grain_selection(i+1,j,k+1)==1)
%                             vertices=[vertices;i j k];
%                         end
%                     end
%                 end
% %                 i
%             end
            % this is more efficient, Oct 30, 2019
            for i=BoundBox(1,1):stepsize(1):BoundBox(1,2) % dipimage object starts index from 0
                for j=BoundBox(2,1):stepsize(2):BoundBox(2,2)
                    for k=BoundBox(3,1):stepsize(3):BoundBox(3,2)
%                         if (grain_selection(i,j,k)==0 && grain_selection(i,j,k+1)==1) ...
%                                 || (grain_selection(i,j,k)==0 && grain_selection(i,j+1,k)==1) ...
%                                 || (grain_selection(i,j,k)==0 && grain_selection(i+1,j,k+1)==1)
                        if ((grain_selection(i,j,k)==0 && grain_selection(i,j,k+1)==1) ...
                                || (grain_selection(i,j,k)==0 && grain_selection(i,j+1,k)==1) ...
                                || (grain_selection(i,j,k)==0 && grain_selection(i+1,j,k)==1) ...
                                || (grain_selection(i,j,k)==0 && grain_selection(i,j+1,k+1)==1) ...
                                || (grain_selection(i,j,k)==0 && grain_selection(i+1,j+1,k)==1) ...
                                || (grain_selection(i,j,k)==0 && grain_selection(i+1,j,k+1)==1) ...
                                || (grain_selection(i,j,k)==1 && grain_selection(i,j,k+1)==0) ...
                                || (grain_selection(i,j,k)==1 && grain_selection(i,j+1,k)==0) ...
                                || (grain_selection(i,j,k)==1 && grain_selection(i+1,j,k)==0) ...
                                || (grain_selection(i,j,k)==1 && grain_selection(i,j+1,k+1)==0) ...
                                || (grain_selection(i,j,k)==1 && grain_selection(i+1,j+1,k)==0) ...
                                || (grain_selection(i,j,k)==1 && grain_selection(i+1,j,k+1)==0))
                            vertices=[vertices;i j k];
                        end
                    end
                end
%                 i
            end
            grain_vertices{jj}(:,1)=((vertices(:,2)+1)-DS.Dimension(1)/2).*DS.VoxSize(1)+center_shift(1); % x [mm]
            grain_vertices{jj}(:,2)=((vertices(:,1)+1)-DS.Dimension(2)/2).*DS.VoxSize(2)+center_shift(2); % y [mm]
            grain_vertices{jj}(:,3)=((vertices(:,3)+1)-DS.Dimension(3)/2).*DS.VoxSize(3)+center_shift(3); % z [mm]
            %}
            
            % new method updated on Feb 19, 2020
            im1 = countneighbours(grain_selection,1,Inf,0);
%             im2=im1>=26;
            im2=im1>=21;
            im3=grain_selection-im2;
            im4=double(im3);
            clear vertices;
            [vertices(:,1),vertices(:,2),vertices(:,3)]=ind2sub(size(im4),find(im4==1)); % indices
            grain_vertices{jj}(:,1)=((vertices(:,1))-DS.Dimension(1)/2).*DS.VoxSize(1)+center_shift(1); % x [mm]
            grain_vertices{jj}(:,2)=((vertices(:,2))-DS.Dimension(2)/2).*DS.VoxSize(2)+center_shift(2); % y [mm]
            grain_vertices{jj}(:,3)=((vertices(:,3))-DS.Dimension(3)/2).*DS.VoxSize(3)+center_shift(3); % z [mm]
            
            SingleGrain=Polyhedron(grain_vertices{jj});
            P_grain(jj)=SingleGrain;
           % Adjust MeshNr for tuning the polyheron size
           % preferably to make sure polyheron size smaller than 12.5 um (at least smaller than 25 um)
%             MeshNr(jj)=round((SingleGrain.volume.*10^9)^(1/3)/mean([pixelysize pixelzsize].*1000));
            MeshNr(jj)=round(log(3*20/((SingleGrain.volume.*10^9)^(1/3)))/(-0.1));
            if MeshNr(jj)>15
                MeshNr(jj)=15;%9
            elseif MeshNr(jj)<2
                MeshNr(jj)=2;
            end
            [Vsub{jj} SubGrain{jj}]=mesh3D(SingleGrain,MeshNr(jj)); % cell number,x,y,z,volume,EqDiameter,EqDiameter of the parent grain, unit based on mm
            cell2Vol(jj)=sum(SubGrain{jj}(:,5))./(DS.nVox(jj)*DS.VoxSize(1)*DS.VoxSize(2)*DS.VoxSize(3)); % total volume
            jj
        end
        GrainColor=gen_graincolor(euler_grains); % get IPF color
        figure;
        subplot(1,2,1);
        grain_union=Union(P_grain);
        grain_union.plot('alpha',0.55);
        xlim([-DS.Dimension(1)/2 DS.Dimension(1)/2].*DS.VoxSize(1)+center_shift(1));
        ylim([-DS.Dimension(2)/2 DS.Dimension(2)/2].*DS.VoxSize(2)+center_shift(2));
        zlim([-DS.Dimension(3)/2 DS.Dimension(3)/2].*DS.VoxSize(3)+center_shift(3));
        xlabel('x (mm)','FontSize',16);
        ylabel('y (mm)','FontSize',16);
        zlabel('z (mm)','FontSize',16);
        set(get(gca,'xlabel'),'rotation',15);
        set(get(gca,'ylabel'),'rotation',-25);
        set(get(gca,'zlabel'),'rotation',90);
        set(gca,'fontsize',16);
%         set(gca,'linewidth',1.5);
        hold off
        box on;
        title('(a) Grain geometry');
        %{ % plot the mesh can be very slow
        subplot(1,2,2);
        hold all;       
        for jj=1:length(grain_measure(DS.SeedID))
           ColorGrain=GrainColor(jj,:);
           for kk=1:length(Vsub{jj}.Set)
%                Vsub{jj}.Set(kk).plot('FaceColor','cyan','alpha',0.15,'edgecolor','k','linewidth',0.75); % Adjust the transparancy
%         %                Vsub{jj}.Set(kk).plot('alpha',0.15,'edgecolor','k','linewidth',0.75); % Adjust the transparancy
%         %                Vsub{jj}.Set(kk).plot('alpha',0.15); % Adjust the transparancy
               Vsub{jj}.Set(kk).plot('color',ColorGrain,'alpha',0.15,'edgecolor','k','linewidth',0.75); % Adjust the transparancy
           end
        end
        Vsub{jj}.plot('FaceColor','cyan','alpha',0.15,'edgecolor','k','linewidth',0.75);
        xlim([-DS.Dimension(1)/2 DS.Dimension(1)/2].*DS.VoxSize(1)+center_shift(1));
        ylim([-DS.Dimension(2)/2 DS.Dimension(2)/2].*DS.VoxSize(2)+center_shift(2));
        zlim([-DS.Dimension(3)/2 DS.Dimension(3)/2].*DS.VoxSize(3)+center_shift(3));
        xlabel('x (mm)','FontSize',16);
        ylabel('y (mm)','FontSize',16);
        zlabel('z (mm)','FontSize',16);
        set(get(gca,'xlabel'),'rotation',15);
        set(get(gca,'ylabel'),'rotation',-25);
        set(get(gca,'zlabel'),'rotation',90);
        set(gca,'fontsize',16);
%         set(gca,'linewidth',1.5);
        hold off;
        box on;
        title('(b) 3D mesh');
        %}
    elseif voxel_flag==1
%         xx=dip_image(DS.GIDvol,'int8');
        xx=dip_image(DS.GIDvol(:,:,1:end-1),'int8');
        Binning=[4 4 3]; % define the binning size in xyz directions, make sure dimensions./Binning = integers
        for jj=1:length(grain_measure.ID)
            yy=xx==jj;
            zz=rebin(yy,Binning);
            k=0;
            for p=0:length(zz(:,1,1))-1
                for q=0:length(zz(1,:,1))-1
                    for r=0:length(zz(1,1,:))-1
                        if zz(p,q,r)==1
                            k=k+1;
                            SubGrain{jj}(k,1)=k;
                            SubGrain{jj}(k,3)=(p-(length(zz(:,1,1))-1)/2).*DS.VoxSize(1)*Binning(1)+center_shift(2); % centroid coordinate y
                            SubGrain{jj}(k,2)=(q-(length(zz(1,:,1))-1)/2).*DS.VoxSize(2)*Binning(2)+center_shift(1); % centroid coordinate x
                            SubGrain{jj}(k,4)=(r-(length(zz(1,1,:))-1)/2).*DS.VoxSize(3)*Binning(3)+center_shift(3); % centroid coordinate z
                            SubGrain{jj}(k,5)=DS.VoxSize(1)*DS.VoxSize(2)*DS.VoxSize(3)*Binning(1)*Binning(2)*Binning(3); % volume
                            SubGrain{jj}(k,6)=2*(3*SubGrain{jj}(k,5)/(4*pi))^(1/3); % EqDiameter
                            SubGrain{jj}(k,7)=2*(3*DS.nVox(jj)*DS.VoxSize(1)*DS.VoxSize(2)*DS.VoxSize(3)/(4*pi))^(1/3); % EqDiameter of the parent grain
                        end
                    end
                end
            end
            jj
            %%%%% if binning is not needed
%             [x,y,z] = ind2sub(size(DS.GIDvol),find(DS.GIDvol == jj));
%             for k=1:length(x) % cell number,x,y,z,volume,EqDiameter,EqDiameter of the parent grain, unit based on mm
%                 BinnedPixel
%                 SubGrain{jj}(k,1)=k;
%                 SubGrain{jj}(k,2)=(x(k)-DS.Dimension(1)/2).*DS.VoxSize(1)+center_shift(1); % centroid coordinate x
%                 SubGrain{jj}(k,3)=(y(k)-DS.Dimension(2)/2).*DS.VoxSize(2)+center_shift(2); % centroid coordinate y
%                 SubGrain{jj}(k,4)=(z(k)-DS.Dimension(3)/2).*DS.VoxSize(3)+center_shift(3); % centroid coordinate z
%                 SubGrain{jj}(k,5)=DS.VoxSize(1)*DS.VoxSize(2)*DS.VoxSize(3); % volume
%                 SubGrain{jj}(k,6)=2*(3*DS.VoxSize(1)*DS.VoxSize(2)*DS.VoxSize(3)/(4*pi))^(1/3); % EqDiameter
%                 SubGrain{jj}(k,7)=2*(3*DS.nVox(jj)*DS.VoxSize(1)*DS.VoxSize(2)*DS.VoxSize(3)/(4*pi))^(1/3); % EqDiameter of the parent grain
%             end
        end
        
    else
        for jj=1:length(grain_measure(DS.SeedID))
            SubGrain{jj}(1)=1;
            SubGrain{jj}(2)=sampos(jj,1); % centroid coordinate x
            SubGrain{jj}(3)=sampos(jj,2); % centroid coordinate y
            SubGrain{jj}(4)=sampos(jj,3); % centroid coordinate z
            SubGrain{jj}(5)=grainvolume(jj)*1e-9; % volume
            SubGrain{jj}(6)=grainsize(jj)*1e-3; % EqDiameter
            SubGrain{jj}(7)=grainsize(jj)*1e-3; % EqDiameter of the parent grain
        end
    end
end

% Generate and calc SF or read list from file. By default, we generate hkl automatically
readhkl = 0;                 % generate reflections
if readhkl == 0
    sg = sglib(space_group_IT_number); % get sysconditions for specific element from the sglib.m
    sysconditions=sg.sysconditions;
end
structfact = 1; % do calculate structure factors

% make frame and store data output
makeframes = 1; % Make tif frames
direc = 'TFT';  % save frames in this directory
prefix = 'grains'; % prefix of frame names default is 'frame'
if ~exist('TFT', 'dir')
   mkdir('TFT'); % TFT folder is to store each output projection
end
if ~exist('DA', 'dir')
   mkdir('DA'); % DA folder is to store data record for each projection
end

% Add background counts and noise
addnoise = 0;            % Do not add Poissonian noise to the images
%addnoise = 1;           % Do
bgint = 10;              % initialize Background level in counts, will be recalculated later

% Define peak shape if diffraction images are to be formed 
% peakshape = 0;         % Represent peak as a spike AxA pixels, A depends on each polyhedron size 
peakshape = 1;       % Apply an anisotropic Gaussian point spread

save('input.mat');

