% Refine the mesh size

clear all;

% % on cluster
% %%%% load dipimage toolbox
% %% After 2017, Linux has included dip_image 2.7 therefore it has installed
% %% the library data and added the path of dipimage.
% % dip_initialise;
% % call for dip_image 2.8 on Linux June 2019
% addpath('/opt/ud/dipimage-2.8/Linuxa64/lib/');
% addpath('/opt/ud/dipimage-2.8/common/dipimage');
% dip_initialise;
% %dipsetpref('imagefilepath','/opt/dipimage-2.3/images');
% addpath('/opt/ud/dipimage-2.8/common/dipimage/demos');
% 
% % start the MPT3 toolbox
% addpath('/home/hfang/Solute_drag_modelling/Linux_cnt_Fe_0_1C_0_5Mn_10Ks');
% tbxmanager restorepath;
% mpt_init;
% % cite using MPT3:
% % M. Herceg, M. Kvasnica, C.N. Jones, and M. Morari. Multi-Parametric Toolbox 3.0. In Proc. of the European Control Conference, pages 502?10, Zurich, Switzerland, July 17?9 2013.
% % load_mtex;
% load('/home/hfang/Linux_LabDCT/100um/Grain100um_400_400_600_input.mat'); % Linux path

% on local computer
% load('C:\Users\hfang\Documents\MATLAB\simul_LabDCT_2019\Input_Al_Dss-14_Dsd14_exp600s_4grains_April6_grainID_1\Input_6grains_MeshNr8.mat');
load('C:\Users\hfang\Documents\MATLAB\simul_LabDCT_2019_v2\Grain100um_400_400_600_input.mat');
% load('C:\Users\hfang\Documents\MATLAB\simul_LabDCT_2019_v2\Grain30um_100_100_150_input.mat');
% load('C:\Users\hfang\Documents\MATLAB\simul_LabDCT_2019_v2\Gradient_100_100_150_input.mat');
% load('C:\Users\hfang\Documents\MATLAB\simul_LabDCT_2019\Input_Al_Dss-14_Dsd14_exp600s_4grains_April6_grainID_1\Input_8grains_MeshNr_vary.mat');

ds0=2;
virtual_input=1;

if virtual_input==1
    %%%%%%%%%%%%%%%% remesh for virtual input
    MeshGrain=length(A_P(:,1)); % define how many grains you want to mesh for simulations
    for jj=1:MeshGrain
    % MeshGrain=53;
    % for jj=MeshGrain
       SingleGrain=Polyhedron(V_microstructure.Set(jj));
       MeshNr(jj)=round(SingleGrain.volume^(1/3)/mean([pixelysize pixelzsize].*1000));
       if MeshNr(jj)>10
           MeshNr(jj)=10;
       elseif MeshNr(jj)<2
           MeshNr(jj)=2;
       end
       [Vsub{jj} SubGrain{jj}]=mesh3D(SingleGrain,MeshNr(jj)); % cell number,x,y,z,volume,EqDiameter,EqDiameter of the parent grain
       SubGrain{jj}(:,2:4)=SubGrain{jj}(:,2:4)./1000; % [x,y,z] [mm]
       jj
    end
    save('ReMesh12_30um.mat');
end

if virtual_input==0
    %%%%%%%%%%%%%%%% remesh for DS file
    for jj=2:length(DS.SeedID)
        SingleGrain=Polyhedron(grain_vertices{jj});
%         P_grain(jj)=SingleGrain;
        MeshNr(jj)=round((SingleGrain.volume.*10^9)^(1/3)/mean([pixelysize pixelzsize].*1000));
        %             MeshNr(jj)=round(log(3*20/((SingleGrain.volume.*10^9)^(1/3)))/(-0.1));
        if MeshNr(jj)>15
            MeshNr(jj)=15;%9
        elseif MeshNr(jj)<6
            MeshNr(jj)=6;
        end
        [Vsub{jj} SubGrain{jj}]=mesh3D(SingleGrain,MeshNr(jj)); % cell number,x,y,z,volume,EqDiameter,EqDiameter of the parent grain, unit based on mm
        cell2Vol(jj)=sum(SubGrain{jj}(:,5))./(DS.nVox(jj)*DS.VoxSize(1)*DS.VoxSize(2)*DS.VoxSize(3)); % total volume
        jj
    end
    save('input_remesh.mat');
end


for k=1:length(SubGrain)
    meanD(k)=mean(SubGrain{k}(:,6));
end
mean(setdiff(meanD,Inf))



