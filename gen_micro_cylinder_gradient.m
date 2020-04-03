% This function is to generate 3D single phase microstructure for
% simulating 2D diffraction patterns
% Using voronoi cell to represent grains
% Haixing Fang, March 2019
% Cylindrical shape

function [A_P V]=gen_micro_cylinder_gradient(Lx,Lz,gr_diameter)

% % illuminated sample size, x-incoming beam; y-lateral; z-vertical
% Lx=100; % determined by sample thickness [um]
% Lz=150; % determined by beam size [um]
Ly=Lx; % determined by either beam size or sample width [um]
Rx=Lx/2;

% average grain diameter as a function of depth
Ldepth=[0:7.5:Lz/3 Lz/3+5:10:Lz*2/3 Lz*2/3:15:Lz];
Ddepth=-5.82373e-03.*Ldepth.^2+1.68125.*Ldepth-0.210396;
Ddepth(1)=gr_diameter; % initialize the grain diamter for the top lay [um]
A=[];
enlarge_flag=0;
gen_flag=1;
for i=1:length(Ldepth)-1
    dgrain(i)=(Ddepth(i)+Ddepth(i+1))/2; 
    dmin(i)=0.25*dgrain(i);
    Vbox(i)=pi*Rx^2*(Ldepth(i+1)-Ldepth(i));
    N_grain(i)=fix(Vbox(i)/(pi*(dgrain(i)^3)/6)); % number of grains
    if N_grain(i)<2
        N_grain(i)=2;
    end
    rou_grain(i)=N_grain(i)/Vbox(i); % Number density of grain [um-3]
    eps=1e-6; % minimum error
    if gen_flag==1;
        %Generate positions potential nucleation sites and plot the 3d voronoi diagram
        theta=rand(1)*360;
        A1=[Rx*rand(1)*cosd(theta) Rx*rand(1)*sind(theta) rand(1)*(Ldepth(i+1)-Ldepth(i))+Ldepth(i)];
        k=1;
        A_depth{i}=A1;
        while k<=N_grain(i)-1
           theta=rand(1)*360;
           A2(1)=Rx.*rand(1)*cosd(theta);
           A2(2)=Rx.*rand(1)*sind(theta);
           A2(3)=(Ldepth(i+1)-Ldepth(i)).*rand([1,1])+Ldepth(i);
           comb=[A_depth{i};A2];
           distance=pdist(comb,'euclidean');
           if min(distance)>dmin(i)
              A_depth{i}=[A_depth{i};A2];
              k=k+1;
           end
        end
        A=[A;A_depth{i}];
        if i==length(Ldepth)-1
            A(:,3)=A(:,3)-Lz/2;
            A(:,3)=-A(:,3);
            % export the predefined grain centroid coordinates
            dlmwrite('SeedingPoints_cylinder_gradient.txt',[A(:,1) A(:,2) A(:,3)],'delimiter',' ');
        end
    else
        % load the predefined grain centroid coordinates
        % DataFolder='C:\Users\hfang\Documents\MATLAB\simul_LabDCT_2019_v2';
		DataFolder=pwd;
        % FileName_prefix='SeedingPoints_cylinder'; % randomly generated
        FileName_prefix='Gradient_100_100_150';
        FileName_pattern='.txt';
        baseFileName=[FileName_prefix FileName_pattern];
        fullFileName=fullfile(DataFolder, baseFileName);
        fileID=fopen(fullFileName,'r');
        A=[];
        while(~feof(fileID))
            textdata=str2num(fgetl(fileID));
            if isempty(textdata)
                continue;
            else
                A=[A; textdata];
            end
        end
        fclose(fileID);
    end   
   i
end

dis=pdist(A,'euclidean'); % Pairwise distance of austenite grain center
[X,Y,Z] = cylinder(1,50);
X=X.*(Rx);
Y=Y.*(Rx);
Z=Z.*(Lz)-Lz/2;
Zstep=50;
for i=1:fix((Lz)/Zstep)-1
    X=[X;X(1,:)];
    Y=[Y;Y(1,:)];
    Z=[Z;ones(1,length(Z(1,:))).*Zstep*i+Z(1,1)];
    if i==fix(Lz/Zstep)-1
        Zz=Z(2,:);
        Z(2:end-1,:)=Z(3:end,:);
        Z(end,:)=Zz;
    end
end
BoundaryVertices=[];
for i=1:length(X(:,1))
    for j=1:length(X(1,:))
        BoundaryVertices=[BoundaryVertices;X(i,j) Y(i,j) Z(i,j)];
    end
end
BoundaryVertices=[BoundaryVertices;0 0 Lz-Lz/2;0 0 -Lz/2];
B=Polyhedron(BoundaryVertices); % Boundary vertices
[V,cells] = mpt_voronoi(A', 'bound', B); % Call for the mpt_voronoi
vertices = V.Set.forEach(@(e) e.V, 'UniformOutput', false);% find the vertices for each cell
new_vertices = cat(1, vertices{:}); % Combine the vertices into one matrix
potential0=unique(round(new_vertices*1e6), 'rows')/1e6; % Get rid of the numerical noise 
to_remove = [-Lx/2 -Ly/2 -Lz/2;Lx/2 -Ly/2 -Lz/2;Lx/2 Ly/2 -Lz/2;-Lx/2 Ly/2 -Lz/2; ...
    -Lx/2 -Ly/2 Lz/2;Lx/2 -Ly/2 Lz/2;Lx/2 Ly/2 Lz/2;-Lx/2 Ly/2 Lz/2];% Remove the corners of the Lx*Ly*Lz cubic
potential=setdiff(potential0, to_remove, 'rows');
if enlarge_flag==1
    for i=1:length(potential)
        for j=1:3
            if j==1
                Lb=Lx;
            else if j==2
                    Lb=Ly;
                else if j==3
                        Lb=Lz;
                    end
                end
            end
        if potential(i,j)>Lb/2
           potential(i,j)=Lb/2+1;
        else if potential(i,j)<-Lb/2
            potential(i,j)=-Lb/2-1;
            end
        end
        end
    end
end
site=[];
% Remove points beyond the Lx*Ly*Lz cubic box
for i=1:length(potential)
  %if all(potential(i,:))&&all(potential(i,:)-Lb-1)
  if enlarge_flag==1
      if (all(potential(i,1)-Lx/2-1)&&all(potential(i,1)+Lx/2+1)) && (all(potential(i,2)-Ly/2-1) ...
              &&all(potential(i,2)+Ly/2+1)) (all(potential(i,3)-Lz/2-1)&&all(potential(i,3)+Lz/2+1))
         eval(['site',num2str(i),'=','potential(i,:)']); % Convert the number to string
         eval(['site=[site;site',num2str(i),'];']);  % Combine site1, site 2,...
      end
  else
        if (all(potential(i,1)-Lx/2)&&all(potential(i,1)+Lx/2)) && (all(potential(i,2)-Ly/2) ...
          &&all(potential(i,2)+Ly/2)) (all(potential(i,3)-Lz/2)&&all(potential(i,3)+Lz/2))
     eval(['site',num2str(i),'=','potential(i,:)']); % Convert the number to string
     eval(['site=[site;site',num2str(i),'];']);  % Combine site1, site 2,...
        end
  end
end
p1=randperm(length(site));
site=site(randperm(length(site)),:); % randomly sort the possible nucleation sites

figure('Name','Grain geometry');
pos = [0.1 0.1 0.85 0.85]; % [left bottom width height]
subplot('Position',pos);
V.plot('alpha', 0.5,'edgecolor','k','linestyle','none','linewidth',0.75,'colororder','fixed'); % Adjust the transparancy
hold all;
theta=0:0.01:2*pi;
plot3((Rx).*cos(theta),(Rx).*sin(theta), ...
    (Lz-Lz/2).*ones(size((Rx).*cos(theta))),'k-');
plot3((Rx).*cos(theta),(Rx).*sin(theta), ...
    (-Lz/2).*ones(size((Rx).*cos(theta))),'k-');
grid off;
xlabel('x (\mum)','FontSize',16);
ylabel('y (\mum)','FontSize',16);
zlabel('z (\mum)','FontSize',16);
axis equal;
axis([-Lx/2 Lx/2 -Ly/2 Ly/2 -Lz/2 Lz/2]);
set(get(gca,'xlabel'),'rotation',30);
set(get(gca,'ylabel'),'rotation',-30);
set(get(gca,'zlabel'),'rotation',90);
view(3);
set(gca,'fontsize',14);
set(gca,'linewidth',1.5);
box on;
% legend('+ centroids of austenite grains','. potential ferrite nucleation sites');
% % legend('+ centroids of austenite grains','. potential ferrite nucleation sites','* randomly selected nucleation sites');
hold off;
% print('GrainGeometry','-dpng');

% 2D slice
figure;
Sxoy = V.slice(3,0); % XOY plane
Sxoz = V.slice(2,0); % XOZ plane
Syoz = V.slice(1,0); % YOZ plane
subplot(2,2,1);
V.plot('alpha', 0.5,'edgecolor','k','linestyle','none','linewidth',0.75,'colororder','fixed'); % Adjust the transparancy
hold all;
plot3((Rx).*cos(theta),(Rx).*sin(theta), ...
    (Lz-Lz/2).*ones(size((Rx).*cos(theta))),'k-');
plot3((Rx).*cos(theta),(Rx).*sin(theta), ...
    (-Lz/2).*ones(size((Rx).*cos(theta))),'k-');
% Sxoz.plot('color','red','alpha',0.2,'linestyle','-','linewidth',1.5);
% Syoz.plot('color','magenta','alpha',0.2,'linestyle','-','linewidth',1.5);
grid off;
xlabel('x (\mum)','FontSize',16);
ylabel('y (\mum)','FontSize',16);
zlabel('z (\mum)','FontSize',16);
axis equal;
axis([-Lx/2 Lx/2 -Ly/2 Ly/2 -Lz/2 Lz/2]);
set(get(gca,'xlabel'),'rotation',30);
set(get(gca,'ylabel'),'rotation',-30);
set(get(gca,'zlabel'),'rotation',90);
view(3);
set(gca,'fontsize',14);
set(gca,'linewidth',1.5);
box on;
subplot(2,2,2);
Sxoy.plot('color','blue','alpha',0.2,'linestyle','-','linewidth',1.5);
xlabel('x (\mum)','FontSize',16);
ylabel('y (\mum)','FontSize',16);
axis([-Lx/2 Lx/2 -Ly/2 Ly/2]);
axis equal;
set(gca,'fontsize',14);
box on;
title('XOY plane');
subplot(2,2,3);
Sxoz.plot('color','red','alpha',0.2,'linestyle','-','linewidth',1.5);
xlabel('x (\mum)','FontSize',16);
ylabel('z (\mum)','FontSize',16);
axis([-Lx/2 Lx/2 -Lz/2 Lz/2]);
axis equal;
set(gca,'fontsize',14);
box on;
title('XOZ plane');
subplot(2,2,4);
Syoz.plot('color','magenta','alpha',0.2,'linestyle','-','linewidth',1.5);
xlabel('y (\mum)','FontSize',16);
ylabel('z (\mum)','FontSize',16);
axis([-Ly/2 Ly/2 -Lz/2 Lz/2]);
axis equal;
set(gca,'fontsize',14);
box on;
title('YOZ plane');

% calculate the average diameter as a function of depth
sliceZ=Lz/2;
slice_step=-2;
i=1;
while sliceZ>=-Lz/2
    Slice=V.slice(3,sliceZ);
    slice_depth(i)=-sliceZ+Lz/2;
    D_grain(i)=mean(2*(Slice.Set.volume/pi).^(1/2));
    sliceZ=sliceZ+slice_step;
    i=i+1;
end

% 2D slice
figure;
Sxoy = V.slice(3,0); % XOY plane
Sxoz = V.slice(2,0); % XOZ plane
Syoz = V.slice(1,0); % YOZ plane
subplot(1,3,1);
V.plot('alpha', 0.5,'edgecolor','k','linestyle','none','linewidth',0.75,'colororder','fixed'); % Adjust the transparancy
hold all;
plot3((Rx).*cos(theta),(Rx).*sin(theta), ...
    (Lz-Lz/2).*ones(size((Rx).*cos(theta))),'k-');
plot3((Rx).*cos(theta),(Rx).*sin(theta), ...
    (-Lz/2).*ones(size((Rx).*cos(theta))),'k-');
% Sxoz.plot('color','red','alpha',0.2,'linestyle','-','linewidth',1.5);
% Syoz.plot('color','magenta','alpha',0.2,'linestyle','-','linewidth',1.5);
grid off;
xlabel('x (\mum)','FontSize',16);
ylabel('y (\mum)','FontSize',16);
zlabel('z (\mum)','FontSize',16);
axis equal;
axis([-Lx/2 Lx/2 -Ly/2 Ly/2 -Lz/2 Lz/2]);
set(get(gca,'xlabel'),'rotation',30);
set(get(gca,'ylabel'),'rotation',-30);
set(get(gca,'zlabel'),'rotation',90);
view(3);
set(gca,'fontsize',14);
set(gca,'linewidth',1.5);
box on;
subplot(1,3,2);
Sxoz.plot('color','red','alpha',0.2,'linestyle','-','linewidth',1.5);
xlabel('x (\mum)','FontSize',16);
ylabel('z (\mum)','FontSize',16);
axis([-Lx/2 Lx/2 -Lz/2 Lz/2]);
axis equal;
set(gca,'fontsize',14);
box on;
% title('XOZ plane');
subplot(1,3,3);
plot(slice_depth,D_grain,'ro-');
rot=-90;
camroll(rot);
xlabel('Distance to top surface (\mum)','FontSize',16);
ylabel('Average grain diameter (\mum)','FontSize',16);
set(get(gca,'xlabel'),'rotation',rot);
set(gca,'XTickLabelRotation',rot);
set(gca,'YTickLabelRotation',rot);
xlim([0 Lz]);
ylim([0 5*(floor(max(D_grain)/5)+1)]);
set(gca,'fontsize',14);
box on;

% Define initial features for each austenite
corr=0;
for i=1:length(V.Set)
   A_P(i,1)=i; % order number
   A_P(i,2)=V.Set(i).Data.voronoi.seed(1,1); % x coordinate
   A_P(i,3)=V.Set(i).Data.voronoi.seed(2,1); % y coordinate
   A_P(i,4)=V.Set(i).Data.voronoi.seed(3,1); % z coordinate
   A_P(i,5)=V.Set(i).volume; % volume [um3]
   if A_P(i,5)==Inf
       cc=Polyhedron(V.Set(i).V);
       A_P(i,5)=cc.volume;
       corr=corr+1
   end
   A_P(i,6)=0; % number of nucleated ferrite in its vertex
   A_P(i,7)=length(V.Set(i).V); % total number of vertex
   A_P(i,8)=0;
   for j=1:length(V.Set(i).V)
       if abs(V.Set(i).V(j,1))>Lx/2 || abs(V.Set(i).V(j,2))>Ly/2 || abs(V.Set(i).V(j,3))>Lz/2
          A_P(i,8)=A_P(i,8)+1; % total number of vertex outside [Lx Ly Lz]
       end
   end
   A_P(i,9)=A_P(i,7)-A_P(i,8); % effective number of vertex(<[Lx Ly Lz])
   
   % randomly generate Euler angles [deg]: ([0 2*pi], [0 pi], [0 2*pi])
   A_P(i,10) = 360*rand(1,1);
   A_P(i,11) = 180*rand(1,1);
   A_P(i,12) = 360*rand(1,1);
   A_P(i,13)=2*(3*A_P(i,5)/(4*pi))^(1/3); % equivalent diameter [um]
end

