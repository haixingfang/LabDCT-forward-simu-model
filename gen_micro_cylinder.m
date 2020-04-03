% This function is to generate 3D single phase microstructure for
% simulating 2D diffraction patterns
% Using voronoi cell to represent grains
% Haixing Fang, March 2019
% Cylindrical shape

function [A_P V]=gen_micro_cylinder(Lx,Lz,dgrain)
% % illuminated sample size, x-incoming beam; y-lateral; z-vertical
% Lx=100; % determined by sample thickness [um]
% Lz=150; % determined by beam size [um]
Ly=Lx; % determined by either beam size or sample width [um]
Rx=Lx/2;

% dgrain=30; % average austenite grain diameter [um]
dmin=0.5*dgrain;     % The minimum distance required between the centroids of grains
% Vbox=Lx*Ly*Lz;
Vbox=pi*Rx^2*Lz;
N_grain=fix(Vbox/(pi*(dgrain^3)/6)); % number of grains
rou_grain=N_grain/Vbox; % Number density of grain [um-3]

enlarge_flag=0; % 0: not to enlarge the box by d_grain/2; 1: enlarge the box by d_grain/2
if enlarge_flag==1
    enlarge_Lx=Lx+2*(0.5*1/rou_grain^(1/3));% enlarged by 2*average neibouring distance
    enlarge_Ly=Ly+2*(0.5*1/rou_grain^(1/3));% enlarged by 2*average neibouring distance
    enlarge_Lz=Lz+2*(0.5*1/rou_grain^(1/3));% enlarged by 2*average neibouring distance
    enlarge_N=fix(rou_grain*enlarge_Lx*enlarge_Ly*enlarge_Lz); % Number of grains in enlargement box (enlarged by 2*average neibouring distance) that have the same number density
    minus_edge=0.5*1/rou_grain^(1/3);
    enlarge_Rx=enlarge_Lx/2;
else
    enlarge_Lx=Lx;
    enlarge_Ly=Ly;
    enlarge_Lz=Lz;
    enlarge_N=N_grain;
    minus_edge=0;
    enlarge_Rx=enlarge_Lx/2;
end
eps=1e-6; % minimum error

gen_flag=1;
if gen_flag==1
    %Generate positions potential nucleation sites and plot the 3d voronoi diagram
    theta=rand(1)*360;
    A1=[enlarge_Rx*rand(1)*cosd(theta) enlarge_Rx*rand(1)*sind(theta) rand(1)*enlarge_Lz-enlarge_Lz/2]-minus_edge;
    i=1;
    A=A1;
    while i<=enlarge_N-1
       theta=rand(1)*360;
       A2(1)=enlarge_Rx.*rand(1)*cosd(theta)-minus_edge;
       A2(2)=enlarge_Rx.*rand(1)*sind(theta)-minus_edge;
       A2(3)=enlarge_Lz.*rand([1,1])-minus_edge-Lz/2;
       comb=[A;A2];
       distance=pdist(comb,'euclidean');
       if min(distance)>dmin
          A=[A;A2];
          i=i+1
       end
    end
    % export the predefined grain centroid coordinates
    dlmwrite('SeedingPoints_cylinder.txt',[A(:,1) A(:,2) A(:,3)],'delimiter',' ');
else
    % load the predefined grain centroid coordinates
%    DataFolder='C:\Users\hfang\Documents\MATLAB\simul_LabDCT_2019_v2';
	DataFolder=pwd;
    % FileName_prefix='SeedingPoints_cylinder'; % randomly generated
%     FileName_prefix='Grain80um_400_400_600';
    FileName_prefix='Grain30um_100_100_150';
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

dis=pdist(A,'euclidean'); % Pairwise distance of austenite grain center
[X,Y,Z] = cylinder(1,50);
X=X.*(enlarge_Rx);
Y=Y.*(enlarge_Rx);
Z=Z.*(enlarge_Lz)-enlarge_Lz/2;
Zstep=50;
for i=1:fix((enlarge_Lz)/Zstep)-1
    X=[X;X(1,:)];
    Y=[Y;Y(1,:)];
    Z=[Z;ones(1,length(Z(1,:))).*Zstep*i+Z(1,1)];
    if i==fix((enlarge_Lz-minus_edge)/Zstep)-1
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
BoundaryVertices=[BoundaryVertices;0 0 enlarge_Lz-minus_edge-Lz/2;0 0 -minus_edge-Lz/2];
B=Polyhedron(BoundaryVertices); % Boundary vertices
[V,cells] = mpt_voronoi(A', 'bound', B); % Call for the mpt_voronoi
vertices = V.Set.forEach(@(e) e.V, 'UniformOutput', false);% find the vertices for each cell
new_vertices = cat(1, vertices{:}); % Combine the vertices into one matrix
potential0=unique(round(new_vertices*1e6), 'rows')/1e6; % Get rid of the numerical noise 
to_remove = [-Lx/2 -Ly/2 -Lz/2;Lx/2 -Ly/2 -Lz/2;Lx/2 Ly/2 -Lz/2;-Lx/2 Ly/2 -Lz/2; ...
    -Lx/2 -Ly/2 Lz/2;Lx/2 -Ly/2 Lz/2;Lx/2 Ly/2 Lz/2;-Lx/2 Ly/2 Lz/2];% Remove the corners of the Lx*Ly*Lz cubic
potential=setdiff(potential0, to_remove, 'rows');
%potential=potential0;
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

%{
figure('Name','Grain geometry');
subplot(1,2,1);
V.plot('alpha', 0.13); % Adjust the transparancy
hold all;
plot3(A(:,1),A(:,2),A(:,3),'b+','LineWidth',2);% Plot the centroids
plot3(potential0(:,1),potential0(:,2),potential0(:,3),'r.','MarkerSize',20); % Plot all the corners
axis([-enlarge_Lx/2 enlarge_Lx/2 -enlarge_Ly/2 enlarge_Ly/2 -enlarge_Lz/2 enlarge_Lz/2]);
xlabel('x (\mum)','FontSize',16);
ylabel('y (\mum)','FontSize',16);
zlabel('z (\mum)','FontSize',16);
set(get(gca,'xlabel'),'rotation',18);
set(get(gca,'ylabel'),'rotation',-25);
set(get(gca,'zlabel'),'rotation',90);
set(gca,'fontsize',14);
set(gca,'linewidth',2);
% legend('+ centroids of austenite grains','. corners of austenite grains');
box on;
grid off;
hold off;

subplot(1,2,2);
V.plot('alpha', 0.13); % Adjust the transparancy
hold all;
plot3(A(:,1),A(:,2),A(:,3),'b+','LineWidth',2);% Plot the centroids
plot3(site(:,1),site(:,2),site(:,3),'r.','MarkerSize',20);% Plot the corners
% plot3(position(:,1),position(:,2),position(:,3),'b*','LineWidth',2);
grid off;
xlabel('x (\mum)','FontSize',16);
ylabel('y (\mum)','FontSize',16);
zlabel('z (\mum)','FontSize',16);
axis([-Lx/2 Lx/2 -Ly/2 Ly/2 -Lz/2 Lz/2]);
% axesLabelsAlign3D;
set(get(gca,'xlabel'),'rotation',18);
set(get(gca,'ylabel'),'rotation',-25);
set(get(gca,'zlabel'),'rotation',90);
set(gca,'fontsize',14);
set(gca,'linewidth',2);
box on;
legend('+ centroids of austenite grains','. potential ferrite nucleation sites','* randomly selected nucleation sites');
hold off;
%}
figure('Name','Grain geometry');
pos = [0.1 0.1 0.85 0.85]; % [left bottom width height]
subplot('Position',pos);
V.plot('alpha', 0.5,'edgecolor','k','linestyle','none','linewidth',0.75,'colororder','fixed'); % Adjust the transparancy
hold all;
theta=0:0.01:2*pi;
% plot3((enlarge_Rx-minus_edge).*cos(theta),(enlarge_Rx-minus_edge).*sin(theta), ...
%     (enlarge_Lz-minus_edge-Lz/2).*ones(size((enlarge_Rx-minus_edge).*cos(theta))),'k-');
% plot3((enlarge_Rx-minus_edge).*cos(theta),(enlarge_Rx-minus_edge).*sin(theta), ...
%     (-minus_edge-Lz/2).*ones(size((enlarge_Rx-minus_edge).*cos(theta))),'k-');
plot3((enlarge_Rx).*cos(theta),(enlarge_Rx).*sin(theta), ...
    (enlarge_Lz-minus_edge-Lz/2).*ones(size((enlarge_Rx).*cos(theta))),'k-');
plot3((enlarge_Rx).*cos(theta),(enlarge_Rx).*sin(theta), ...
    (-minus_edge-Lz/2).*ones(size((enlarge_Rx).*cos(theta))),'k-');

% plot3([-(enlarge_Rx-minus_edge) -(enlarge_Rx-minus_edge)], [0 0], ...
%     [(-minus_edge-Lz/2) (enlarge_Lz-minus_edge-Lz/2)],'k-');
% plot3([(enlarge_Rx-minus_edge) (enlarge_Rx-minus_edge)], [0 0], ...
%     [(-minus_edge-Lz/2) (enlarge_Lz-minus_edge-Lz/2)],'k-');
% plot3([0 0],[-(enlarge_Rx-minus_edge) -(enlarge_Rx-minus_edge)], ...
%     [(-minus_edge-Lz/2) (enlarge_Lz-minus_edge-Lz/2)],'k-');
% plot3([0 0],[(enlarge_Rx-minus_edge) (enlarge_Rx-minus_edge)], ...
%     [(-minus_edge-Lz/2) (enlarge_Lz-minus_edge-Lz/2)],'k-');

% hold all;
% for kk=1:length(V.Set)
% %    V.Set(kk).plot('FaceColor','cyan','alpha',0.15,'edgecolor','k','linewidth',0.75); % Adjust the transparancy
% %    V.Set(kk).plot('alpha',0.15,'edgecolor','k','linewidth',0.75); % Adjust the transparancy
%    V.Set(kk).plot('alpha',0.15); % Adjust the transparancy
% end
grid off;
xlabel('x (\mum)','FontSize',16);
ylabel('y (\mum)','FontSize',16);
zlabel('z (\mum)','FontSize',16);
axis equal;
axis([-enlarge_Lx/2 enlarge_Lx/2 -enlarge_Ly/2 enlarge_Ly/2 -enlarge_Lz/2 enlarge_Lz/2]);
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

% % 2D slice
% figure;
% S = V.slice(3,200);
% S.plot('color','blue','alpha',0.2,'linestyle','-','linewidth',1.5);

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

% figure;
% aa = V.convexHull;
% aa.plot('wire',true,'linewidth',3,'linestyle','--')

% % N_p contains properties of the vertex
% for i=1:length(site)
%     N_p(i,1)=i;%site number
%     N_p(i,2)=site(i,1);%x coordinate
%     N_p(i,3)=site(i,2);%y coordinate
%     N_p(i,4)=site(i,3);%z coordinate
%     N_PR{i}=[];
%     for j=1:length(V.Set)
%       for k=1:length(V.Set(j).V)
%        if abs(site(i,1)-V.Set(j).V(k,1))<=eps && abs(site(i,2)-V.Set(j).V(k,2))<=eps && abs(site(i,3)-V.Set(j).V(k,3))<=eps
%           N_PR{i}=[N_PR{i}; j k]; % derive the corresponding relationship between vertex and austenite grains
%        end
%       end
%     end
%     if ~isempty(N_PR{i})
%         for m=1:length(N_PR{i}(:,1))
%             N_PR{i}(m,3)=sqrt((site(i,1)-A_P(N_PR{i}(m,1),2))^2+(site(i,2)-A_P(N_PR{i}(m,1),3))^2+(site(i,3)-A_P(N_PR{i}(m,1),4))^2);% identical distance to neighbouring centers
%         end
%         N_p(i,5)=mean(N_PR{i}(:,3)); % average distance between vertex and austenite centers [um]
%     else
%         N_p(i,5)=0; % average distance between vertex and austenite centers [um]
%     end
% end
% for j=1:length(V.Set)
%     A_PR{j}=[];
%     for k=1:length(V.Set(j).V)
%         for i=1:length(N_p(:,1))
%            if abs(N_p(i,2)-V.Set(j).V(k,1))<=eps && abs(N_p(i,3)-V.Set(j).V(k,2))<=eps && abs(N_p(i,4)-V.Set(j).V(k,3))<=eps
% %               A_PR{j}=[A_PR{j};N_p(i,1)]; % restore the site order number
%               A_PR{j}(k)=N_p(i,1); % restore the site order number
%            end
%         end
%     end
% end


% model = createpde(1);
% geometryFromEdges(model,@lshapeg);

