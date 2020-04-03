
function [Vsub SubGrain]=mesh3D(SingleGrain,MeshNr)
% % DipImage toolbox
% addpath('C:\Program Files\dipimage_2.9_win64\dip\common\dipimage');
% dip_initialise;
% dipsetpref('imagefilepath','C:\Program Files\dipimage_2.9_win64\dip\images');

% % iso2mesh toolbox
% addpath('C:\Users\fangh\Documents\MATLAB\iso2mesh');

% % This is for testing
% SingleGrain=Polyhedron(V.Set(1,1));
% MeshNr=8;

flag_figure=0; % 1 - plot; 0 - no plot

xi=linspace(min(SingleGrain.V(:,1)),max(SingleGrain.V(:,1)),MeshNr);
yi=linspace(min(SingleGrain.V(:,2)),max(SingleGrain.V(:,2)),MeshNr);
zi=linspace(min(SingleGrain.V(:,3)),max(SingleGrain.V(:,3)),MeshNr);
slack=(xi(2)-xi(1));
xi(1)=xi(1)+slack/2;
xi(end)=xi(end)-slack/2;
slack=(yi(2)-yi(1));
yi(1)=yi(1)+slack/2;
yi(end)=yi(end)-slack/2;
slack=(zi(2)-zi(1));
zi(1)=zi(1)+slack/2;
zi(end)=zi(end)-slack/2;

GridPointsIn=[];
GridPointsOut=[];
for i=1:length(xi)
    for j=1:length(yi)
        for k=1:length(zi)
%             xi(i)=xi(i)+(rand(1)*slack/4-slack/8);
%             yi(j)=yi(j)+(rand(1)*slack/4-slack/8);
%             zi(k)=zi(k)+(rand(1)*slack/4-slack/8);
            if SingleGrain.contains([xi(i) yi(j) zi(k)]')==1;
                GridPointsIn=[GridPointsIn;xi(i) yi(j) zi(k)];
            else
                GridPointsOut=[GridPointsOut;xi(i) yi(j) zi(k)];
            end
        end
    end
end
GridPointsIn=GridPointsIn(randperm(length(GridPointsIn)),:);

% GridPointsIn=SingleGrain.grid(MeshNr);
if ~isempty(GridPointsIn)
    [Vsub,Cellsub] = mpt_voronoi(GridPointsIn', 'bound', SingleGrain);
else
    Vsub=PolyUnion(SingleGrain);
end
for i=1:length(Vsub.Set)
    SubGrain(i,1)=i;
    SubGrain(i,2)=Vsub.Set(i,1).Data.voronoi.seed(1); % centroid coordinate x
    SubGrain(i,3)=Vsub.Set(i,1).Data.voronoi.seed(2); % centroid coordinate y
    SubGrain(i,4)=Vsub.Set(i,1).Data.voronoi.seed(3); % centroid coordinate z
    SubGrain(i,5)=Vsub.Set(i,1).volume; % volume
    if SubGrain(i,5)==Inf
        cc=Polyhedron(Vsub.Set(i,1).V);
        SubGrain(i,5)=cc.volume;
    end
    SubGrain(i,6)=2*(3*Vsub.Set(i,1).volume/(4*pi))^(1/3); % EqDiameter
    SubGrain(i,7)=2*(3*SingleGrain.volume/(4*pi))^(1/3); % EqDiameter of the parent grain
end
if flag_figure==1
    figure;
    subplot(2,2,1);
    hold all;
    SingleGrain.plot('FaceColor','cyan','alpha', 0.15,'edgecolor','w','linewidth',1); % Adjust the transparancy
    plot3(SingleGrain.V(:,1),SingleGrain.V(:,2),SingleGrain.V(:,3),'r.','MarkerSize',14);
    plot3(GridPointsIn(:,1),GridPointsIn(:,2),GridPointsIn(:,3),'b.','MarkerSize',11);
    % plot3(GridPointsOut(:,1),GridPointsOut(:,2),GridPointsOut(:,3),'b+','MarkerSize',6);
    box on;
    subplot(2,2,2);
    Vsub.plot('FaceColor','cyan','alpha',0.15,'edgecolor','k','linewidth',0.75);
    box on;
    subplot(2,2,3);
    hold all;
    for i=1:length(Vsub.Set)
%         Vsub.Set(i).plot('alpha', 0.1+i/100); % Adjust the transparancy
        Vsub.Set(i).plot('alpha', 0.1); % Adjust the transparancy
    end
    box on;
    subplot(2,2,4);
    histogram(SubGrain(:,6));
    xlabel('EqDiameter');
    ylabel('Number of cells');
    set(gca,'fontsize',14);
    set(gca,'linewidth',1.2);
    box on;
end







