close all
clc
clear

close all
clc
clear
NumberOfPoints = 10000;
t = linspace(0,2*pi,NumberOfPoints);


numberLum = 60;
numberBas = 51;
t_l = linspace(0,2*pi,numberLum);
t_b = linspace(0,2*pi,numberBas);
tspaced = [t,t,t];
L = 1;
A = 1.8; %0 0.25 1.2 1.7
r = L*5;
n = 4;
Rl = 1.2*L;
Centre = [0,0];
RotConst = -0.1;
LumScaleR = 0.82;



NumberOfGhosts = 1000;
GhostDistLum = 0.5;
GhostDistBas = 0.2;

tg = linspace(0,2*pi,NumberOfGhosts);


for i = 1:length(t_l)
    
    xRl(i) = LumScaleR.*R(A,r,n,t_l(i)+RotConst).*cos(t_l(i) + RotConst)+ Centre(1);
    yRl(i) = LumScaleR.*R(A,r,n,t_l(i)+RotConst).*sin(t_l(i) + RotConst)+ Centre(2); 
end

for i = 1:length(t_b)
    
    xRb(i) = 1*R(A,r,n,t_b(i)).*cos(t_b(i))+ Centre(1);
    yRb(i) = 1*R(A,r,n,t_b(i)).*sin(t_b(i))+ Centre(2);
end

InnerGhostNodes = [];
OuterGhostNodes = [];

for i = 1:NumberOfGhosts
    InnerGhostNodes = [InnerGhostNodes; GhostDistLum.*R(1.5*A,r,n,tg(i)).*[ cos(tg(i))+ Centre(1),sin(tg(i))+ Centre(2) ] ];
    OuterGhostNodes = [OuterGhostNodes; (1 + GhostDistBas).*R(0.5*A,r,n,tg(i)).*[ cos(tg(i))+ Centre(1),sin(tg(i))+ Centre(2) ] ];
end

FullCells = [[xRb',yRb'];[xRl',yRl']];
FullCellsWithGhost = [FullCells;InnerGhostNodes;OuterGhostNodes];


DT_orginal = delaunayTriangulation(FullCells);
[VOrg,rOrg] = voronoiDiagram(DT_orginal);
numOfBoundedPolygons = length(rOrg);


DT = delaunayTriangulation(FullCellsWithGhost);
[V,r] = voronoiDiagram(DT);
Colors = lines(2);
figure;
hold on
for i = 1:numOfBoundedPolygons
    if i <= numberLum
        poly = V(r{i},:);
        patch(poly(:,1),poly(:,2),[1,1,1],'linewidth',2)
    else
        
        poly = V(r{i},:);
        patch(poly(:,1),poly(:,2),[1,1,1],'linewidth',2)
    end
end

%scatter(xRb,yRb,'filled')
hold on
%scatter(xRl,yRl,'filled')
%scatter(InnerGhostNodes(:,1),InnerGhostNodes(:,2),'filled')
%scatter(OuterGhostNodes(:,1),OuterGhostNodes(:,2),'filled')

axis equal
axis off

%%
function out = R(A,r,n,t)

out  =  r + A.*sin(n.*t);

end

function out = Rp(A,r,n,t)

out  =  A.*n.*cos(n.*t);

end

function out = Rpp(A,r,n,t)

out  =  -A.*(n.^2).*sin(n.*t);

end