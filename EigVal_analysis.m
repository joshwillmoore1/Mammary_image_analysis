close all
clc
clear

load('PlantSeg_pipeline_outputs/R2_cleared/Cleared_R2_comp_membrane0015_PS_pipeOutput')
w1 = 0.01;
w2 = 1;

%%

AdjacencyMat = ConstructAdjacencyMatrix(CellData,w1,w2,1);

plotlinesX1 = [];
plotlinesX2 = [];
plotlinesY1 = [];
plotlinesY2 = [];
CentroidsX_lum = [];
CentroidsY_lum = [];
CentroidsX_basal = [];
CentroidsY_basal = [];



SizeOfCellData = size(CellData);
NumOfCells = SizeOfCellData(1);


for ii = 1:NumOfCells
    for jj = 1:NumOfCells
        
        
        if ismember(CellData(jj,1).LabelId,CellData(ii,1).AdjacentCells)
            
            plotlinesX1 = [plotlinesX1 , CellData(ii,1).Centroid(1)];
            plotlinesX2 = [plotlinesX2 , CellData(jj,1).Centroid(1)];
            plotlinesY1 = [plotlinesY1 , CellData(ii,1).Centroid(2)];
            plotlinesY2 = [plotlinesY2 , CellData(jj,1).Centroid(2)];
            
            
        end
        
    end
    
    if CellData(ii,1).cellType == 'Basal'
        CentroidsX_basal = [CentroidsX_basal;CellData(ii,1).Centroid(1)];
        CentroidsY_basal = [CentroidsY_basal;CellData(ii,1).Centroid(2)];
    elseif CellData(ii,1).cellType == 'Luminal'
        CentroidsX_lum = [CentroidsX_lum;CellData(ii,1).Centroid(1)];
        CentroidsY_lum = [CentroidsY_lum;CellData(ii,1).Centroid(2)];
    end
end





%% remove those cells with only homotypic connections
RemoveHomotypicCells = [];

for i = 1:NumOfCells
    
    
    AdjCells = cast(CellData(i,1).AdjacentCells,"double");
    
    ThisCellType = CellData(i,1).cellType;
    SameTypeCount = 0;
    
    for j = 1:length(AdjCells)
        for k = 1:NumOfCells
            
            if cast(CellData(k,1).LabelId,"double")  == AdjCells(j) && CellData(k,1).cellType == ThisCellType
                SameTypeCount = SameTypeCount + 1;
            end
            
        end
        
    end
    
    
    if cast(SameTypeCount,'double') == cast(length(AdjCells),'double')
        RemoveHomotypicCells = [RemoveHomotypicCells,i];
    end
    
end
AdjacencyMatHetero = AdjacencyMat;
AdjacencyMatHetero(RemoveHomotypicCells,:) = [];
AdjacencyMatHetero(:,RemoveHomotypicCells) = [];



%% check for laminar pattern eigvector  - do this in a separate script
%get row index of luminal cells
clc
close all

figure
Boundaries = imshow(CellBoundaries);
alpha(Boundaries,0.5)
hold on
plot([plotlinesX1;plotlinesX2 ],[plotlinesY1;plotlinesY2],'linewidth',1.5, 'color','k')
scatter(CentroidsX_basal,CentroidsY_basal,50,[0 0.4470 0.7410],'filled')
scatter(CentroidsX_lum,CentroidsY_lum,50,[0.8500 0.3250 0.0980],'filled')
axis equal
axis off
xlim([200,690])
title("Connectivity graph")




%%

TypeRowIndexSgn = ones(length(AdjacencyMat(:,1)),1);
CellNum = length(AdjacencyMat(:,1));

for i = 1:CellNum
    if CellData(i,1).cellType == "Basal"
        TypeRowIndexSgn(i) = -1;
    end
end

TransformRowMat = diag(TypeRowIndexSgn);
TransformRowMatHetero = TransformRowMat;
TransformRowMatHetero(RemoveHomotypicCells,:) =[]; 
TransformRowMatHetero(:,RemoveHomotypicCells) =[]; 

[EigVec, EigVals] = eig(AdjacencyMatHetero);

for i = 1:length(AdjacencyMatHetero(:,1))
   
    CheckSigns = TransformRowMatHetero*EigVec(:,i);
     
    PosEigsCount = length(AdjacencyMatHetero(:,1)) - sum(CheckSigns>0);
    NegEigsCount = length(AdjacencyMatHetero(:,1))- sum(CheckSigns<0);
    
    
    if sum(CheckSigns>0) == length(AdjacencyMatHetero(:,1))|| sum(CheckSigns<0) ==length(AdjacencyMatHetero(:,1))
        disp("A monotonic laminar pattern transformation exists")
        
    end
end


