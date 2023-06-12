function uniqueCellsToMeanCurvatures = ClosestCellCurvatureCount(cellData,shape_details,AdjMat,w1,w2)

X= shape_details.XY(:,1);
Y = shape_details.XY(:,2);
C = shape_details.curvature'*1;

 
    for j = 1:length(cellData(:,1))
        Centroids(j,:) = cellData(j,1).Centroid;
    end
    
    for i = 1:length(C)
        
        TempXY = [X(i),Y(i)];
        
        [minDist,minDistInd] = min(vecnorm(Centroids - TempXY,2,2));
   
        CellToCurv(i,:) = [minDistInd,C(i)];
        
        
    end
    
    UniqueClosestCells = unique(CellToCurv(:,1));
    
    uniqueCellsToMeanCurvatures = [];
    for i = 1:length(UniqueClosestCells)
        
        HomoCount = sum(AdjMat(UniqueClosestCells(i),:) == w1);
        HeteroCount = sum(AdjMat(UniqueClosestCells(i),:) == w2);
       
        CurvatureValmin = min( CellToCurv(CellToCurv(:,1) == UniqueClosestCells(i),2) );
        CurvatureValmax = max( CellToCurv(CellToCurv(:,1) == UniqueClosestCells(i),2) );
         if (abs(CurvatureValmin) > abs(CurvatureValmax))
            CurvatureVal  = CurvatureValmin;
         else
             CurvatureVal = CurvatureValmax;
         end
        
        %CurvatureVal = mean( CellToCurv(CellToCurv(:,1) == UniqueClosestCells(i),2) );
      
        
        uniqueCellsToMeanCurvatures(i,:) = [UniqueClosestCells(i),...
            CurvatureVal,...
            HomoCount,HeteroCount ];
        
    end
    
    
end