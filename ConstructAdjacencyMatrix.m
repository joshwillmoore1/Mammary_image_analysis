function A = ConstructAdjacencyMatrix(CellDataInput,w1,w2,RowStochBool)

%this function assumes the cell data structure

%cellData->(LabelID,Area,PixelLocations,Perimeter,Centroid,cellType,AdjacenctCells)
SizeOfCellData = size(CellDataInput);
NumOfCells = SizeOfCellData(1);


A = zeros(NumOfCells,NumOfCells);

for ii = 1:NumOfCells
    for jj = 1:NumOfCells
        
        if ismember(CellDataInput(jj,1).LabelId,CellDataInput(ii,1).AdjacentCells)
            
            if CellDataInput(ii,1).cellType == CellDataInput(jj,1).cellType
                A(ii,jj) =  w1;
                
            else
                A(ii,jj) =  w2;
            end
            
        end
        
        
    end
    
    if RowStochBool == 1
    %normalise for row-stochastic
    A(ii,:) = (1./(sum(A(ii,:)))).*(A(ii,:));
    
    end
    
end

end