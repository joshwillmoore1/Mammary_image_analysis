function [ BasalCellHomoCount,BasalCellHeteroCount,...
    LuminalCellHomoCount,LuminalCellHeteroCount,RemoveCells] = TypeDependentConnectivtityCounter(cellData,AdjMat,w1,w2)
    

    BasalCellHomoCount = [];
    BasalCellHeteroCount = [];
    
    LuminalCellHomoCount = [];
    LuminalCellHeteroCount = [];
    RemoveCells = [];
    
    
    
        
    for i = 1:length(cellData(:,1))
        
        if cellData(i).Area ~= 0
            
            if cellData(i).cellType == "Basal"
                BasalCellHomoCount = [BasalCellHomoCount,sum(AdjMat(i,:) == w1)];
                BasalCellHeteroCount = [BasalCellHeteroCount,sum(AdjMat(i,:) == w2)];
                
            elseif cellData(i).cellType == "Luminal"
                
                LuminalCellHomoCount  = [LuminalCellHomoCount ,sum(AdjMat(i,:) == w1)];
                LuminalCellHeteroCount = [LuminalCellHeteroCount,sum(AdjMat(i,:) == w2)];
                
            end
            
        else
            RemoveCells = [RemoveCells,i];
            
        end
        
    end
    
    
    
    
end