close all
clc
clear


PlotHistograms = 0;

clims = [-5e-2, 5e-2];
w1 = 0.5;
w2 = 1;
%BoundaryConst = 143.102;
BoundaryConst = 5;
XSCALE_dapi2018_scaled = 0.225;
XSCALE_dapi_GFP2_scaled = 0.225;
XSCALE_LargeN1_scaled = 0.225;
XSCALE_R2_cleared = 0.2203746;
XSCALE_R3_1_cleared = 0.2328895;

plotR2cleared = 1;
plotR3_1cleared = 1;
plot2018Dapicleared = 1;
plotGFP2cleared = 1;
plotLargeN1cleared = 1;

%ONLY IF ALL ARE ON ^^^
GroupedCurvatureAnalysis = 1;


oragnoidArea = [];
LumCellArea = [];
BasCellArea = [];


R2_basal_count = 0;
R2_luminal_count = 0;

R3_1_basal_count = 0;
R3_1_luminal_count = 0;

r2018_basal_count = 0;
r2018_luminal_count = 0;

GFP2_basal_count = 0;
GFP2_luminal_count = 0;

LN1_basal_count = 0;
LN1_luminal_count = 0;


%% R2 Cleared


if plotR2cleared  == 1
    SliceOfInterest = 16;
    
    
    PlantSegData_R2 = load('PlantSeg_pipeline_outputs/R2_cleared/Cleared_R2_comp_membrane0015_PS_pipeOutput.mat');
    
    
    
    
    
    twoDSliceR2 = imread("../PlantSeg/Cleared_R2/Cleared_R2_comp_membrane0015.tif");
    
    BinImageR2 = imbinarize(twoDSliceR2,0.1);
    FilterImageR2 = medfilt2(BinImageR2);
    FillImageR2 = imfill(FilterImageR2,'holes');
    SmoothImR2 = bwareaopen(FillImageR2,50);
    N = 9;
    kernel = ones(N, N) / N^3;
    blurryImage = conv2(double(SmoothImR2), kernel, 'same');
    newBinaryImageR2 = blurryImage>0.1 ;
    newBinaryImageR2 = imdilate(newBinaryImageR2, ones(3));
    
    oragnoidArea = [oragnoidArea,bwarea(newBinaryImageR2).*XSCALE_R2_cleared.^2];
    %check the input image
    % figure;
    % subplot(1,2,1)
    % imshow(twoDSliceR2)
    % subplot(1,2,2)
    % imshow(newBinaryImageR2)
    % hold off
    
    
    boundaryPoint = round(BoundaryConst/XSCALE_R2_cleared)    ;   %12, 10 number of boundary points curvature is found over
    curvatureThresh = 0.06;     %0.06 the maximum allowed value of the curvature measure
    bp_tangent = 10;            % number of boundary points the tangent angle is found over
    interpdmin = 0.3;           % the minimum number of pixels seperating boundary points after interpolation
    loopclose = 1;   % 0 - if open boundaries | 1 - if closed boundaries
    
    % Find the curvature
    [shape_details_R2, Icurv_R2] = curvature2D(newBinaryImageR2',...
        boundaryPoint, curvatureThresh, bp_tangent, interpdmin, loopclose);
    shape_details_R2.curvature = shape_details_R2.curvature.*(1/XSCALE_R2_cleared);
    %Plot curvature
    X_R2 = shape_details_R2.XY(:,1);
    Y_R2 = shape_details_R2.XY(:,2);
    Z_R2 = zeros(size(X_R2));
    C_R2 = shape_details_R2.curvature'*1;
    % Plot
    
    figure;
    imshow(1-PlantSegData_R2.Bound_image);
    hold on
    p3 = plot([PlantSegData_R2.plotlinesX1;PlantSegData_R2.plotlinesX2 ],...
        [PlantSegData_R2.plotlinesY1;PlantSegData_R2.plotlinesY2],'linewidth',1.5, 'color','k');
    scatter(PlantSegData_R2.CentroidsX_basal,PlantSegData_R2.CentroidsY_basal,50,[0 0.4470 0.7410],'filled')
    scatter(PlantSegData_R2.CentroidsX_lum,PlantSegData_R2.CentroidsY_lum,50,[0.8500 0.3250 0.0980],'filled')
    
    
    figure;
    s1_R2 = subplot(1,2,1);
    imshow(1-PlantSegData_R2.Bound_image);
    hold on
    p_R2 = plot(X_R2,Y_R2,'linewidth',3);
    
    %// modified jet-colormap
    ColorMap_R2 = uint8(jet(length(X_R2(:,1)))*255);
    Indices4color_R2 = linspace(clims(1),clims(2),length(X_R2(:,1)));
    colourIndices_R2 = zeros(length(X_R2(:,1)),3);
    
    for i = 1:length(X_R2(:,1))
        
        [TempMin,TempInd] = min(abs(Indices4color_R2 -C_R2(i) ));
        colourIndices_R2(i,:) = ColorMap_R2(TempInd,:);
    end
    
    cd_R2 = [colourIndices_R2 uint8(ones(length(X_R2(:,1)),1))].'; %'
    
    drawnow
    set(p_R2.Edge, 'ColorBinding','interpolated', 'ColorData',cd_R2)
    
    drawnow
    p3 = plot([PlantSegData_R2.plotlinesX1;PlantSegData_R2.plotlinesX2 ],...
        [PlantSegData_R2.plotlinesY1;PlantSegData_R2.plotlinesY2],'linewidth',1.5, 'color','k');
    scatter(PlantSegData_R2.CentroidsX_basal,PlantSegData_R2.CentroidsY_basal,50,[0 0.4470 0.7410],'filled')
    scatter(PlantSegData_R2.CentroidsX_lum,PlantSegData_R2.CentroidsY_lum,50,[0.8500 0.3250 0.0980],'filled')
    
    
    plot((773:1000)-500, 800.*ones(1,228),'-k')
    hold off
    title("R2  Cleared")
    
    %
    s2_R2 = subplot(1,2,2);
    imshow(Icurv_R2)
    hold on
    surf([X_R2(:) X_R2(:)], [Y_R2(:) Y_R2(:)], [Z_R2(:) Z_R2(:)], [C_R2 C_R2], ...  % Reshape and replicate data
        'FaceColor', 'none', ...    % Don't bother filling faces with color
        'EdgeColor', 'interp', ...  % Use interpolated color for edges
        'LineWidth', 5);            % Make a thicker line
    hold off
    cmap = jet;
    colormap(s2_R2,cmap);
    cb_R2 = colorbar;  % Add a colorbar
    cb_R2.Label.String = 'Curvature ($\mu m^{-1}$)';
    cb_R2.Label.Interpreter = 'latex';
    cb_R2.AxisLocation = 'in';
    cbPos_R2 = get(cb_R2,'Position');
    cb_R2.TickLabelInterpreter = 'latex';
    cb_R2.Position = [cbPos_R2(1)-0.37 cbPos_R2(2)*0.9 cbPos_R2(3).*0.7 cbPos_R2(4)*1.2];
    s1Pos_R2 = get(s1_R2,'position');
    s2Pos_R2 = get(s2_R2,'position');
    s2Pos_R2(3:4) = [s1Pos_R2(3:4)];
    set(s2_R2,'position',s2Pos_R2);
    caxis([clims(1),clims(2)])
    drawnow
    
    %set(get(s2_R2,'Children'),'Visible','off');
    
    if PlotHistograms == 1
        figure;
        histogram(rmoutliers(1./C_R2).*XSCALE_R2_cleared,'normalization','probability');
        ylabel("Probability")
        xlabel("Principle radius ($\mu m$)")
        
        xline(25,'--','linewidth',2)
        
        title("R2  Cleared")
        
    end
    
    
    
    % connectivity analysis
    AdjacencyMat_R2 = ConstructAdjacencyMatrix(PlantSegData_R2.CellData,w1,w2,0);
    
[ BasalCellHomoCount_R2,...
    BasalCellHeteroCount_R2,...
    LuminalCellHomoCount_R2,...
    LuminalCellHeteroCount_R2,...
    RemoveCells_R2] = TypeDependentConnectivtityCounter(PlantSegData_R2.CellData,AdjacencyMat_R2,w1,w2);

 %remove non existent cells
    AdjacencyMat_R2(RemoveCells_R2,:) = [];
    AdjacencyMat_R2(:,RemoveCells_R2) = [];
    PlantSegData_R2.CellData(RemoveCells_R2,:) = [];

    
    for i = 1:length(PlantSegData_R2.CellData(:,1))
        
        
        if PlantSegData_R2.CellData(i,1).cellType == "Luminal"
            R2_luminal_count = R2_luminal_count + 1;
        else
            R2_basal_count =  R2_basal_count + 1;
        end
        
        if PlantSegData_R2.CellData(i,:).Area > 0
            
            if PlantSegData_R2.CellData(i,1).cellType == "Luminal"
                LumCellArea = [LumCellArea, PlantSegData_R2.CellData(i,:).Area*XSCALE_R2_cleared.^2];
                
            elseif  PlantSegData_R2.CellData(i,1).cellType == "Basal"
                BasCellArea = [BasCellArea, PlantSegData_R2.CellData(i,:).Area*XSCALE_R2_cleared.^2];
                
            end
        end
        
    end

    
%     figure;
%     subplot(2,2,1)
%     histogram(BasalCellHomoCount_R2,'Normalization','probability')
%     hold on
%     ylabel("Probability")
%     xlabel("Number of Connections")
%     title('Homotypic connections')
%     legend('Basal')
%     
%     
%     subplot(2,2,3)
%     histogram(LuminalCellHomoCount_R2,'Normalization','probability','facecolor',[0.8500 0.3250 0.0980])
%     hold on
%     legend('Luminal')
%     ylabel("Probability")
%     xlabel("Number of Connections")
%     
%     subplot(2,2,2)
%     histogram(BasalCellHeteroCount_R2,'Normalization','probability')
%     ylabel("Probability")
%     xlabel("Number of Connections")
%     legend('Basal')
%     title('Heteroypic connections')
%     
%     subplot(2,2,4)
%     histogram(LuminalCellHeteroCount_R2,'Normalization','probability','facecolor',[0.8500 0.3250 0.0980])
%     legend('Luminal')
%     
%     ylabel("Probability")
%     xlabel("Number of Connections")
    
    
    %local curvature analysis
    
    uniqueCellsToMeanCurvatures_R2 = ClosestCellCurvatureCount(PlantSegData_R2.CellData,...
        shape_details_R2,AdjacencyMat_R2,w1,w2);

    
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% R3_1 CLEARED

if plotR3_1cleared == 1
    
    PlantSegData_R3_1 = load('PlantSeg_pipeline_outputs/R3_1_cleared/Cleared_R2_comp_membrane009_PS_pipeOutput.mat');
   
   
   
    
    twoDSliceR3_1 = imread("../PlantSeg/Cleared_R3_1/Cleared_R3_1_membrane_composite0009.tif");
    BinImageR3_1 = imbinarize(twoDSliceR3_1,0.1);
    FilterImageR3_1 = medfilt2(BinImageR3_1);
    FillImageR3_1 = imfill(FilterImageR3_1,'holes');
    SmoothImR3_1 = bwareaopen(FillImageR3_1,50);
    N = 9;
    kernel = ones(N, N) / N^3;
    blurryImageR3_1 = conv2(double(SmoothImR3_1), kernel, 'same');
    newBinaryImageR3_1 = blurryImageR3_1>0.1 ;
    
    oragnoidArea = [oragnoidArea,bwarea(newBinaryImageR3_1).*XSCALE_R3_1_cleared.^2];
    %check the input image
    % figure;
    % subplot(1,2,1)
    % imshow(twoDSliceR3_1)
    % subplot(1,2,2)
    % imshow(newBinaryImageR3_1)
    % hold off
    
    boundaryPoint = round(BoundaryConst/XSCALE_R3_1_cleared)      %12, 10 number of boundary points curvature is found over
    curvatureThresh = 0.06;     %0.06 the maximum allowed value of the curvature measure
    bp_tangent = 10;            % number of boundary points the tangent angle is found over
    interpdmin = 0.3;           % the minimum number of pixels seperating boundary points after interpolation
    loopclose = 1;   % 0 - if open boundaries | 1 - if closed boundaries
    
    % Find the curvature
    [shape_details_R3_1, Icurv_R3_1] = curvature2D(newBinaryImageR3_1', boundaryPoint, curvatureThresh, bp_tangent, ...
        interpdmin, loopclose);
    shape_details_R3_1.curvature = shape_details_R3_1.curvature.*(1/XSCALE_R3_1_cleared);
    %Plot curvature
    
    X_R3_1 = shape_details_R3_1.XY(:,1);
    Y_R3_1 = shape_details_R3_1.XY(:,2);
    Z_R3_1 = zeros(size(X_R3_1));
    C_R3_1 = shape_details_R3_1.curvature';
    % Plot
    
    figure;
    s1_R3_1 = subplot(1,2,1);
    imshow(imcomplement(PlantSegData_R3_1.Bound_image));
    hold on
    p_R3_1 = plot(X_R3_1,Y_R3_1,'linewidth',3);
    
    %// modified jet-colormap
    ColorMap_R3_1 = uint8(jet(length(X_R3_1(:,1)))*255);
    Indices4color_R3_1 = linspace(clims(1),clims(2),length(X_R3_1(:,1)));
    colourIndices_R3_1 = zeros(length(X_R3_1(:,1)),3);
    
    for i = 1:length(X_R3_1(:,1))
        
        [TempMin,TempInd] = min(abs(Indices4color_R3_1 -C_R3_1(i) ));
        colourIndices_R3_1(i,:) = ColorMap_R3_1(TempInd,:);
    end
    
    cd_R3_1 = [colourIndices_R3_1 uint8(ones(length(X_R3_1(:,1)),1))].'; %'
    
    drawnow
    set(p_R3_1.Edge, 'ColorBinding','interpolated', 'ColorData',cd_R3_1)
    
    p3 = plot([PlantSegData_R3_1.plotlinesX1;PlantSegData_R3_1.plotlinesX2 ],[PlantSegData_R3_1.plotlinesY1;PlantSegData_R3_1.plotlinesY2],'linewidth',1.5, 'color','k');
    scatter(PlantSegData_R3_1.CentroidsX_basal,PlantSegData_R3_1.CentroidsY_basal,50,[0 0.4470 0.7410],'filled')
    scatter(PlantSegData_R3_1.CentroidsX_lum,PlantSegData_R3_1.CentroidsY_lum,50,[0.8500 0.3250 0.0980],'filled')
      plot((773:1000)-300, 500.*ones(1,228),'-k')
    
    hold off
    title("R3 1 Cleared")
    
    s2_R3_1 = subplot(1,2,2);
    imshow(Icurv_R3_1)
    hold on
    surf([X_R3_1(:) X_R3_1(:)], [Y_R3_1(:) Y_R3_1(:)], [Z_R3_1(:) Z_R3_1(:)], [C_R3_1 C_R3_1], ...  % Reshape and replicate data
        'FaceColor', 'none', ...    % Don't bother filling faces with color
        'EdgeColor', 'interp', ...  % Use interpolated color for edges
        'LineWidth', 3);            % Make a thicker line
    hold off
    cmap = jet;
    colormap(s2_R3_1,cmap);
    cb_R3_1 = colorbar;  % Add a colorbar
    cb_R3_1.Label.String = 'Curvature ($\mu m^{-1}$)';
    cb_R3_1.Label.Interpreter = 'latex';
    cb_R3_1.AxisLocation = 'in';
    cbPos_R3_1 = get(cb_R3_1,'Position');
    cb_R3_1.TickLabelInterpreter = 'latex';
    cb_R3_1.Position = [cbPos_R3_1(1)-0.37 cbPos_R3_1(2)*0.9 cbPos_R3_1(3).*0.7 cbPos_R3_1(4)*1.2];
    s1Pos_R3_1 = get(s1_R3_1,'position');
    s2Pos_R3_1 = get(s2_R3_1,'position');
    s2Pos_R3_1(3:4) = [s1Pos_R3_1(3:4)];
    set(s2_R3_1,'position',s2Pos_R3_1);
    caxis([clims(1),clims(2)])
    drawnow
    set(get(s2_R3_1,'Children'),'Visible','off');
    
    if PlotHistograms == 1
        figure;
        histogram(rmoutliers(1./C_R3_1).*XSCALE_R3_1_cleared,'normalization','probability');
        ylabel("Probability")
        xlabel("Principle radius ($\mu m$)")
        
        xline(25,'--','linewidth',2)
        title("R3 1 Cleared")
        
    end
    
    
    
    % connectivity analysis
    AdjacencyMat_R3_1 = ConstructAdjacencyMatrix(PlantSegData_R3_1.CellData,w1,w2,0);
    
[ BasalCellHomoCount_R3_1,...
    BasalCellHeteroCount_R3_1,...
    LuminalCellHomoCount_R3_1,...
    LuminalCellHeteroCount_R3_1,...
    RemoveCells_R3_1] = TypeDependentConnectivtityCounter(PlantSegData_R3_1.CellData,AdjacencyMat_R3_1,w1,w2);


    %remove non existent cells
    AdjacencyMat_R3_1(RemoveCells_R3_1,:) = [];
    AdjacencyMat_R3_1(:,RemoveCells_R3_1) = [];
    PlantSegData_R3_1.CellData(RemoveCells_R3_1,:) = [];
    
    
    for i = 1:length(PlantSegData_R3_1.CellData(:,1))
        
        if PlantSegData_R3_1.CellData(i,:).cellType == "Luminal"
            R3_1_luminal_count = R3_1_luminal_count + 1;
        else
            R3_1_basal_count = R3_1_basal_count + 1;
        end
        
        
        if PlantSegData_R3_1.CellData(i,:).Area > 0
            if PlantSegData_R3_1.CellData(i,1).cellType == "Luminal"
                LumCellArea = [LumCellArea, PlantSegData_R3_1.CellData(i,:).Area*XSCALE_R3_1_cleared.^2];
                
            elseif  PlantSegData_R3_1.CellData(i,1).cellType == "Basal"
                BasCellArea = [BasCellArea, PlantSegData_R3_1.CellData(i,:).Area*XSCALE_R3_1_cleared.^2];
                
            end
        end
        
    end
    
%     figure;
%     subplot(2,2,1)
%     histogram(BasalCellHomoCount_R3_1,'Normalization','probability')
%     hold on
%     ylabel("Probability")
%     xlabel("Number of Connections")
%     title('Homotypic connections')
%     legend('Basal')
%     
%     
%     subplot(2,2,3)
%     histogram(LuminalCellHomoCount_R3_1,'Normalization','probability','facecolor',[0.8500 0.3250 0.0980])
%     hold on
%     legend('Luminal')
%     ylabel("Probability")
%     xlabel("Number of Connections")
%     
%     subplot(2,2,2)
%     histogram(BasalCellHeteroCount_R3_1,'Normalization','probability')
%     ylabel("Probability")
%     xlabel("Number of Connections")
%     legend('Basal')
%     title('Heteroypic connections')
%     
%     subplot(2,2,4)
%     histogram(LuminalCellHeteroCount_R3_1,'Normalization','probability','facecolor',[0.8500 0.3250 0.0980])
%     legend('Luminal')
%     
%     ylabel("Probability")
%     xlabel("Number of Connections")
%     
    
    
    
      %local curvature analysis
    
    uniqueCellsToMeanCurvatures_R3_1 = ClosestCellCurvatureCount(PlantSegData_R3_1.CellData,...
        shape_details_R3_1,AdjacencyMat_R3_1,w1,w2);

    
    
    
end

%%%%%%%%%%%%%%%%%%% END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 2018 dapi tom

if plot2018Dapicleared == 1
    
    PlantSegData_2018 = load('PlantSeg_pipeline_outputs/N1_GFP_2018_small_cleared/20180615_S7_DAPI GFP TOM cc3_Cy5_combined0011_edited2_PS_pipeOutput.mat');
    

    
    twoDSlice_2018 = imread("../PlantSeg/2018_dapi_tom/20180615_S7_DAPI GFP TOM cc3_Cy5_combined0011_SCALED.tif");
    twoDSlice_2018 = twoDSlice_2018';
    BinImage_2018 = imbinarize(twoDSlice_2018,0.1);
    FilterImage_2018 = medfilt2(BinImage_2018);
    FillImage_2018 = imfill(FilterImage_2018,'holes');
    SmoothIm_2018 = bwareaopen(FillImage_2018,50);
    N = 9;
    kernel = ones(N, N) / N^3;
    blurryImage_2018 = conv2(double(SmoothIm_2018), kernel, 'same');
    newBinaryImage_2018 = blurryImage_2018>0.1 ;
    
    oragnoidArea = [oragnoidArea,bwarea(newBinaryImage_2018).*XSCALE_dapi2018_scaled.^2];
    %check the input image
    % figure;
    % subplot(1,2,1)
    % imshow(twoDSliceR3_1)
    % subplot(1,2,2)
    % imshow(newBinaryImageR3_1)
    % hold off
    
    boundaryPoint = round(BoundaryConst/XSCALE_dapi2018_scaled)         %12, 10 number of boundary points curvature is found over
    curvatureThresh = 0.06;     %0.06 the maximum allowed value of the curvature measure
    bp_tangent = 10;            % number of boundary points the tangent angle is found over
    interpdmin = 0.3;           % the minimum number of pixels seperating boundary points after interpolation
    loopclose = 1;   % 0 - if open boundaries | 1 - if closed boundaries
    
    % Find the curvature
    [shape_details_2018, Icurv_2018] = curvature2D(newBinaryImage_2018, boundaryPoint, curvatureThresh, bp_tangent, ...
        interpdmin, loopclose);
   
    shape_details_2018.curvature =shape_details_2018.curvature'*1.*(1/XSCALE_dapi2018_scaled);
    
    %Plot curvature
    X_2018 = shape_details_2018.XY(:,1);
    Y_2018 = shape_details_2018.XY(:,2);
    Z_2018 = zeros(size(X_2018));
    C_2018 = shape_details_2018.curvature'*1;
    % Plot

    figure;
    s1_2018 = subplot(1,2,1);
    imshow(imcomplement(PlantSegData_2018.Bound_image)');
    hold on
    p_2018 = plot(X_2018,Y_2018,'linewidth',3);
    
    %// modified jet-colormap
    ColorMap_2018 = uint8(jet(length(X_2018(:,1)))*255);
    Indices4color_2018 = linspace(clims(1),clims(2),length(X_2018(:,1)));
    colourIndices_2018 = zeros(length(X_2018(:,1)),3);
    
    for i = 1:length(X_2018(:,1))
        
        [TempMin,TempInd] = min(abs(Indices4color_2018 -C_2018(i) ));
        colourIndices_2018(i,:) = ColorMap_2018(TempInd,:);
    end
    
    cd_2018 = [colourIndices_2018 uint8(ones(length(X_2018(:,1)),1))].'; %'
    
    drawnow
    set(p_2018.Edge, 'ColorBinding','interpolated', 'ColorData',cd_2018)
    
       
    p3 = plot([PlantSegData_2018.plotlinesX1;PlantSegData_2018.plotlinesX2 ],...
        [PlantSegData_2018.plotlinesY1;PlantSegData_2018.plotlinesY2],'linewidth',1.5, 'color','k');
    scatter(PlantSegData_2018.CentroidsX_basal,PlantSegData_2018.CentroidsY_basal,50,[0 0.4470 0.7410],'filled')
    scatter(PlantSegData_2018.CentroidsX_lum,PlantSegData_2018.CentroidsY_lum,50,[0.8500 0.3250 0.0980],'filled')
    
    plot((778:1000)-300, 700.*ones(1,223),'-k')
    
    hold off
    title("Dapi Tom 2018 Cleared")
    
    s2_2018 = subplot(1,2,2);
    imshow(Icurv_2018)
    hold on
    surf([X_2018(:) X_2018(:)], [Y_2018(:) Y_2018(:)], [Z_2018(:) Z_2018(:)], [C_2018 C_2018], ...  % Reshape and replicate data
        'FaceColor', 'none', ...    % Don't bother filling faces with color
        'EdgeColor', 'interp', ...  % Use interpolated color for edges
        'LineWidth', 3);            % Make a thicker line
    hold off
    cmap = jet;
    colormap(s2_2018,cmap);
    cb_2018 = colorbar;  % Add a colorbar
    cb_2018.Label.String = 'Curvature ($\mu m^{-1}$)';
    cb_2018.Label.Interpreter = 'latex';
    cb_2018.AxisLocation = 'in';
    cbPos_2018 = get(cb_2018,'Position');
    cb_2018.TickLabelInterpreter = 'latex';
    cb_2018.Position = [cbPos_2018(1)-0.37 cbPos_2018(2)*0.9 cbPos_2018(3).*0.7 cbPos_2018(4)*1.2];
    s1Pos_2018 = get(s1_2018,'position');
    s2Pos_2018 = get(s2_2018,'position');
    s2Pos_2018(3:4) = [s1Pos_2018(3:4)];
    set(s2_2018,'position',s2Pos_2018);
    caxis([clims(1),clims(2)])
    drawnow
    
    
    set(get(s2_2018,'Children'),'Visible','off');
    
    
    if PlotHistograms == 1
        figure;
        histogram(abs(rmoutliers(1./C_2018)).*XSCALE_dapi2018_scaled,'normalization','probability');
        ylabel("Probability")
        xlabel("Principle radius ($\mu m$)")
        
        xline(25,'--','linewidth',2)
        title("Dapi Tom 2018")
        
    end
    
    
       % connectivity analysis
    AdjacencyMat_2018 = ConstructAdjacencyMatrix(PlantSegData_2018.CellData,w1,w2,0);
    
[ BasalCellHomoCount_2018,...
    BasalCellHeteroCount_2018,...
    LuminalCellHomoCount_2018,...
    LuminalCellHeteroCount_2018,...
    RemoveCells_2018] = TypeDependentConnectivtityCounter(PlantSegData_2018.CellData,AdjacencyMat_2018,w1,w2);


    %remove non existent cells
    AdjacencyMat_2018(RemoveCells_2018,:) = [];
    AdjacencyMat_2018(:,RemoveCells_2018) = [];
    PlantSegData_2018.CellData(RemoveCells_2018,:) = [];
    
    
    for i = 1:length(PlantSegData_2018.CellData(:,1))
        
        if PlantSegData_2018.CellData(i,:).Area > 0
            if PlantSegData_2018.CellData(i,1).cellType == "Luminal"
                LumCellArea = [LumCellArea, PlantSegData_2018.CellData(i,:).Area*XSCALE_dapi2018_scaled.^2];
                r2018_luminal_count = r2018_luminal_count + 1;
            elseif  PlantSegData_2018.CellData(i,1).cellType == "Basal"
                BasCellArea = [BasCellArea, PlantSegData_2018.CellData(i,:).Area*XSCALE_dapi2018_scaled.^2];
                r2018_basal_count = r2018_basal_count + 1;
            end
            
        end
        
    end
  
%     figure;
%     subplot(2,2,1)
%     histogram(BasalCellHomoCount_2018,'Normalization','probability')
%     hold on
%     ylabel("Probability")
%     xlabel("Number of Connections")
%     title('Homotypic connections')
%     legend('Basal')
%     
%     
%     subplot(2,2,3)
%     histogram(LuminalCellHomoCount_2018,'Normalization','probability','facecolor',[0.8500 0.3250 0.0980])
%     hold on
%     legend('Luminal')
%     ylabel("Probability")
%     xlabel("Number of Connections")
%     
%     subplot(2,2,2)
%     histogram(BasalCellHeteroCount_2018,'Normalization','probability')
%     ylabel("Probability")
%     xlabel("Number of Connections")
%     legend('Basal')
%     title('Heteroypic connections')
%     
%     subplot(2,2,4)
%     histogram(LuminalCellHeteroCount_2018,'Normalization','probability','facecolor',[0.8500 0.3250 0.0980])
%     legend('Luminal')
%     
%     ylabel("Probability")
%     xlabel("Number of Connections")
    
     %local curvature analysis
    
    uniqueCellsToMeanCurvatures_2018 = ClosestCellCurvatureCount(PlantSegData_2018.CellData,...
        shape_details_2018,AdjacencyMat_2018,w1,w2);

    
    
end

%%%%%%%%%%%%%%%%%%%% END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% GFP 2 

if plotGFP2cleared == 1
    
    
    
    PlantSegData_GFP2 = load('PlantSeg_pipeline_outputs/N1_TOM_GFP_2/C1-DAPI_GFP_TOM_2-1-2_adj_predictions_gasp_average_19_edited2_PS_pipeOutput.mat');
    

    twoDSlice_GFP2 = imread("../PlantSeg/DAPI_GFP_2_N1/C1-DAPI_GFP_TOM_2-1-2_adj.tif");
    twoDSlice_GFP2 = twoDSlice_GFP2';
    BinImage_GFP2 = imbinarize(twoDSlice_GFP2,0.1);
    FilterImage_GFP2 = medfilt2(BinImage_GFP2);
    FillImage_GFP2 = imfill(FilterImage_GFP2,'holes');
    SmoothIm_GFP2 = bwareaopen(FillImage_GFP2,100);
    N = 11;
    kernel = ones(N, N) / N^3;
    blurryImage_GFP2 = conv2(double(SmoothIm_GFP2), kernel, 'same');
    newBinaryImage_GFP2 = blurryImage_GFP2>0.07 ;
    newBinaryImage_GFP2 = bwareaopen(newBinaryImage_GFP2,100);
    oragnoidArea = [oragnoidArea,bwarea(newBinaryImage_GFP2).*XSCALE_dapi_GFP2_scaled.^2];
     %check the input image
%     figure;
%     subplot(1,2,1)
%     imshow(twoDSlice_GFP2)
%     subplot(1,2,2)
%     imshow(newBinaryImage_GFP2)
%     hold off
    
    boundaryPoint = round(BoundaryConst/XSCALE_dapi_GFP2_scaled);         %12, 10 number of boundary points curvature is found over
    curvatureThresh = 0.06;     %0.06 the maximum allowed value of the curvature measure
    bp_tangent = 10;            % number of boundary points the tangent angle is found over
    interpdmin = 0.3;           % the minimum number of pixels seperating boundary points after interpolation
    loopclose = 1;   % 0 - if open boundaries | 1 - if closed boundaries
    
    % Find the curvature
    [shape_details_GFP2, Icurv_GFP2] = curvature2D(newBinaryImage_GFP2, boundaryPoint, curvatureThresh, bp_tangent, ...
        interpdmin, loopclose);
    
    shape_details_GFP2.curvature = shape_details_GFP2.curvature.*(1/XSCALE_dapi_GFP2_scaled);
    
    %Plot curvature
    X_GFP2 = shape_details_GFP2.XY(:,1);
    Y_GFP2 = shape_details_GFP2.XY(:,2);
    Z_GFP2 = zeros(size(X_GFP2));
    C_GFP2 = shape_details_GFP2.curvature'*1;
    % Plot

    figure;
    s1_GFP2 = subplot(1,2,1);
    imshow(imcomplement(PlantSegData_GFP2.Bound_image)');
    hold on
    p_GFP2 = plot(X_GFP2,Y_GFP2,'linewidth',3);
    
    %// modified jet-colormap
    ColorMap_GFP2 = uint8(jet(length(X_GFP2(:,1)))*255);
    Indices4color_GFP2 = linspace(clims(1),clims(2),length(X_GFP2(:,1)));
    colourIndices_GFP2 = zeros(length(X_GFP2(:,1)),3);
    
    for i = 1:length(X_GFP2(:,1))
        
        [TempMin,TempInd] = min(abs(Indices4color_GFP2 -C_GFP2(i) ));
        colourIndices_GFP2(i,:) = ColorMap_GFP2(TempInd,:);
    end
    
    cd_GFP2 = [colourIndices_GFP2 uint8(ones(length(X_GFP2(:,1)),1))].'; %'
    
    drawnow
    set(p_GFP2.Edge, 'ColorBinding','interpolated', 'ColorData',cd_GFP2)
    
       
    p3 = plot([PlantSegData_GFP2.plotlinesX1;PlantSegData_GFP2.plotlinesX2 ],...
        [PlantSegData_GFP2.plotlinesY1;PlantSegData_GFP2.plotlinesY2],'linewidth',1.5, 'color','k');
    scatter(PlantSegData_GFP2.CentroidsX_basal,PlantSegData_GFP2.CentroidsY_basal,50,[0 0.4470 0.7410],'filled')
    scatter(PlantSegData_GFP2.CentroidsX_lum,PlantSegData_GFP2.CentroidsY_lum,50,[0.8500 0.3250 0.0980],'filled')
    
    
    hold off
    title("Dapi Tom GFP2 Cleared")
    
    s2_GFP2 = subplot(1,2,2);
    imshow(Icurv_GFP2)
    hold on
    surf([X_GFP2(:) X_GFP2(:)], [Y_GFP2(:) Y_GFP2(:)], [Z_GFP2(:) Z_GFP2(:)], [C_GFP2 C_GFP2], ...  % Reshape and replicate data
        'FaceColor', 'none', ...    % Don't bother filling faces with color
        'EdgeColor', 'interp', ...  % Use interpolated color for edges
        'LineWidth', 3);            % Make a thicker line
    hold off
    cmap = jet;
    colormap(s2_GFP2,cmap);
    cb_GFP2 = colorbar;  % Add a colorbar
    cb_GFP2.Label.String = 'Curvature ($\mu m^{-1}$)';
    cb_GFP2.Label.Interpreter = 'latex';
    cb_GFP2.AxisLocation = 'in';
    cbPos_GFP2 = get(cb_GFP2,'Position');
    cb_GFP2.TickLabelInterpreter = 'latex';
    cb_GFP2.Position = [cbPos_GFP2(1)-0.37 cbPos_GFP2(2)*0.9 cbPos_GFP2(3).*0.7 cbPos_GFP2(4)*1.2];
    s1Pos_GFP2 = get(s1_GFP2,'position');
    s2Pos_GFP2 = get(s2_GFP2,'position');
    s2Pos_GFP2(3:4) = [s1Pos_GFP2(3:4)];
    set(s2_GFP2,'position',s2Pos_GFP2);
    caxis([clims(1),clims(2)])
    drawnow
    
    
    set(get(s2_GFP2,'Children'),'Visible','off');
    
    
         % connectivity analysis
    AdjacencyMat_GFP2 = ConstructAdjacencyMatrix(PlantSegData_GFP2.CellData,w1,w2,0);
    
[ BasalCellHomoCount_GFP2,...
    BasalCellHeteroCount_GFP2,...
    LuminalCellHomoCount_GFP2,...
    LuminalCellHeteroCount_GFP2,...
    RemoveCells_GFP2] = TypeDependentConnectivtityCounter(PlantSegData_GFP2.CellData,AdjacencyMat_GFP2,w1,w2);


    %remove non existent cells
    AdjacencyMat_GFP2(RemoveCells_GFP2,:) = [];
    AdjacencyMat_GFP2(:,RemoveCells_GFP2) = [];
    PlantSegData_GFP2.CellData(RemoveCells_GFP2,:) = [];
    
    
    for i = 1:length(PlantSegData_GFP2.CellData(:,1))
        
        
        
        if PlantSegData_GFP2.CellData(i,:).Area > 0
            if PlantSegData_GFP2.CellData(i,1).cellType == "Luminal"
                LumCellArea = [LumCellArea, PlantSegData_GFP2.CellData(i,:).Area*XSCALE_dapi_GFP2_scaled.^2];
                GFP2_luminal_count = GFP2_luminal_count + 1;
            elseif  PlantSegData_GFP2.CellData(i,1).cellType == "Basal"
                BasCellArea = [BasCellArea, PlantSegData_GFP2.CellData(i,:).Area*XSCALE_dapi_GFP2_scaled.^2];
                GFP2_basal_count = GFP2_basal_count + 1;
            end
            
        end
        
    end
%   
%     figure;
%     subplot(2,2,1)
%     histogram(BasalCellHomoCount_GFP2,'Normalization','probability')
%     hold on
%     ylabel("Probability")
%     xlabel("Number of Connections")
%     title('Homotypic connections')
%     legend('Basal')
%     
%     
%     subplot(2,2,3)
%     histogram(LuminalCellHomoCount_GFP2,'Normalization','probability','facecolor',[0.8500 0.3250 0.0980])
%     hold on
%     legend('Luminal')
%     ylabel("Probability")
%     xlabel("Number of Connections")
%     
%     subplot(2,2,2)
%     histogram(BasalCellHeteroCount_GFP2,'Normalization','probability')
%     ylabel("Probability")
%     xlabel("Number of Connections")
%     legend('Basal')
%     title('Heteroypic connections')
%     
%     subplot(2,2,4)
%     histogram(LuminalCellHeteroCount_GFP2,'Normalization','probability','facecolor',[0.8500 0.3250 0.0980])
%     legend('Luminal')
%     
%     ylabel("Probability")
%     xlabel("Number of Connections")
    
     %local curvature analysis
    
    uniqueCellsToMeanCurvatures_GFP2 = ClosestCellCurvatureCount(PlantSegData_GFP2.CellData,...
        shape_details_GFP2,AdjacencyMat_GFP2,w1,w2);
    
    
    
end




%% Large N1

if plotLargeN1cleared == 1
    
    
    
    PlantSegData_LN1 = load('PlantSeg_pipeline_outputs/LargeN1/LargeN1Result_0009_scaled_predictions_gasp_average_edited_PS_pipeOutput.mat');
    

    
        
    twoDSlice_LN1 = imread("../PlantSeg/LargeN1_2018_EDITED/Result_0009_scaled.tif");
    twoDSlice_LN1 = twoDSlice_LN1';
    BinImage_LN1 = imbinarize(twoDSlice_LN1,0.1);
    FilterImage_LN1 = medfilt2(BinImage_LN1);
    FillImage_LN1 = imfill(FilterImage_LN1,'holes');
    SmoothIm_LN1 = bwareaopen(FillImage_LN1,100);
    N = 13;
    kernel = ones(N, N) / N^3;
    blurryImage_LN1 = conv2(double(SmoothIm_LN1), kernel, 'same');
    newBinaryImage_LN1 = blurryImage_LN1>0.07 ;
    newBinaryImage_LN1 = bwareaopen(newBinaryImage_LN1,100);
    oragnoidArea = [oragnoidArea,bwarea(newBinaryImage_LN1).*XSCALE_LargeN1_scaled.^2];
     %check the input image
%     figure;
%     subplot(1,2,1)
%     imshow(twoDSlice_LN1)
%     subplot(1,2,2)
%     imshow(newBinaryImage_LN1)
%     hold off
   
    boundaryPoint = round(BoundaryConst/XSCALE_LargeN1_scaled);         %12, 10 number of boundary points curvature is found over
    curvatureThresh = 0.06;     %0.06 the maximum allowed value of the curvature measure
    bp_tangent = 10;            % number of boundary points the tangent angle is found over
    interpdmin = 0.3;           % the minimum number of pixels seperating boundary points after interpolation
    loopclose = 1;   % 0 - if open boundaries | 1 - if closed boundaries
    
    % Find the curvature
    [shape_details_LN1, Icurv_LN1] = curvature2D(newBinaryImage_LN1, boundaryPoint, curvatureThresh, bp_tangent, ...
        interpdmin, loopclose);
    
    shape_details_LN1.curvature = shape_details_LN1.curvature.*(1/XSCALE_LargeN1_scaled);
    
    %Plot curvature
    X_LN1 = shape_details_LN1.XY(:,1);
    Y_LN1 = shape_details_LN1.XY(:,2);
    Z_LN1 = zeros(size(X_LN1));
    C_LN1 = shape_details_LN1.curvature'*1;
    % Plot

    figure;
    s1_LN1 = subplot(1,2,1);
    imshow(imcomplement(PlantSegData_LN1.Bound_image)');
    hold on
    p_LN1 = plot(X_LN1,Y_LN1,'linewidth',3);
    
    %// modified jet-colormap
    ColorMap_LN1 = uint8(jet(length(X_LN1(:,1)))*255);
    Indices4color_LN1 = linspace(clims(1),clims(2),length(X_LN1(:,1)));
    colourIndices_LN1 = zeros(length(X_LN1(:,1)),3);
    
    for i = 1:length(X_LN1(:,1))
        
        [TempMin,TempInd] = min(abs(Indices4color_LN1 -C_LN1(i) ));
        colourIndices_LN1(i,:) = ColorMap_LN1(TempInd,:);
    end
    
    cd_LN1 = [colourIndices_LN1 uint8(ones(length(X_LN1(:,1)),1))].'; %'
    
    drawnow
    set(p_LN1.Edge, 'ColorBinding','interpolated', 'ColorData',cd_LN1)
    drawnow
       plot((778:1000)+200, 1800.*ones(1,223),'-k')
       
    p3 = plot([PlantSegData_LN1.plotlinesX1;PlantSegData_LN1.plotlinesX2 ],...
        [PlantSegData_LN1.plotlinesY1;PlantSegData_LN1.plotlinesY2],'linewidth',1.5, 'color','k');
    scatter(PlantSegData_LN1.CentroidsX_basal,PlantSegData_LN1.CentroidsY_basal,50,[0 0.4470 0.7410],'filled')
    scatter(PlantSegData_LN1.CentroidsX_lum,PlantSegData_LN1.CentroidsY_lum,50,[0.8500 0.3250 0.0980],'filled')
    
    
    hold off
    title("Dapi Tom LN1 Cleared")
    
    s2_LN1 = subplot(1,2,2);
    imshow(Icurv_LN1)
    hold on
    surf([X_LN1(:) X_LN1(:)], [Y_LN1(:) Y_LN1(:)], [Z_LN1(:) Z_LN1(:)], [C_LN1 C_LN1], ...  % Reshape and replicate data
        'FaceColor', 'none', ...    % Don't bother filling faces with color
        'EdgeColor', 'interp', ...  % Use interpolated color for edges
        'LineWidth', 3);            % Make a thicker line
    hold off
    cmap = jet;
    colormap(s2_LN1,cmap);
    cb_LN1 = colorbar;  % Add a colorbar
    cb_LN1.Label.String = 'Curvature ($\mu m^{-1}$)';
    cb_LN1.Label.Interpreter = 'latex';
    cb_LN1.AxisLocation = 'in';
    cbPos_LN1 = get(cb_LN1,'Position');
    cb_LN1.TickLabelInterpreter = 'latex';
    cb_LN1.Position = [cbPos_LN1(1)-0.37 cbPos_LN1(2)*0.9 cbPos_LN1(3).*0.7 cbPos_LN1(4)*1.2];
    s1Pos_LN1 = get(s1_LN1,'position');
    s2Pos_LN1 = get(s2_LN1,'position');
    s2Pos_LN1(3:4) = [s1Pos_LN1(3:4)];
    set(s2_LN1,'position',s2Pos_LN1);
    caxis([clims(1),clims(2)])
    drawnow
    
    
    set(get(s2_LN1,'Children'),'Visible','off');
    
    
         % connectivity analysis
    AdjacencyMat_LN1 = ConstructAdjacencyMatrix(PlantSegData_LN1.CellData,w1,w2,0);
    
[ BasalCellHomoCount_LN1,...
    BasalCellHeteroCount_LN1,...
    LuminalCellHomoCount_LN1,...
    LuminalCellHeteroCount_LN1,...
    RemoveCells_LN1] = TypeDependentConnectivtityCounter(PlantSegData_LN1.CellData,AdjacencyMat_LN1,w1,w2);


    %remove non existent cells
    AdjacencyMat_LN1(RemoveCells_LN1,:) = [];
    AdjacencyMat_LN1(:,RemoveCells_LN1) = [];
    PlantSegData_LN1.CellData(RemoveCells_LN1,:) = [];
    
    for i = 1:length(PlantSegData_LN1.CellData(:,1))
        
        if PlantSegData_LN1.CellData(i,:).Area > 0
            if PlantSegData_LN1.CellData(i,1).cellType == "Luminal"
                LumCellArea = [LumCellArea, PlantSegData_LN1.CellData(i,:).Area*XSCALE_LargeN1_scaled.^2];
                LN1_luminal_count = LN1_luminal_count + 1;
            elseif  PlantSegData_LN1.CellData(i,1).cellType == "Basal"
                BasCellArea = [BasCellArea, PlantSegData_LN1.CellData(i,:).Area*XSCALE_LargeN1_scaled.^2];
                LN1_basal_count = LN1_basal_count + 1;
            end
        end
    end
  
%     figure;
%     subplot(2,2,1)
%     histogram(BasalCellHomoCount_LN1,'Normalization','probability')
%     hold on
%     ylabel("Probability")
%     xlabel("Number of Connections")
%     title('Homotypic connections')
%     legend('Basal')
%     
%     
%     subplot(2,2,3)
%     histogram(LuminalCellHomoCount_LN1,'Normalization','probability','facecolor',[0.8500 0.3250 0.0980])
%     hold on
%     legend('Luminal')
%     ylabel("Probability")
%     xlabel("Number of Connections")
%     
%     subplot(2,2,2)
%     histogram(BasalCellHeteroCount_LN1,'Normalization','probability')
%     ylabel("Probability")
%     xlabel("Number of Connections")
%     legend('Basal')
%     title('Heteroypic connections')
%     
%     subplot(2,2,4)
%     histogram(LuminalCellHeteroCount_LN1,'Normalization','probability','facecolor',[0.8500 0.3250 0.0980])
%     legend('Luminal')
%     
%     ylabel("Probability")
%     xlabel("Number of Connections")
    
     %local curvature analysis
    
    uniqueCellsToMeanCurvatures_LN1 = ClosestCellCurvatureCount(PlantSegData_LN1.CellData,...
        shape_details_LN1,AdjacencyMat_LN1,w1,w2);
    
    
    
end







%% Grouped curvatature analysis

%close all

if GroupedCurvatureAnalysis == 1
 
    close all
    
    
    CurvConnectTotalExtreme = zeros(4,4);
    CurvConnectTotalMean = zeros(4,4);
    CurvConnectTotalFreq = zeros(4,4);
    CurvRatioDependence = [];
    CompleteCurvatures  = [uniqueCellsToMeanCurvatures_R2;
        uniqueCellsToMeanCurvatures_R3_1;uniqueCellsToMeanCurvatures_2018;uniqueCellsToMeanCurvatures_GFP2;
        uniqueCellsToMeanCurvatures_LN1];
    
    CompleteCurvatures  = [uniqueCellsToMeanCurvatures_LN1];
    
    
    CurvDist = struct;
    
    for i = 1:4
        for j = 1:4
            TempC = [];
            for k = 1:length(CompleteCurvatures(:,1))
                if CompleteCurvatures(k,3) == i && CompleteCurvatures(k,4) == j
                    TempC = [TempC ,CompleteCurvatures(k ,2)];
                end
            end
         
            
            
           
            
            if ~isempty(TempC)
                
                CurvConnectTotalFreq(i,j) = length(TempC);
                
                TempMax = max(TempC);
                TempMin = min(TempC);
                
                if abs(TempMax) >= abs(TempMin)
                    
                    CurvConnectTotalExtreme(i,j) = TempMax;
                else
                    CurvConnectTotalExtreme(i,j) = TempMin;
                end
            else
                CurvConnectTotalExtreme(i,j) = NaN;
                CurvConnectTotalFreq(i,j) = NaN;
            end
            
            if length(TempC)>=2
                
                
                
             CurvConnectTotalMean(i,j) = median(TempC);
             
            else
                 CurvConnectTotalMean(i,j) = NaN;
            end
                 
             
             if length(TempC)>=3 && i <= 3
                 TempC = rmoutliers(TempC);
                 stdMean = std(TempC);
                 upperCI = prctile(TempC,75);
                 lowerCI = prctile(TempC,25);
                CurvRatioDependence = [CurvRatioDependence; [i/j,CurvConnectTotalMean(i,j),upperCI,lowerCI, i ]];
                
      
             end
             
             
             if i == 2 && length(TempC)>=1
                % TempC = rmoutliers(TempC);
                 CurvDist(j).Ratio = i/j;
                 CurvDist(j).Dist = TempC;
             end
           
        end
        
    end
    
  
    figure('Renderer', 'painters', 'Position', [10 10 0.08 0.13])
    
    
     h = heatmap(CurvConnectTotalFreq);
     %h.XLabel = "Heterotypic neighbours";
     %h.YLabel = "Homotypic neighbours";
     h.NodeChildren(3).XAxis.Label.Interpreter = 'latex';
     h.NodeChildren(3).YAxis.Label.Interpreter = 'latex';
     h.NodeChildren(3).Title.Interpreter = 'latex';
     colormap('parula')
     colorbar off
     title("R")
     
 
     
     
    
    
     MaxNum = max([ length(CurvDist(1).Dist) ,length(CurvDist(2).Dist),...
      length(CurvDist(3).Dist),length(CurvDist(4).Dist) ]  );
  
   N2_1  = [ CurvDist(1).Dist' ; NaN(MaxNum-length(CurvDist(1).Dist),1)];
       N2_2  = [ CurvDist(2).Dist' ; NaN(MaxNum-length(CurvDist(2).Dist),1)];
       N2_3  = [ CurvDist(3).Dist' ; NaN(MaxNum-length(CurvDist(3).Dist),1)];
       N2_4  = [ CurvDist(4).Dist' ; NaN(MaxNum-length(CurvDist(4).Dist),1)];
       
       xx = linspace(0,5.5,100);
       
    
 figure('Renderer', 'painters', 'Position', [10 10 0.08 0.14])
    
       
  % plot(xx,-0.003731.*xx + 0.0337,'--r')
  
   violinplot([N2_1,N2_2,N2_3,N2_4],["1","2","3","4"])
   %ylabel("Average Curvature $(\mu m^{-1})$")
   %xlabel("Number of heterotypic connections, $n_{2}$")
   title("LN1")
   
   box on
   xlim([0.5,4.5])
   ylim([-0.1,0.1])
   
   
   
   
   
   %% cell areas
   BasCellArea = BasCellArea(1,BasCellArea~=0);
   LumCellArea = LumCellArea(1,LumCellArea~=0);
   
   BasCellArea_plot = [ BasCellArea';NaN(length(LumCellArea)-length(BasCellArea),1)];
   LumCellArea_plot = LumCellArea';
   
   
   mean(BasCellArea)
   mean(LumCellArea)
   
   [h,p,ci,stats] = ttest2(BasCellArea',LumCellArea')
   
   figure;
   % violinplot([BasCellArea_plot,LumCellArea_plot],["Basal","Luminal"]);
   
   
   Xdat = [1,2];
   Ydat = [mean(BasCellArea),mean(LumCellArea)]
   bar(Xdat,Ydat,'LineWidth',2,'FaceColor',0.75.*[1 1 1])
   
   
   SEM_bas = std(BasCellArea)/sqrt(length(BasCellArea));               % Standard Error
   ts_bas = tinv([0.025  0.975],length(BasCellArea)-1);      % T-Score
   CI_bas = ts_bas*SEM_bas   ;              % Confidence Intervals
   
   SEM_lum = std(LumCellArea)/sqrt(length(LumCellArea));               % Standard Error
   ts_lum = tinv([0.025  0.975],length(LumCellArea)-1);      % T-Score
   CI_lum =  ts_lum*SEM_lum    ;
   
   hold on
   er = errorbar(Xdat,Ydat,[CI_bas(1),CI_lum(1)],[CI_bas(2),CI_lum(2)],'linewidth',2);
   er.Color = [0 0 0];
   er.LineStyle = 'none';
   er.CapSize = 12;
   
   
   ylabel("Cell area $(\mu m^{2})$")
   xlabel("Cell type")
   
   
   %% organoid connectivity
   
   
     connectCountR2 = zeros(7,7);
    for i = 1:length(BasalCellHeteroCount_R2)
        if BasalCellHomoCount_R2(i) > 0 && BasalCellHeteroCount_R2(i) > 0
            connectCountR2(BasalCellHomoCount_R2(i), BasalCellHeteroCount_R2(i)) = connectCountR2(BasalCellHomoCount_R2(i), BasalCellHeteroCount_R2(i)) + 1;
        end
    end
    
    for i = 1:length(LuminalCellHeteroCount_R2)
        if LuminalCellHeteroCount_R2(i) > 0  && LuminalCellHomoCount_R2(i) > 0
            
            connectCountR2(LuminalCellHomoCount_R2(i), LuminalCellHeteroCount_R2(i)) = connectCountR2(LuminalCellHomoCount_R2(i),LuminalCellHeteroCount_R2(i)) + 1;
        end
    end
    connectCountR2(connectCountR2 == 0) = NaN;
     figure('Renderer', 'painters', 'Position', [10 10 0.08 0.14])
    
    
     h = heatmap(connectCountR2);
     h.XLabel = "Heterotypic neighbours";
     h.YLabel = "Homotypic neighbours";
     h.NodeChildren(3).XAxis.Label.Interpreter = 'latex';
     h.NodeChildren(3).YAxis.Label.Interpreter = 'latex';
     h.NodeChildren(3).Title.Interpreter = 'latex';
     colormap('jet')
     colorbar off
     
     
          connectCountR3_1 = zeros(7,7);
    for i = 1:length(BasalCellHeteroCount_R3_1)
        if BasalCellHomoCount_R3_1(i) > 0 && BasalCellHeteroCount_R3_1(i) > 0
            
            connectCountR3_1(BasalCellHomoCount_R3_1(i), BasalCellHeteroCount_R3_1(i)) = connectCountR3_1(BasalCellHomoCount_R3_1(i), BasalCellHeteroCount_R3_1(i)) + 1;
        end
    end
    
    for i = 1:length(LuminalCellHeteroCount_R3_1)
        if LuminalCellHeteroCount_R3_1(i) > 0  && LuminalCellHomoCount_R3_1(i) > 0
                LuminalCellHomoCount_R3_1(i)
            connectCountR3_1(LuminalCellHomoCount_R3_1(i), LuminalCellHeteroCount_R3_1(i)) = connectCountR3_1(LuminalCellHomoCount_R3_1(i),LuminalCellHeteroCount_R3_1(i)) + 1;
        end
    end
    connectCountR3_1(connectCountR3_1 == 0) = NaN;
     figure('Renderer', 'painters', 'Position', [10 10 0.08 0.14])
    
    
     h = heatmap(connectCountR3_1);
     h.XLabel = "Heterotypic neighbours";
     h.YLabel = "Homotypic neighbours";
     h.NodeChildren(3).XAxis.Label.Interpreter = 'latex';
     h.NodeChildren(3).YAxis.Label.Interpreter = 'latex';
     h.NodeChildren(3).Title.Interpreter = 'latex';
     colormap('jet')
     colorbar off
     
     
             connectCount2018 = zeros(8,8);
    for i = 1:length(BasalCellHeteroCount_2018)
        if BasalCellHomoCount_2018(i) > 0 && BasalCellHeteroCount_2018(i) > 0
            
            connectCount2018(BasalCellHomoCount_2018(i), BasalCellHeteroCount_2018(i)) = connectCount2018(BasalCellHomoCount_2018(i), BasalCellHeteroCount_2018(i)) + 1;
        end
    end
    
    for i = 1:length(LuminalCellHeteroCount_2018)
        if LuminalCellHeteroCount_2018(i) > 0  && LuminalCellHomoCount_2018(i) > 0
                LuminalCellHomoCount_2018(i)
            connectCount2018(LuminalCellHomoCount_2018(i), LuminalCellHeteroCount_2018(i)) = connectCount2018(LuminalCellHomoCount_2018(i),LuminalCellHeteroCount_2018(i)) + 1;
        end
    end
    connectCount2018(connectCount2018 == 0) = NaN;
     figure('Renderer', 'painters', 'Position', [10 10 0.08 0.14])
    
    
     h = heatmap(connectCount2018);
     h.XLabel = "Heterotypic neighbours";
     h.YLabel = "Homotypic neighbours";
     h.NodeChildren(3).XAxis.Label.Interpreter = 'latex';
     h.NodeChildren(3).YAxis.Label.Interpreter = 'latex';
     h.NodeChildren(3).Title.Interpreter = 'latex';
     colormap('jet')
     colorbar off
   
   
end




%% population ratio
subplot(5,1,1)
pie([R2_basal_count,R2_luminal_count])
colormap(lines(2)); %<- Set properties

subplot(5,1,2)
pie([R3_1_basal_count,R3_1_luminal_count])

subplot(5,1,3)
pie([r2018_basal_count,r2018_luminal_count])

subplot(5,1,4)
pie([GFP2_basal_count,GFP2_luminal_count])

subplot(5,1,5)
pie([LN1_basal_count,LN1_luminal_count])
