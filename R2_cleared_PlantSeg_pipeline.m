close all
clc
clear


PixelAreaTol = 1e4;
SliceOfInterest = 16;

%tune for cell identity labelling
lumActivationTol = 0.1;
basActivationTol = 0.01;
SmoothingSigma = 1.5;
manuallyEditCellTypes = 1;
ListOfManualEdits = [563.168,700.468,"Basal"];

%choose to validate the image
ValidateImageSeg = 1;

%display segmentation results
PlotAllOutputs = 1;

%weights for the adjacency matrix
AdjacencyKeneral = 4;
w1 = 0.01;
w2 = 1;

%%%%%%%%%%%% Output directory - LOOK HERE %%%%%%%%%%%%%%%%%%%%
SaveOutputs = 0;
OutputDirect = 'PlantSeg_pipeline_outputs/R2_cleared';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%input directories

if SliceOfInterest <= 10
    
    CombinedSegFile =strcat('../PlantSeg/Cleared_R2/PS_2DUNET_DS3_16/confocal_2D_unet_bce_dice_ds3x/GASP/PostProcessing/Cleared_R2_comp_membrane000',num2str(SliceOfInterest-1),'_predictions_gasp_average.h5');
    CombinedBoundFile = strcat('../PlantSeg/Cleared_R2/PS_2DUNET_DS3_16/confocal_2D_unet_bce_dice_ds3x/PostProcessing/Cleared_R2_comp_membrane000',num2str(SliceOfInterest-1),'_predictions.h5');
    
else
    CombinedSegFile =strcat('../PlantSeg/Cleared_R2/PS_2DUNET_DS3_16/confocal_2D_unet_bce_dice_ds3x/GASP/PostProcessing/Cleared_R2_comp_membrane00',num2str(SliceOfInterest-1),'_predictions_gasp_average.h5');
    CombinedBoundFile = strcat('../PlantSeg/Cleared_R2/PS_2DUNET_DS3_16/confocal_2D_unet_bce_dice_ds3x/PostProcessing/Cleared_R2_comp_membrane00',num2str(SliceOfInterest-1),'_predictions.h5');
    
end
PS_seg = hdf5info(CombinedSegFile);
PS_bound = hdf5info(CombinedBoundFile);


Raw_image = hdf5read(PS_seg.GroupHierarchy.Datasets(1));
Seg_image = hdf5read(PS_seg.GroupHierarchy.Datasets(2));
Bound_image = hdf5read(PS_bound.GroupHierarchy.Datasets(1));


%read in image data for luminal and basal identity - run a gaussian over it
FullImageStack_bas = imread('../PlantSeg/Cleared_R2/Cleared_R2_basal.tif',SliceOfInterest);
FullImageStack_lum = imread('../PlantSeg/Cleared_R2/Cleared_R2_luminal.tif',SliceOfInterest);

FullImageStack_bas = imgaussfilt(rescale(FullImageStack_bas),SmoothingSigma);
FullImageStack_lum = imgaussfilt(rescale(FullImageStack_lum),SmoothingSigma);

%Read in validation image if required 
if ValidateImageSeg == 1
    
    Val_image = imread(strcat('../PlantSeg/Cleared_R2/Cleared_labels/Label_',num2str(SliceOfInterest),'.tif'));
    
end

%% Segementation feature extraction into cell dataset

NumberOfLabels = max(Seg_image,[],'all');
CellAreas = regionprops(Seg_image,'area');
CellPerimeters = regionprops(Seg_image,'perimeter');
CellCentroid = regionprops(Seg_image,'centroid');
CellPixelList = regionprops(Seg_image,'PixelList');
CellMajorAxis = regionprops(Seg_image,'MajorAxisLength');


%filter out any cells that are too large - i.e background and lumen
RemoveLabels = [];
for i = 1:NumberOfLabels
    if CellAreas(i,1).Area > PixelAreaTol
        RemoveLabels = [RemoveLabels,i];
    end
end

%Set the removed labels to NaN
ListOfLabel = 1:NumberOfLabels;
ListOfLabel(RemoveLabels) = [];

%create Data structure for each cell
CellData = struct;
[ImageSizeX,ImageSizeY] = size(Seg_image);

for ii = 1:length(ListOfLabel)
    % cell ids
    CellData(ii,1).LabelId = ListOfLabel(ii);
    CellData(ii,1).Area = CellAreas(ListOfLabel(ii)).Area;
    CellData(ii,1).PixelLocations = CellPixelList(ListOfLabel(ii)).PixelList;
    CellData(ii,1).Perimeter = CellPerimeters(ListOfLabel(ii)).Perimeter;
    CellData(ii,1).Centroid = CellCentroid(ListOfLabel(ii)).Centroid;
end

%Determine if a cell if basal or luminal


for ii = 1:length(ListOfLabel)
    
    PixelsToCheck = CellData(ii,1).PixelLocations;
    SumOfBasalActivePixels = 0;
    SumOfLuminalActivePixels = 0;
    
    for i = 1:length(PixelsToCheck(:,1))
        
        if FullImageStack_bas(PixelsToCheck(i,1),PixelsToCheck(i,2)) > basActivationTol
            SumOfBasalActivePixels = SumOfBasalActivePixels + 1;
        end
        
        if FullImageStack_lum(PixelsToCheck(i,1),PixelsToCheck(i,2)) > lumActivationTol
            SumOfLuminalActivePixels = SumOfLuminalActivePixels + 1;
        end
        
    end
    
    if SumOfLuminalActivePixels > SumOfBasalActivePixels
        CellData(ii,1).cellType = "Luminal";
    else
        CellData(ii,1).cellType = "Basal";
    end
    
    if manuallyEditCellTypes == 1
        
        for xx = 1:length(ListOfManualEdits(:,1))
            if ~isnan(CellData(ii,1).Centroid)
                if norm(cast(ListOfManualEdits(xx,1:2),'double')- CellData(ii,1).Centroid) < 1e-2
                    
                    CellData(ii,1).cellType = ListOfManualEdits(xx,3);
                end
            end
        end
        
    end
    
    
end

%Visually check correct identities

CellIdentitiesMap = zeros(ImageSizeX,ImageSizeY,'uint8');

for ii = 1:length(ListOfLabel)
    
    PixelsToCheck = CellData(ii,1).PixelLocations;
    for i = 1:length(PixelsToCheck(:,1))
        if CellData(ii,1).cellType == "Luminal"
            
            CellIdentitiesMap(PixelsToCheck(i,1),PixelsToCheck(i,2)) = 2;
        else
            CellIdentitiesMap(PixelsToCheck(i,1),PixelsToCheck(i,2)) = 1;
        end
    end
    
end


%% Validation using pre-defined labels

if ValidateImageSeg == 1
    %Jaccard index for cell identities
    InterCount = 0;
    UnionCount = 0;
    
    InterCount_identities = 0;
    UnionCount_identities = 0;
    
    for i = 1:ImageSizeX
        for j = 1:ImageSizeY
            
            % for generic pxiels
            if CellIdentitiesMap(i,j) > 0 || Val_image(i,j) >0
                UnionCount = UnionCount + 1;
                
                if CellIdentitiesMap(i,j) > 0 && Val_image(i,j) >0
                    InterCount = InterCount + 1;
                    
                end
                
            end
            
            
            if CellIdentitiesMap(i,j) == 1 || Val_image(i,j)  == 1
                UnionCount_identities = UnionCount_identities + 1;
                
                if CellIdentitiesMap(i,j) == 1 && Val_image(i,j) == 1
                    InterCount_identities = InterCount_identities + 1;
                    
                end
            end
            
            if CellIdentitiesMap(i,j) == 2 || Val_image(i,j)  == 2
                UnionCount_identities = UnionCount_identities + 1;
                
                if CellIdentitiesMap(i,j) == 2 && Val_image(i,j) == 2
                    InterCount_identities = InterCount_identities + 1;
                    
                end
            end
            
        end
    end
    
    Jaccard_Index_Pixels = InterCount/UnionCount
    Jaccard_Identities = InterCount_identities/UnionCount_identities;
    
    CellIdentityAccuracy = Jaccard_Identities/Jaccard_Index_Pixels
end

%% Network Construction

%March through the image around the boundaries of the in labelled cell
%record the labels for each cell

  
%Finding the connecting regions
for ii = 1:numel(ListOfLabel)
    DontConnect = [RemoveLabels,CellData(ii,1).LabelId];
    AdjCells = [];
    AdjCells = [AdjCells;unique(Seg_image(bwdist(Seg_image == CellData(ii,1).LabelId,'chessboard')<=AdjacencyKeneral))];
    AdjCells = [AdjCells;unique(Seg_image(bwdist(Seg_image == CellData(ii,1).LabelId,'cityblock')<=AdjacencyKeneral))];
    AdjCells = unique(AdjCells);
    AdjCells = setdiff(AdjCells,DontConnect );
    CellData(ii,1).AdjacentCells = AdjCells';
end


% Construct adjacency matrix
%AdjacencyMat = zeros(length(ListOfLabel),length(ListOfLabel));
AdjacencyMat = ConstructAdjacencyMatrix(CellData,w1,w2,1);

plotlinesX1 = [];
plotlinesX2 = [];
plotlinesY1 = [];
plotlinesY2 = [];
CentroidsX_lum = [];
CentroidsY_lum = [];
CentroidsX_basal = [];
CentroidsY_basal = [];


for ii = 1:length(ListOfLabel)
    for jj = 1:length(ListOfLabel)
        
        
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

%% final processing of outputs before visualing

%smooth mark boundaries using convolution

Mask = boundarymask(Seg_image);
image_binary = single(Mask);
kernalSize=3;  % Decide as per your requirements
kernel=ones(kernalSize)/kernalSize^2;
result=conv2(single(image_binary),kernel,'same');
result=result>0.5;
image_binary(~result)=0; 


CellBoundaries = imcomplement(image_binary);

%seg_image with large cells removed
for i = 2:numel(RemoveLabels)
    Seg_image(Seg_image == RemoveLabels(i)) = RemoveLabels(1);
end
%% Plot results 
if PlotAllOutputs == 1
    close all
    
    ImIndentity = label2rgb(CellIdentitiesMap','lines');
    Seg_labels = label2rgb(Seg_image-1,'jet',[0,0,0],'shuffle');
    Mask = boundarymask(Seg_image);
    
    if ValidateImageSeg == 1
        ImValid= label2rgb(Val_image','lines');
        figure
        subplot(1,3,1)
        imshow(FullImageStack_bas' )
        title('Input image: basal')
        subplot(1,3,2)
        imshow(FullImageStack_lum' )
        
        title('Input image: luminal')
        subplot(1,3,3)
        imshow(ImValid)
        %imshow(labeloverlay(ImIndentity,Mask,'Transparency',0))
        title('Input image: segmented identities')
        
    end
    
    
    figure
    subplot(2,2,1);
    imshow(Raw_image)
    title('Input image')
    
    subplot(2,2,2);
    B_im = imshow(Bound_image);
    title('U-net: boundary probabilities')
    
    subplot(2,2,3);
    %imshow(labeloverlay(Seg_labels,Mask,'Transparency',0))
    imshow(Seg_labels)
    title('GASP: segmentation')
    
    subplot(2,2,4)
    Boundaries = imshow(CellBoundaries);
    alpha(Boundaries,0.5)
    hold on
    plot([plotlinesX1;plotlinesX2 ],[plotlinesY1;plotlinesY2],'linewidth',1.5, 'color','k')
    scatter(CentroidsX_basal,CentroidsY_basal,50,[0 0.4470 0.7410],'filled')
    scatter(CentroidsX_lum,CentroidsY_lum,50,[0.8500 0.3250 0.0980],'filled')
    
    axis equal
    axis off
    title("Connectivity graph")
    



    figure;
        ImValid_dark = ImValid;
        ImIndentity_dark = ImIndentity;
    for i = 1:1024
        for j = 1:1024
            if ImValid_dark(i,j,1) == 255 && ImValid_dark(i,j,2) == 255 && ImValid_dark(i,j,3) == 255
                    ImValid_dark(i,j,1) = 0;
                    ImValid_dark(i,j,2) = 0;
                    ImValid_dark(i,j,3) = 0;
            end

            if ImIndentity_dark(i,j,1) == 255 && ImIndentity_dark(i,j,2) == 255 && ImIndentity_dark(i,j,3) == 255
                    ImIndentity_dark(i,j,1) = 0;
                    ImIndentity_dark(i,j,2) = 0;
                    ImIndentity_dark(i,j,3) = 0;
            end
        end
    end

       

    subplot(1,2,1)
    imshow(ImValid_dark)

    subplot(1,2,2)
    imshow(imoverlay(ImIndentity_dark,Mask,'black'))

end


%% Save outputs for analysis
if SaveOutputs == 1
    outputFileName = strcat('/Cleared_R2_comp_membrane00',num2str(SliceOfInterest-1),'_PS_pipeOutput');
    save(strcat(OutputDirect,outputFileName), 'CellData','Seg_image','CellBoundaries',...
        'Bound_image','plotlinesX1','plotlinesX2','plotlinesY1','plotlinesY2','CentroidsX_basal','CentroidsY_basal','CentroidsX_lum','CentroidsY_lum')
    
end




