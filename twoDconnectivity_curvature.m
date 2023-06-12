close all

close all
clc
clear

%get this from the metadata in the tif
x_pixel_scale= 0.2203746;
y_pixel_scale = 0.2203746;
z_pixel_scale = 2.9982000;
w1 = 1;
w2 = 1;

%Curvature paramaters
boundaryPoint = 10;         %12, 10 number of boundary points curvature is found over 
curvatureThresh = 0.1;     %0.06 the maximum allowed value of the curvature measure
bp_tangent = 10;            % number of boundary points the tangent angle is found over 
interpdmin = 0.15;           % the minimum number of pixels seperating boundary points after interpolation
loopclose = 1;    

sliceOfInterest = 13;
ResizeRatio = 1;


Files=dir(fullfile('PlantSeg_pipeline_outputs/R2_cleared/*.mat'));

for i = 1:length(Files)
    
    CellDataZ(i) = load(strcat('PlantSeg_pipeline_outputs/R2_cleared/Cleared_R2_comp_membrane00',num2str(i),'_PS_pipeOutput'));
    
    Mask(:,:,i) = imresize(imcomplement(CellDataZ(i).CellBoundaries),ResizeRatio);
    FilledMask(:,:,i) = imbinarize( imgaussfilt(imfill(Mask(:,:,i),'holes'),3));
    
    
end
BoundFilledMask = cast(FilledMask,'double');
BoundFilledMask=imfill(BoundFilledMask,'holes');



%Curvature Calculations
CurvatureCal = struct;
for i = 1:length(Files)
[shape_details, Icurv] = curvature2D(FilledMask(:,:,i),...
    boundaryPoint, curvatureThresh, bp_tangent, interpdmin, loopclose);

CurvatureCal(i,1).ShapeDetails = shape_details;
CurvatureCal(i,1).Image = Icurv;
end


plotlinesX1 = [];
plotlinesX2 = [];
plotlinesY1 = [];
plotlinesY2 = [];
CentroidsX_lum = [];
CentroidsY_lum = [];
CentroidsX_basal = [];
CentroidsY_basal = [];


sizeOfData = size(CellDataZ(sliceOfInterest).CellData);


%
for ii = 1:sizeOfData(1)
    for jj = 1:sizeOfData(1)
        
        
        if ismember(CellDataZ(sliceOfInterest).CellData(jj,1).LabelId,CellDataZ(sliceOfInterest).CellData(ii,1).AdjacentCells)
            
            plotlinesX1 = [plotlinesX1 , CellDataZ(sliceOfInterest).CellData(ii,1).Centroid(1)];
            plotlinesX2 = [plotlinesX2 , CellDataZ(sliceOfInterest).CellData(jj,1).Centroid(1)];
            plotlinesY1 = [plotlinesY1 , CellDataZ(sliceOfInterest).CellData(ii,1).Centroid(2)];
            plotlinesY2 = [plotlinesY2 , CellDataZ(sliceOfInterest).CellData(jj,1).Centroid(2)];
            
            
        end
        
    end
    
    if CellDataZ(sliceOfInterest).CellData(ii,1).cellType == 'Basal'
        
        CentroidsX_basal = [CentroidsX_basal;CellDataZ(sliceOfInterest).CellData(ii,1).Centroid(1)];
        CentroidsY_basal = [CentroidsY_basal;CellDataZ(sliceOfInterest).CellData(ii,1).Centroid(2)];
    elseif CellDataZ(sliceOfInterest).CellData(ii,1).cellType == 'Luminal'
        
        CentroidsX_lum = [CentroidsX_lum;CellDataZ(sliceOfInterest).CellData(ii,1).Centroid(1)];
        CentroidsY_lum = [CentroidsY_lum;CellDataZ(sliceOfInterest).CellData(ii,1).Centroid(2)];
    end
end




% Plot curvature
figure;
hold on
for i = 1:1:length(Files)
X = x_pixel_scale.*CurvatureCal(i,1).ShapeDetails.XY(:,1);
Y = y_pixel_scale.*CurvatureCal(i,1).ShapeDetails.XY(:,2);
Z = i.*z_pixel_scale.*ones(size(X));
C = CurvatureCal(i,1).ShapeDetails.curvature'*1;

XX = [X(:) X(:)];
YY = [Y(:) Y(:)];
ZZ = [Z(:) Z(:)];
CC  = [C C];
% Plot
%imshow(Icurv(i))
surf(XX, YY, ZZ, CC, ...  % Reshape and replicate data
 'FaceColor', 'none', ...    % Don't bother filling faces with color
 'EdgeColor', 'interp', ...  % Use interpolated color for edges
 'LineWidth', 5);            % Make a thicker line


% 
% plot([plotlinesX1;plotlinesX2 ],[plotlinesY1;plotlinesY2],'linewidth',1.5, 'color','k')
%     scatter(CentroidsX_basal,CentroidsY_basal,50,[0 0.4470 0.7410],'filled')
%     scatter(CentroidsX_lum,CentroidsY_lum,50,[0.8500 0.3250 0.0980],'filled')
%     


end

cmap = jet;
colormap(cmap);
cb = colorbar;  % Add a colorbar
cb.Label.String = 'Curvature';
view(3)
axis equal