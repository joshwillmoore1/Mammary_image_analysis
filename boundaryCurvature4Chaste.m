close all
clc
clear
NumberOfPoints = 1000;
t = linspace(0,2*pi,NumberOfPoints);
t_red = linspace(0,2*pi,75);
tspaced = [t,t,t];
L = 1;
A = 0; %0 0.25 1.2 1.7
r = L*5;
n = 5;
Rl = 1.2*L;
Centre = [0,0];

for i = 1:length(t_red)
   xR(i) = 0.6.*R(A,r,n,t_red(i)).*cos(t_red(i))+ Centre(1);
    yR(i) = 0.6.*R(A,r,n,t_red(i)).*sin(t_red(i))+ Centre(2); 
end

for i = 1:(length(t))
    
    xI(i) = Rl.*cos(t(i));
    yI(i) = Rl.*sin(t(i));
    
end


boundarygap = 11;
for i = 1:(length(t)+1)
    
    t0 = tspaced(NumberOfPoints + i);
    t1 = tspaced(NumberOfPoints + i-boundarygap);
    t2 = tspaced(NumberOfPoints + i+boundarygap);
    
    x(i) =  R(A,r,n,t0).*cos(t0) + Centre(1);
    y(i) =  R(A,r,n,t0).*sin(t0)+ Centre(2);

    %point 1
    x1(i) =  R(A,r,n,t1).*cos(t1)+ Centre(1);
    y1(i) =  R(A,r,n,t1).*sin(t1)+ Centre(2);
    
    xp1(i) = Rp(A,r,n,t1).*cos(t1) - R(A,r,n,t1).*sin(t1) + Centre(1);
    yp1(i) = Rp(A,r,n,t1).*sin(t1) + R(A,r,n,t1).*cos(t1)+ Centre(2);
    
    %normalise 1
    TangentNorm1(i) = norm([xp1(i),yp1(i)],2);
    xpN1(i) = xp1(i)/TangentNorm1(i);
    ypN1(i) = yp1(i)/TangentNorm1(i);
    
    %point 2
    x2(i) =  R(A,r,n,t2).*cos(t2)+ Centre(1);
    y2(i) =  R(A,r,n,t2).*sin(t2)+ Centre(2);
    xp2(i) = Rp(A,r,n,t2).*cos(t2) - R(A,r,n,t2).*sin(t2)+ Centre(1);
    yp2(i) = Rp(A,r,n,t2).*sin(t2) + R(A,r,n,t2).*cos(t2)+ Centre(2);
    
    %normalise 2
    TangentNorm2(i) = norm([xp2(i),yp2(i)],2);
    xpN2(i) = xp2(i)/TangentNorm2(i);
    ypN2(i) = yp2(i)/TangentNorm2(i);
    
    
    grad1 = ypN1(i)/xpN1(i);
    grad2 = ypN2(i)/xpN2(i);
    
    %if the tangent lines are parallel then set the curvature to 0
    if abs(grad1 - grad2) < 1e-5
        
        k(i) = 0;
        
    else
        
        %else find the intersection point of the normals at the point
        NormalX1(i) = -ypN1(i) + x1(i);
        NormalY1(i) = xpN1(i) + y1(i);
        
        %grad and intercept
        
        M1(i) = (NormalY1(i) - y1(i))./(NormalX1(i) - x1(i));
        C1(i) = NormalY1(i) - M1(i).*NormalX1(i);
        
        
        NormalX2(i) = -ypN2(i) + x2(i);
        NormalY2(i) = xpN2(i) + y2(i);
        
        M2(i) = (NormalY2(i) - y2(i))./(NormalX2(i) - x2(i));
        C2(i) = NormalY2(i) - M2(i).*NormalX2(i);
        
        CentreX(i) = (C2(i) - C1(i))/(M1(i) - M2(i));
        CentreY(i) = M1(i)*CentreX(i) + C1(i);
        
        
        Rad(i) = norm([CentreX(i),CentreY(i)]- [x1(i),y1(i)  ] )*7.92;
       
        
        
        %take a ray through the centre of the circle and one of the
        %boundary points - if at least one is inside the shape then the
        %curvature is positive
        
        cent2boundX = linspace(CentreX(i),x2(i),10);
        
        midpointX = (x1(i) + x2(i))/2;
        midpointY = (y1(i) + y2(i))/2;
        tempT = atan2(midpointY, midpointX);
        
        if R(A,r,n,tempT) > norm([midpointX  ,midpointY ])
            
            k(i) = 1/Rad(i);
        else
            k(i) = -1/Rad(i);
        end
        
        

        
        

        
        
        
        
    end
    
    
    
   [MinDist,MinIndx] = min( vecnorm( [xR',yR'] - [x(i),y(i)],2,2));
   
   closestInd(i) = MinIndx;
    
    
    
    
end
z = zeros(1,length(t)+1);

close all
max(k)
min(k)
figure('renderer','painters')
surf([x',x'],[y',y'],[z',z'], [k', k'],'FaceColor', 'none', ...    % Don't bother filling faces with color
    'EdgeColor', 'interp', ...  % Use interpolated color for edges
    'LineWidth', 5)
hold on
colormap('jet')
c = colorbar('location','north','Visible','off');
c.Label.Interpreter = 'latex';
c.TickLabelInterpreter = 'latex';

caxis([-0.3,0.2])
UniClosId = unique(closestInd);
meanPlot = [];
for i = 1:length(UniClosId)
    
    
    TempX = xR(UniClosId(i));
    TempY = yR(UniClosId(i));
    
    
    MeanK = median(k(closestInd == UniClosId(i)));
    meanPlot = [meanPlot; [TempX,TempY,MeanK]];
    
    
end



 %scatter(meanPlot(:,1),meanPlot(:,2),100,meanPlot(:,3),'filled')
 % scatter(xI,yI,100,'filled')
%colorbar
view(2)
grid off
axis equal
axis off



function out = R(A,r,n,t)

out  =  r + A.*sin(n.*t);

end

function out = Rp(A,r,n,t)

out  =  A.*n.*cos(n.*t);

end

function out = Rpp(A,r,n,t)

out  =  -A.*(n.^2).*sin(n.*t);

end