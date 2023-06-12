close all
clc
clear


R = 0.5;

x  = linspace(-2,2,50);
y  = linspace(-2,2,50);
z  = linspace(-2,2,50);

Accept = zeros(50,50,50);

for i = 1:50
    for j = 1:50
        for k = 1:50
            
            if WeightedEuclidDist3D([0,0,0],[x(i),y(j),z(k)]) <= R
                
                Accept(i,j,k) = 1;
                
            end
            
        end
        
    end
end


isosurface(Accept)
axis equal