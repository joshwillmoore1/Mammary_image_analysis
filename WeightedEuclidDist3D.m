function Distance = WeightedEuclidDist3D(X,Y)
    
    
    Weights = [1.5,1.5,0.6];
   
    for i = 1:length(Y(:,1))
     Distance(i) = sqrt(sum((Weights.*(X - Y(i,:))).^2));

    end
    
end