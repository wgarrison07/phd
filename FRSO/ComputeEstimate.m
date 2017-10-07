function est = ComputeEstimate(pos, coef, powerMat)

%Initialize Outputs
dim = length(pos);
est = zeros(dim+1, 1);


%Compute Data
for i = 1:length(coef)
   prod = repmat(coef(i),dim+1,1); 
   for j = 1:dim
       if powerMat(i,j) == 0
           prod(1+j) = 0;
       else
           prod(1) = prod(1) * (pos(j)^powerMat(i,j));
           prod(1+j) = prod(1+j) * (powerMat(i,j)*pos(j)^(powerMat(i,j)-1));
       end
   end
   est = est + prod;
end

end