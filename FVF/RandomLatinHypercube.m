function points = RandomLatinHypercube(varRanges, numPoints)

%Initialize the points
dim = size(varRanges, 2);
points = zeros(numPoints, dim);

for i = 1:dim
    perm = GetPermutationSequence(numPoints);
    for j = 1:numPoints
        points(j,i) = (perm(j) - 1) / (numPoints - 1);
    end
    points(:,i) = points(:,i)*(varRanges(2,i) - varRanges(1,i)) + varRanges(1,i);
end

end

function sequence = GetPermutationSequence(numPoints)
    sequence = 1:numPoints;
    for i = 1:numPoints-1
        %Generate random number between i and n
        j = ceil(rand()*(numPoints - i)) + i;
        temp = sequence(i);
        sequence(i) = sequence(j);
        sequence(j) = temp;
    end
end
