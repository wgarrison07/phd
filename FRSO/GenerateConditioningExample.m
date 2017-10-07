close all;
dim = 20;
points = RandomLatinHypercube([0,0;1,1],dim);
figure();
plot(points(:,1), points(:,2), 'b.', 'MarkerSize', 25);
Phi = zeros(dim,dim);
for i = 1:dim
    for j = 1:dim
        r = norm(points(i,:)-points(j,:));
        Phi(i, j) = r;
    end
end
sigsq = var(Phi(:));
Phi = exp(-Phi.^2/(2*sigsq));
fprintf(1, 'Matrix Condition: %d\n', rcond(Phi));

hold on;
points = [points; points(1,:) + [0.00,0.01]; points(1,:) + [-0.008, -0.008]; points(1,:) + [0.008, -0.008]];
plot(points(dim+1:dim+3,1), points(dim+1:dim+3,2), 'r.', 'MarkerSize', 25);
dim = dim + 3;
Phi = zeros(dim, dim);
for i = 1:dim
    for j = 1:dim
        r = norm(points(i,:)-points(j,:));
        Phi(i, j) = exp(-r^2);
    end
end
fprintf(1, 'Matrix Condition: %d\n', rcond(Phi));
