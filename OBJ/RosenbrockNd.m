function val = RosenbrockNd(input)

val = 0;
for i = 1:length(input)-1
   val = val + 100*(input(i+1) - input(i)^2)^2 + (input(i) - 1)^2;
end

end