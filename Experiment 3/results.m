% code for finding the results from the comparison for the x-coordinate
x1 = sqrt(sum(xstar1.^2,2));
x2 = sqrt(sum(xstar2.^2,2));
x4 = sqrt(sum(xstar4.^2,2));

c = [1, -1, 0; 1 0 -1; 0, 1, -1];
cont = [x1,x2,x4]*c';
ts = tinv([0.025, 0.975], 99);
mean(cont)
std(cont)
[mean(cont,1);mean(cont,1)] + std(cont,1).*ts'/10
