% results for the contrast for the y value
y = [ystar1, ystar2, ystar4];
c = [1, -1, 0; 1 0 -1; 0, 1, -1];
cont = [x1,x2,x4]*c';
ts = tinv([0.025, 0.975], 99);

mean(cont)
std(cont)
[mean(cont,1);mean(cont,1)] + std(cont,1).*ts'/10