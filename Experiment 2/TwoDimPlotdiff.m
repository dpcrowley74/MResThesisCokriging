function diff = TwoDimPlotdiff(model)
% function for plotting the 2 dimnensional surrogate function from model
[x,y] = meshgrid(-1:0.02:1,-1:0.02:1);
zhat = zeros(101);
z = zeros(101);
for i = 1:101
    for j = 1:101
        zhat(i,j) = pred ([x(i,j),y(i,j)],model,0);
        z(i,j) = rastrigin([x(i,j),y(i,j)],10000);
    end
end
diff = abs(zhat - z);
figure(1)
surf(x,y,diff);
figure(2)
surf(x,y,z);
figure(3)
surf(x,y,zhat);
end