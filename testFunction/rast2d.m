[x,y] = meshgrid(-1:0.01:1,-1:0.01:1);
z = zeros(201);
z = zeros(201);
for i = 1:201
    for j = 1:201
        zl(i,j) = rastrigin([x(i,j),y(i,j)],1000, 1);
        z(i,j) = rastrigin([x(i,j), y(i,j)], 10000 ,1);
    end
end

figure(1)
surf(x,y,z)
xlabel("x");
ylabel("y");
zlabel("z");
title("High Fidelity Function")

figure(2)
surf(x,y,zl)
xlabel("x");
ylabel("y");
zlabel("z");
title("Low Fidelity Function")