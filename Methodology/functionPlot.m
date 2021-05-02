x = 0:0.01:1;
ye = testfunction(x);
yc = testfunction(x,1);

figure(1)
plot(x, ye, x, yc, '--')
title("Test Function")
xlabel('x')
ylabel('y')
legend('$f_e(x)$','$f_c(x)$', 'Location','northwest', 'Interpreter', 'latex')
