x = 0:0.01:1;
y = testfunction(x);
xc = [0,0.25,0.5,0.75,1];
yc = testfunction(xc,1);
xe = [0,0.5,1];
ye = testfunction(xe);
figure(1)
plot(x, y, xc, yc, 'x', xe,ye,'*')
title("Test Function with Samples")
xlabel('x')
ylabel('y')
legend('$f_e(x)$','$s_c$','$s_e$', 'Location','northwest', 'Interpreter', 'latex')