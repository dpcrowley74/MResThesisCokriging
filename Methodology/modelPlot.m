rng(12)
x = 0:0.01:1;
y = testfunction(x);
xc = 0:0.25:1;
xc = xc';
yc = testfunction(xc,1);
xe = xc([1,3,5]);
ye = testfunction(xe);
model = Cokriging(xe,ye,xc,yc);
model1 = Cokriging(xe,ye,xe,ye);
yhat = zeros(101,1);
ykrig = zeros(101,1);
for i = 1:101
    yhat(i) = pred(x(i),model);
    ykrig(i) = pred(x(i),model1,1);
end
figure(1)
plot(x, y, x, yhat, x, ykrig, xc, yc, 'x', xe,ye,'*')
title("Test Function with Models based of Sample")
xlabel('x')
ylabel('y')
legend('$f_e(x)$','$\hat{f}_e(x)$','Kriging','$s_c$','$s_e$', 'Location','northwest', 'Interpreter', 'latex')