rng(1)
x = 0:0.01:1;
y = testfunction(x);
xc = 0:0.25:1;
xc = xc';
yc = testfunction(xc,1);
xe = xc([1,2,3,5]);
ye = testfunction(xe);
model = Cokriging(xe,ye,xc,yc);
yhat = zeros(101,1);
hiImp = zeros(101,1);
lowImp = zeros(101,1);
for i = 1:101
    yhat(i) = pred(x(i),model);
    hiImp(i) = VarImp(x(i),model,0,0,0);
    lowImp(i) = VarImp(x(i),model,1,0,0);
end
ycz = zeros(size(xc));
yez = zeros(size(xe));
figure(1)
plot(x, hiImp, x, lowImp, xc, ycz, 'x', xe,yez,'*')
title("Improvement function and Samples")
xlabel('x')
ylabel('y')
legend('Expected Improvement for High Fidleity', 'Expected Improvement for Low Fidelity','$s_c$','$s_e$', 'Location','northwest', 'Interpreter', 'latex')

figure(2)
plot(x, y, x, yhat, xe,ye, '*',xc,yc,'x')
title("Objective Function and New Model")
xlabel("x")
ylabel("y")
legend('$f_e(x)$','$\hat{f}_e(x)$','$s_c$','$s_e$', 'Location','northwest', 'Interpreter', 'latex')