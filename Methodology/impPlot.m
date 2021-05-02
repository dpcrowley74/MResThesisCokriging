rng(1)
x = 0:0.01:1;
y = testfunction(x);
xc = 0:0.25:1;
xc = xc';
yc = testfunction(xc,1);
xe = xc([1,3,5]);
ye = testfunction(xe);
model = Cokriging(xe,ye,xc,yc);
yhat = zeros(101,1);
hiImp = zeros(101,1);
lowImp = zeros(101,1);
for i = 1:101
    hiImp(i) = VarImp(x(i),model,0,0,0);
    lowImp(i) = VarImp(x(i),model,1,0,0);
end
ycz = zeros(size(xc));
yez = zeros(size(xe));
figure(1)
plot( x, hiImp, x, lowImp, xc, ycz, 'x', xe,yez,'*')
title("Test Function with Models based of Sample")
xlabel('x')
ylabel('y')
legend('Expected Improvement for High Fidleity', 'Expected Improvement for Low Fidelity','$s_c$','$s_e$', 'Location','northeast', 'Interpreter', 'latex')