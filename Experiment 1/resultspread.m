function [y,x]= resultspread(ystar, ystarvar, dim)
% function for generating the plots for result graphs
x = initialdesign(100,dim,-1,1);
y = rastrigin(x, 10000);

figure(1)
histogram(y)
hold on
plot([min(ystar),max(ystar)], zeros(2,1),'*','linewidth', 10)
plot([min(ystarvar),max(ystarvar)], zeros(2,1),'*', 'linewidth', 10)
legend('Latin HyperCube Sample', 'Single Fidelity Range', 'Variable Fidelity Range', 'location', 'NorthWest')
hold off

