figure
a=0.05:0.05:5;
plot(a,revo_alpha,'b','LineWidth',1.5)
hold on 
plot([0,5],[3.3479,3.3479],'r--','LineWidth',1.5)
box on
grid on
set(gca,'FontSize',16,'Fontname', 'Times New Roman')
xla=xlabel('$\alpha$');
yla=ylabel('$H_2$-norm upper bound');%Upper bound to the $H_2$ norm 
set(xla,'interpreter','latex')
set(yla,'interpreter','latex')
legend('Slacking-variable method','Quadratic-stability-based method')