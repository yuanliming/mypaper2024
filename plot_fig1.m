figure
x=1:1:16;
scatter(x,h2norm_Ksv,'s','b','filled')
box on
hold on
scatter(x,h2norm_Kqs,'d','r','filled')

%xlim([1 256])
set(gca,'FontSize',16,'Fontname', 'Times New Roman')
leg=legend('$K_{\rm SV}$','$K_{\rm QS}$');
%legend('Slacking-variable method','Quadratic stability based method')
xla=xlabel('Extreme systems');
yla=ylabel('$H_2$ cost');
set(leg,'interpreter','latex')
set(xla,'interpreter','latex')
set(yla,'interpreter','latex')
% plot([0,16],[3.2014,3.2014],'b','LineWidth',1.5)
% plot([0,16],[3.3479,3.3479],'r','LineWidth',1.5)
grid on