n = 500;
m = n/4; s= n/20;
data1 = load('data.mat');
data1 = data1.data;
subplot(1,1,1)
mediumgray = [0.6 0.6 0.6];
mediumred = [0.6350 0.0780 0.1840];
mediumgreen = [0.4660 0.6740 0.1880];
mediumblue = [0 0.4470 0.7410];
mediumblack = [0.2 0.2 0.2];
mediumyellow = [0.8500 0.3250 0.0980];

theta1 = data1(1,:);
theta2 = data1(2,:);
theta3 = data1(3,:);


hold off


plot(0.1:0.1:1,log10(theta1),'Color',mediumgreen,'LineWidth',3),hold on
plot(0.1:0.1:1,log10(theta2),'Color',mediumblue,'LineWidth',3),hold on
plot(0.1:0.1:1,log10(theta3),'Color',mediumred,'LineWidth',3),hold on

set(gca,'linewidth',3,'fontsize',18,'fontname','Times');
legend('\theta = 0.5', '\theta = 1', '\theta = 1.5',...
    'Location','NorthEast','Fontsize', 18);
%title('Dual objective')
set(gca,'linewidth',3,'fontsize',18,'fontname','Times');
xl = xlabel('\alpha'); 
%set(xl,'Fontsize',30)
yl = ylabel('log(Error) at 200 iterations');
%set(yl,'Fontsize',30)
hold off

tit = title([sprintf('m = %d, n = %d, s = %d',m,n,s)]);
grid on;
set(tit,'FontSize',18);