clear;clc;
global L a b 
L = 1; rho = 1.25/L;
% L = 1/2; rho = 2;
a = rho/2*L^2; b = L*sqrt(1 - rho^2/4*L^2);

num_iter = 30;     % Number of iterations
x_1 = 10; y_1 = 10;   % Initial point


%% SA-MGDA
% gd_norm_sq = zeros(1,num_iter);
% x = cell(1,num_iter); y = cell(1,num_iter);
% x{1} = x_1; y{1} = y_1;
% for k = 1:num_iter
%     x{k+1} = x{k} - 1/L*gd_x(x{k},y{k});
%     y{k+1} = (y{k} + 2/L*b*x{k+1} - 1/L*gd_y(x{k},y{k}))/(1+2/L*a);
%     gd_norm_sq(k) = gd_norm(x{k+1},y{k+1});
% end
% gd_norm_sq_best = best(gd_norm_sq);

%% PD
Q = -0.5;
gd_norm_pd1 = zeros(1,num_iter);
x = cell(1,num_iter); y = cell(1,num_iter);
x{1} = x_1; y{1} = y_1;
theta = 0.5;
for k = 1:num_iter
    x{k+1} = x{k} - 1/L*gd_x(x{k},y{k});
    y{k+1} = y{k} + 1/L*gd_y(x{k},y{k}) + Q*1/L*(x{k}-x{k+1});
    r_k = [x{k+1},y{k+1}];
    u_k = [x{k},y{k}];
    M_ur = [1/L*(x{k}-x{k+1}) , Q*(x{k}-x{k+1}) + 1/L*(y{k}-y{k+1})]-...
        [gd_x(x{k},y{k})-gd_x(x{k+1},y{k+1}) , -gd_y(x{k},y{k})+gd_y(x{k+1},y{k+1})];
    P_k = u_k - M_ur*((u_k-r_k)*M_ur')/norm(M_ur,2)^2;
    u_kk = (1-theta)*u_k + theta*P_k;
    x{k+1} = u_kk(1); y{k+1} = u_kk(2);

%     gd_norm_pd1(k) = gd_norm(x{k+1},y{k+1});
gd_norm_pd1(k) = norm(u_kk-r_k);
end
gd_norm_pd_best1 = best(gd_norm_pd1);

%% PD
Q = -0.5;
gd_norm_pd2 = zeros(1,num_iter);
x = cell(1,num_iter); y = cell(1,num_iter);
x{1} = x_1; y{1} = y_1;
theta = 1;
for k = 1:num_iter
    x{k+1} = x{k} - 1/L*gd_x(x{k},y{k});
    y{k+1} = y{k} + 1/L*gd_y(x{k},y{k}) + Q*1/L*(x{k}-x{k+1});
    r_k = [x{k+1},y{k+1}];
    u_k = [x{k},y{k}];
    M_ur = [1/L*(x{k}-x{k+1}) , Q*(x{k}-x{k+1}) + 1/L*(y{k}-y{k+1})]-...
        [gd_x(x{k},y{k})-gd_x(x{k+1},y{k+1}) , -gd_y(x{k},y{k})+gd_y(x{k+1},y{k+1})];
    P_k = u_k - M_ur*((u_k-r_k)*M_ur')/norm(M_ur,2)^2;
    u_kk = (1-theta)*u_k + theta*P_k;
    x{k+1} = u_kk(1); y{k+1} = u_kk(2);

    %gd_norm_pd2(k) = gd_norm(x{k+1},y{k+1});
    gd_norm_pd2(k) = norm(u_kk-r_k);
end
gd_norm_pd_best2 = best(gd_norm_pd2);

%% PD
Q = -0.5;
gd_norm_pd3 = zeros(1,num_iter);
x = cell(1,num_iter); y = cell(1,num_iter);
x{1} = x_1; y{1} = y_1;
theta = 1.5;
for k = 1:num_iter
    x{k+1} = x{k} - 1/L*gd_x(x{k},y{k});
    y{k+1} = y{k} + 1/L*gd_y(x{k},y{k}) + Q*1/L*(x{k}-x{k+1});
    r_k = [x{k+1},y{k+1}];
    u_k = [x{k},y{k}];
    M_ur = [1/L*(x{k}-x{k+1}) , Q*(x{k}-x{k+1}) + 1/L*(y{k}-y{k+1})]-...
        [gd_x(x{k},y{k})-gd_x(x{k+1},y{k+1}) , -gd_y(x{k},y{k})+gd_y(x{k+1},y{k+1})];
    P_k = u_k - M_ur*((u_k-r_k)*M_ur')/norm(M_ur,2)^2;
    u_kk = (1-theta)*u_k + theta*P_k;
    x{k+1} = u_kk(1); y{k+1} = u_kk(2);

    %gd_norm_pd3(k) = gd_norm(x{k+1},y{k+1});
    gd_norm_pd3(k) = norm(u_kk-r_k);
end
gd_norm_pd_best3 = best(gd_norm_pd3);

%% Plot
mediumgray = [0.6 0.6 0.6];
mediumred = [0.6350 0.0780 0.1840];
mediumgreen = [0.4660 0.6740 0.1880];
mediumblue = [0 0.4470 0.7410];
mediumblack = [0.2 0.2 0.2];
mediumyellow = [0.8500 0.3250 0.0980];
mediumpurple = [0.4940 0.1840 0.5560];
mediumorange = [0.9290 0.6940 0.1250];

plot(1:num_iter,log10(gd_norm_pd_best1),'Color',mediumblue,'LineWidth',3);hold on 
plot(1:num_iter,log10(gd_norm_pd_best2),'Color',mediumgreen,'LineWidth',3);hold on
plot(1:num_iter,log10(gd_norm_pd_best3),'Color',mediumred,'LineWidth',3);hold on
% plot(1:maxiter,log10(meanerr_pdQ),'Color',mediumyellow,'LineWidth',2);
set(gca,'linewidth',3,'fontsize',18,'fontname','Times');
legend('Q=-0.5, \theta=0.5','Q=-0.5, \theta=1','Q=-0.5, \theta=1.5','Fontsize',18,'Location','NorthEast');
% legend('Generalized CP-PPA','Non-PPP,Q=-\alpha A', 'Non-PPP,Q=(\alpha)A','Fontsize',18,'Location','NorthEast');
set(gca,'linewidth',3,'fontsize',18,'fontname','Times');
xl = xlabel('iteration'); 
yl = ylabel('log($\|u_k-r_k\|^2$)','interpreter','latex');
hold off

%% Plot
% color = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250],...
%     [0.4940, 0.1840, 0.5560],[0.4660, 0.6740, 0.1880]};
% 
% linewidth=2;
% figure
% hold on
% %plot(iter,gd_norm_sq_best,'-x','MarkerIndices',20,'Linewidth',linewidth)
% p1 = plot(1:num_iter,gd_norm_pd_best1,'-x','MarkerIndices',20,'Linewidth',linewidth);
% p2 = plot(1:num_iter,gd_norm_pd_best2,'-x','MarkerIndices',20,'Linewidth',linewidth);
% p3 = plot(1:num_iter,gd_norm_pd_best3,'-x','MarkerIndices',20,'Linewidth',linewidth);
% legend('Q=-0.5, $\theta$=1.5','Q=0, $\theta$=1.5','Q=0.5, $\theta$=1.5','Fontsize',12,'Location','southeast','interpreter','latex')
% set(gca, 'XScale', 'log')
% set(gca, 'YScale', 'log')
% ylim([1e-20,10])
% xlim([1,1e3])
% ylabel('$\|Bu_k\|^2$','interpreter','latex')
% xlabel('Iterations','interpreter','latex')
% set(gca,'FontSize',15)
% hold off

%% Functions
function val = clip(x)
    global lb ub
    val = max(lb,min(ub,x));
end

function vec = best(x)
    vec = zeros(1,length(x));
    for i = 1:length(x)
        vec(i) = min(x(1:i));
    end
end

function grad = gd_x(x,y)
    global a b
    grad = a*x + b*y;
end

function grad = gd_y(x,y)
    global a b
    grad = b*x - a*y;
end

function value = gd_norm(x,y)
    value = norm(gd_x(x,y))^2 + norm(gd_y(x,y))^2;
end