clear;clc;
global L a b 
L = 1; rho = 1.25/L;
% L = 1/2; rho = 2;
a = rho/2*L^2; b = L*sqrt(1 - rho^2/4*L^2);

num_iter = 80;     % Number of iterations
x_1 = 10; y_1 = 10;   % Initial point

maxiter = 1000;
% SA-MGDA
gd_norm_sq = zeros(num_iter,maxiter);
t_gd = zeros(num_iter,maxiter);
x = cell(1,num_iter); y = cell(1,num_iter);
x{1} = x_1; y{1} = y_1;
for i = 1:maxiter
    for k = 1:num_iter
        tic

        x{k+1} = x{k} - 1/L*gd_x(x{k},y{k});
        y{k+1} = (y{k} + 2/L*b*x{k+1} - 1/L*gd_y(x{k},y{k}))/(1+2/L*a);
        t_gd(k+1,i) = t_gd(k,i) + toc;

        % gd_norm_sq(k) = gd_norm(x{k+1},y{k+1});
        gd_norm_sq(k,i) = sqrt(x{k+1}^2 + y{k+1}^2);
    end
end
t_gd = mean(t_gd,2);
gd_norm_sq = mean(gd_norm_sq,2);

%% PD
%Q = -0.55;
Q = -0.6;
gd_norm_pd1 = zeros(num_iter,maxiter);
t_pd1 = zeros(num_iter,maxiter);
x = cell(1,num_iter); y = cell(1,num_iter);
x{1} = x_1; y{1} = y_1;
theta = 1.1;% 1.15
for i = 1:maxiter
    for k = 1:num_iter
        tic

        delta_x = gd_x(x{k},y{k});
        delta_y = gd_y(x{k},y{k});

        x{k+1} = x{k} - 1/L*delta_x;
        bian_x = x{k}-x{k+1};

        y{k+1} = y{k} + 1/L*delta_y + Q*1/L*bian_x;
        r_k = [x{k+1},y{k+1}];
        u_k = [x{k},y{k}];
        M_ur = [1/L*bian_x , Q*bian_x + 1/L*(y{k}-y{k+1})]-...
            [delta_x-gd_x(x{k+1},y{k+1}) , -delta_y+gd_y(x{k+1},y{k+1})];
        u_kk = u_k - theta*M_ur*((u_k-r_k)*M_ur')/norm(M_ur,2)^2;

        t_pd1(k+1,i) = t_pd1(k,i) + toc;
        x{k+1} = u_kk(1); y{k+1} = u_kk(2);
        % gd_norm_pd1(k) = gd_norm(x{k+1},y{k+1});
        gd_norm_pd1(k,i) = sqrt(x{k+1}^2 + y{k+1}^2);
        % gd_norm_pd1(k,i) = norm(u_kk-u_k);
    end
end
gd_norm_pd1 = mean(gd_norm_pd1,2);
t_pd1 = mean(t_pd1,2);
% 
%% PD
Q = -0.5;% -0.7
gd_norm_pd2 = zeros(num_iter,maxiter);
t_pd2 = zeros(num_iter,maxiter);

x = cell(1,num_iter); y = cell(1,num_iter);
x{1} = x_1; y{1} = y_1;
theta = 1.5;% 1
for i = 1:maxiter
    for k = 1:num_iter
        tic

        delta_x = gd_x(x{k},y{k});
        delta_y = gd_y(x{k},y{k});

        x{k+1} = x{k} - 1/L*delta_x;
        bian_x = x{k}-x{k+1};
        y{k+1} = y{k} + 1/L*delta_y + Q*1/L*bian_x;
        r_k = [x{k+1},y{k+1}];
        u_k = [x{k},y{k}];

        M_ur = [1/L*bian_x , Q*bian_x + 1/L*(y{k}-y{k+1})]-[delta_x-gd_x(x{k+1},y{k+1}) , -delta_y+gd_y(x{k+1},y{k+1})];
        u_kk = u_k - theta*M_ur*((u_k-r_k)*M_ur')/norm(M_ur,2)^2;

        t_pd2(k+1,i) = t_pd2(k,i) + toc;
        x{k+1} = u_kk(1); y{k+1} = u_kk(2);
        % gd_norm_pd2(k) = gd_norm(x{k+1},y{k+1});
        gd_norm_pd2(k,i) = sqrt(x{k+1}^2 + y{k+1}^2);
        %gd_norm_pd2(k,i) = norm(u_kk-r_k);
    end
end
gd_norm_pd2 = mean(gd_norm_pd2,2);
t_pd2 = mean(t_pd2,2);
% 
% 
%% Plot
mediumgray = [0.6 0.6 0.6];
mediumred = [0.6350 0.0780 0.1840];
mediumgreen = [0.4660 0.6740 0.1880];
mediumblue = [0 0.4470 0.7410];
mediumblack = [0.2 0.2 0.2];
mediumyellow = [0.8500 0.3250 0.0980];
mediumpurple = [0.4940 0.1840 0.5560];
mediumorange = [0.9290 0.6940 0.1250];

subplot(1,2,1)
plot(1:num_iter,log10(gd_norm_pd1),'Color',mediumgreen,'LineWidth',3);hold on
plot(1:num_iter,log10(gd_norm_sq),'Color',mediumblue,'LineWidth',3);hold on 

plot(1:num_iter,log10(gd_norm_pd2),'Color',mediumred,'LineWidth',3);hold on
% plot(1:maxiter,log10(meanerr_pdQ),'Color',mediumyellow,'LineWidth',2);
set(gca,'linewidth',3,'fontsize',18,'fontname','Times');
legend('Non-PPP,Q=-0.6,\theta=1.1','SA-MGDA','Non-PPP,Q=-0.5,\theta=1.5','Fontsize',18,'Location','NorthEast');
%legend('Generalized CP-PPA','Non-PPP,Q=-\alpha A','Fontsize',18,'Location','NorthEast');
set(gca,'linewidth',3,'fontsize',18,'fontname','Times');
xl = xlabel('iteration'); 
xlim([0 40]);
yl = ylabel('log($\|u_k\|^2$)','interpreter','latex');
hold off

subplot(1,2,2)
plot(t_pd1(2:end),log10(gd_norm_pd1),'Color',mediumgreen,'LineWidth',3);hold on
plot(t_gd(2:end),log10(gd_norm_sq),'Color',mediumblue,'LineWidth',3);hold on 

plot(t_pd2(2:end),log10(gd_norm_pd2),'Color',mediumred,'LineWidth',3);hold on
% plot(1:maxiter,log10(meanerr_pdQ),'Color',mediumyellow,'LineWidth',2);
set(gca,'linewidth',3,'fontsize',18,'fontname','Times');
legend('Non-PPP,Q=-0.6,\theta=1.1','SA-MGDA','Non-PPP,Q=-0.5,\theta=1.5','Fontsize',18,'Location','NorthEast');
% legend('Generalized CP-PPA','Non-PPP,Q=-\alpha A', 'Fontsize',18,'Location','NorthEast');
set(gca,'linewidth',3,'fontsize',18,'fontname','Times');
xl = xlabel('time(s)'); 
xlim([0 0.0006]);
yl = ylabel('log($\|u_k\|^2$)','interpreter','latex');
hold off


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