clear;clc;
n = 500; v = 10; maxiter = 1500;
m = n/4; sparse = n/20;

alpha = 0.5;

S = @(z,lambda) max(abs(z)-lambda,0).*sign(z);

num_repeats = 100;

theta = 1;

err_generalize = zeros(maxiter,num_repeats);
t_generalize = zeros(maxiter,num_repeats);
err_pd = zeros(maxiter,num_repeats);
t_pd = zeros(maxiter,num_repeats);
err_pdQ = zeros(maxiter,num_repeats);
t_pdQ = zeros(maxiter,num_repeats);

for repeats = 1:num_repeats
rng(69462991)
    A = randn(m,n);
    xhat = sparserandn(n,sparse);  % true solution
    b = A*xhat;
    L = norm(A'*A,2);

    r = sqrt(L*(1-alpha+alpha^2))/v; 
    s = sqrt(L*(1-alpha+alpha^2))*v;

    % gerneralized CP-PPA
    x = zeros(n,1); y = zeros(m,1);
    for k = 2:maxiter
        tic
        % gerneralized CP-PPA
        xn = S(x+A'*y/r,1/r); % x_{k+1}
        xbar = xn + alpha*(xn-x);
        ybar = y - (A*xbar-b)/s;
        yn = ybar - (1-alpha)*A*(xn-x)/s; % y_{k+1}
        % err_generalize(k,repeats) = norm(xn-x)/norm(xn);
        err_generalize(k,repeats) = norm(xn-xhat)/norm(xhat);
        x = xn; y = yn;

        t_generalize(k,repeats) = t_generalize(k-1,repeats) + toc;
        %       % CP-PPA
        %       xn = S(x-A'*y/r,1/r); % x_{k+1}
        %       y = y-(A*(2*xn-x)-b)/s;
        %       err_geralize(k,num_repeats) = norm(xn-x)/norm(xn);
        %       x = xn;
    end
    

    % PD
    x = zeros(n,1); y = zeros(m,1);
    Q = (alpha)*A;
    for k = 2:maxiter
        tic
        % Framework
        xn = S(x+A'*y/r,1/r); % x_{k+1}
        yn = y - ((A-Q)*x+Q*xn-b)/s;
        an = [r*(x-xn)+A'*(y-yn);(Q-A)*(x-xn)+s*(y-yn)];
        un = [x;y];
        rn = [xn;yn];
        un1 = un - theta*(un-rn)'*an*an/norm(an,2)^2;

        % err_pd(k,repeats) = norm(un1(1:n)-x)/norm(un1(1:n));
        err_pd(k,repeats) = norm(un1(1:n)-xhat)/norm(xhat);
        x = un1(1:n); y = un1((n+1):(n+m));

        t_pd(k,repeats) = t_pd(k-1,repeats) + toc;
    end

    x = zeros(n,1); y = zeros(m,1);
    Q = -(alpha+1)*A;
    for k = 2:maxiter
        tic
        % Framework
        xn = S(x+A'*y/r,1/r); % x_{k+1}
        yn = y - ((A-Q)*x+Q*xn-b)/s;
        an = [r*(x-xn)+A'*(y-yn);(Q-A)*(x-xn)+s*(y-yn)];
        un = [x;y];
        rn = [xn;yn];
        un1 = un - theta*(un-rn)'*an*an/norm(an,2)^2;

        % err_pd(k,repeats) = norm(un1(1:n)-x)/norm(un1(1:n));
        err_pdQ(k,repeats) = norm(un1(1:n)-xhat)/norm(xhat);
        x = un1(1:n); y = un1((n+1):(n+m));

        t_pdQ(k,repeats) = t_pdQ(k-1,repeats) + toc;
    end


end


meanerr_generalize = median(err_generalize,2);
% meanerr_generalize = best(meanerr_generalize);
meanerr_pd = median(err_pd,2);
% meanerr_pd = best(meanerr_pd);
meanerr_pdQ = median(err_pdQ,2);

mean_t_generalize = median(t_generalize,2);
mean_t_pd = median(t_pd,2);
mean_t_pdQ = median(t_pdQ,2);

mediumgray = [0.6 0.6 0.6];
mediumred = [0.6350 0.0780 0.1840];
mediumgreen = [0.4660 0.6740 0.1880];
mediumblue = [0 0.4470 0.7410];
mediumblack = [0.2 0.2 0.2];
mediumyellow = [0.8500 0.3250 0.0980];
mediumpurple = [0.4940 0.1840 0.5560];
mediumorange = [0.9290 0.6940 0.1250];

subplot(1,2,1)

plot(1:maxiter,log10(meanerr_generalize),'Color',mediumblue,'LineWidth',3);hold on 
plot(1:maxiter,log10(meanerr_pd),'Color',mediumgreen,'LineWidth',3);hold on
plot(1:maxiter,log10(meanerr_pdQ),'Color',mediumyellow,'LineWidth',3);
set(gca,'linewidth',3,'fontsize',18,'fontname','Times');
% legend('Generalized CP-PPA','Non-PPP','Fontsize',18,'Location','NorthEast');
legend('Generalized CP-PPA','Non-PPP,Q=\alpha A', 'Non-PPP,Q=-(\alpha+1)A','Fontsize',18,'Location','NorthEast');
set(gca,'linewidth',3,'fontsize',18,'fontname','Times');
xl = xlabel('iteration'); 
yl = ylabel('log(Error)');
hold off

subplot(1,2,2)
plot(mean_t_generalize,log10(meanerr_generalize),'Color',mediumblue,'LineWidth',3);hold on 
plot(mean_t_pd,log10(meanerr_pd),'Color',mediumgreen,'LineWidth',3);
plot(mean_t_pdQ,log10(meanerr_pdQ),'Color',mediumyellow,'LineWidth',3);
set(gca,'linewidth',3,'fontsize',18,'fontname','Times');
%legend('Generalized CP-PPA','Non-PPP','Fontsize',18,'Location','NorthEast');
legend('Generalized CP-PPA','Non-PPP,Q=\alpha A', 'Non-PPP,Q=-(\alpha+1)A','Fontsize',18,'Location','NorthEast');
set(gca,'linewidth',3,'fontsize',18,'fontname','Times');
xl = xlabel('time(s)'); 
yl = ylabel('log(Error)');
hold off


