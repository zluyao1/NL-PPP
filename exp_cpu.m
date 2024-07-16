function data_CPU = exp_cpu(v,n,maxiter,erro,num_repeats)

m = n/4; sparse = n/20; alpha = 0.5; theta = 1;

S = @(z,lambda) max(abs(z)-lambda,0).*sign(z);

err_generalized = zeros(maxiter,num_repeats);
err_pd = zeros(maxiter,num_repeats);
err_pdQ = zeros(maxiter,num_repeats);

cpu_generalized = zeros(num_repeats,1);
it_generalized = zeros(num_repeats,1);
cpu_pd = zeros(num_repeats,1);
it_pd = zeros(num_repeats,1);
cpu_pdQ = zeros(num_repeats,1);
it_pdQ = zeros(num_repeats,1);

rng(69462991)
for repeats = 1:num_repeats
    % Set up instance
    rng(69462991)
    A = randn(m,n);
    xhat = sparserandn(n,sparse);  % true solution
    b = A*xhat;
    L = norm(A'*A,2);

    r = sqrt(L*(1-alpha+alpha^2))/v; 
    s = sqrt(L*(1-alpha+alpha^2))*v;

    % gerneralized CP-PPA
    x = zeros(n,1); y = zeros(m,1);
    tic
    for k = 2:maxiter
        % gerneralized CP-PPA
        xn = S(x+A'*y/r,1/r); % x_{k+1}
        xbar = xn + alpha*(xn-x);
        ybar = y - (A*xbar-b)/s;
        yn = ybar - (1-alpha)*A*(xn-x)/s; % y_{k+1}
        % err_generalize(k,repeats) = norm(xn-x)/norm(xn);
        err_generalized(k,repeats) = norm(xn-xhat)/norm(xhat);
        x = xn; y = yn;

        if err_generalized(k,repeats) < erro
            cpu_generalized(repeats) = toc;
            it_generalized(repeats) = k;
            break;
        else
            cpu_generalized(repeats) = toc;
            it_generalized(repeats) = k;
        end
    end

    % Non_ppA, Q=\alpha A
    x = zeros(n,1); y = zeros(m,1);
    Q = (alpha)*A;

    tic
    for k = 1:maxiter

        xn = S(x+A'*y/r,1/r); % x_{k+1}
        yn = y - ((A-Q)*x+Q*xn-b)/s;
        an = [r*(x-xn)+A'*(y-yn);(Q-A)*(x-xn)+s*(y-yn)];
        un = [x;y];
        rn = [xn;yn];
        un1 = un - theta*(un-rn)'*an*an/norm(an,2)^2;

        err_pd(k,repeats) = norm(un1(1:n)-xhat)/norm(xhat);
        x = un1(1:n); y = un1((n+1):(n+m));

        if err_pd(k,repeats) < erro
            cpu_pd(repeats) = toc;
            it_pd(repeats) = k;
            break;
        else
            cpu_pd(repeats) = toc;
            it_pd(repeats) = k;
        end

    end

    % Non_ppA, Q=-(\alpha+1) A
    x = zeros(n,1); y = zeros(m,1);
    Q = -(1+alpha)*A;

    tic
    for k = 1:maxiter

        xn = S(x+A'*y/r,1/r); % x_{k+1}
        yn = y - ((A-Q)*x+Q*xn-b)/s;
        an = [r*(x-xn)+A'*(y-yn);(Q-A)*(x-xn)+s*(y-yn)];
        un = [x;y];
        rn = [xn;yn];
        un1 = un - theta*(un-rn)'*an*an/norm(an,2)^2;

        err_pdQ(k,repeats) = norm(un1(1:n)-xhat)/norm(xhat);
        x = un1(1:n); y = un1((n+1):(n+m));

        if err_pdQ(k,repeats) < erro
            cpu_pdQ(repeats) = toc;
            it_pdQ(repeats) = k;
            break;
        else
            cpu_pdQ(repeats) = toc;
            it_pdQ(repeats) = k;
        end

    end
end  % end for loop over repeats

data_CPU = struct('cpu_generalized',cpu_generalized,'it_generalized',it_generalized,...
    'cpu_pd',cpu_pd,'it_pd',it_pd,...
    'cpu_pdQ',cpu_pdQ,'it_pdQ',it_pdQ);

    %res_esrk
    %semilogy(1:maxiter,res_esrk)
    %pause
end