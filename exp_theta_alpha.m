function value = exp_theta_alpha(num_repeats,maxiter,n,alpha,theta)
v = 10; 
m = n/4; sparse = n/20;

S = @(z,lambda) max(abs(z)-lambda,0).*sign(z);

err_pd = zeros(maxiter,num_repeats);

rng(69462991)
for repeats = 1:num_repeats
    % Set up instance
    A = randn(m,n);
    xhat = sparserandn(n,sparse);  % true solution
    b = A*xhat;
    L = norm(A'*A,2);

    r = sqrt(L*(1-alpha+alpha^2))/v; 
    s = sqrt(L*(1-alpha+alpha^2))*v;

    % Initialize methods
    x = zeros(n,1); y = zeros(m,1);
    Q = (alpha)*A;
    for k = 1:maxiter

        xn = S(x+A'*y/r,1/r); % x_{k+1}
        yn = y - ((A-Q)*x+Q*xn-b)/s;
        an = [r*(x-xn)+A'*(y-yn);(Q-A)*(x-xn)+s*(y-yn)];
        un = [x;y];
        rn = [xn;yn];
        un1 = un - theta*(un-rn)'*an*an/norm(an,2)^2;

        % err_pd(k,repeats) = norm(un1(1:n)-x)/norm(un1(1:n));
        err_pd(k,repeats) = norm(un1(1:n)-xhat)/norm(xhat);
        x = un1(1:n); y = un1((n+1):(n+m));

    end
end  % end for loop over repeats
value = median(err_pd(maxiter,:),2); 
end