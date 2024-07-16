data = zeros(3,10);
n = 500;
num_repeats = 10;
maxiter = 200;
for i = 1:10
    alpha = 0.1*i;
    for j = 1:3
        theta = 0.5*j;
        data(j,i) = exp_theta_alpha(num_repeats,maxiter,n,alpha,theta);
    end
end

save data data