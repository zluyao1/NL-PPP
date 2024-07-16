
v = 10; n = 5000; erro = 1e-6; maxiter = 2000; num_repeats = 5;
data = exp_cpu(v,n,maxiter,erro,num_repeats);
med_cpu = [median(data.cpu_generalized),median(data.cpu_pd),median(data.cpu_pdQ)]

med_it = [median(data.it_generalized),median(data.it_pd),median(data.it_pdQ)]