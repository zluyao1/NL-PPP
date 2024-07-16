% a = imread('fig1.bmp');
% whos a
% b = im2double(a);
% whos b
%% method 1
data_type = 1;  % 1 = gaussian noise, 2 = salt & pepper noise, 
% 3 = real image and guassian noise, 4 = real image and pepper noise

% clear_img = imread('fig3.bmp'); % fig1, fig2, fig3, fig4
% clear_img = im2double(clear_img);
% % add noise
% img = clear_img + sigma * randn(size(clear_img));
[clear_img, img] = gen_data(data_type, 0.1);
%% method 2
clear;clc;
%I = imread('Cameraman.tif');
I = imread('Peppers.tif');
% clear_img = imread('fig3.bmp'); % fig1, fig2, fig3
% I = im2double(clear_img);
doubleI = double(I);
doubleI = doubleI(:);
clear_img = im2double(I);
sigma = 0.1; % 0.05, 0.1
img = clear_img + sigma * randn(size(clear_img));
num_steps = 100;

[H, W] = size(img);
N = H * W;

%% precomputed
nabla = make_derivatives_mine(H, W);
divop = nabla';
% divop = make_divop(H, W);

tol_pred = zeros(1, num_steps);
tol_prim = zeros(1, num_steps);
tol_gener = zeros(1, num_steps);
toll_pred = zeros(1, num_steps);
toll_prim = zeros(1, num_steps);
toll_gener = zeros(1, num_steps);
criterion_pred = zeros(1, num_steps);
criterion_prim = zeros(1, num_steps);
criterion_gener = zeros(1, num_steps);

%% predict-correct 
u = img(:);
p = zeros(N * 2, 1);
L = sqrt(8);

alpha = 0.5; 
theta = 2;
% v = 1;
% r = sqrt(L*(1-alpha+alpha^2))/v; 
% s = sqrt(L*(1-alpha+alpha^2))*v;
r = 0.01;
s = 1/(L^2*r);
lambda = 16;

Q = (alpha)*nabla;
%Q = (alpha + sqrt(3)/2 + 1.5)*nabla;%1.5
%Q = sqrt(3)*sigma*divop;
u_tilde = zeros(N,1);
for step = 1:num_steps

    % predict-correct
    u_tilde = ( r*u - divop*p + lambda*img(:) )/(lambda+r);
    
    pp = p + (nabla * u + Q* (u-u_tilde))/s;
    p_len = sqrt(pp(1:N).^2 + pp(N+1:end).^2);
    p_len = max(1, p_len);
    p_tilde = pp ./ repmat(p_len, 2, 1);

    a = [r*(u-u_tilde) - divop*(p-p_tilde);
        (Q+nabla)*(u-u_tilde) + s*(p-p_tilde)];
    uk = [u;p];
    rk = [u_tilde;p_tilde];
    uk = uk - theta*(uk-rk)'*a*a/norm(a,2)^2;

    u = uk(1:N);
    p = uk((N+1):end);

    toll_pred(step) = norm(u-clear_img(:))/norm(clear_img(:));

%         pp = p + tau * nabla * u;
%         p_len = sqrt(pp(1:N).^2 + pp(N+1:end).^2);
%         p_len = max(1, p_len);
%         p_tilde = pp ./ repmat(p_len, 2, 1);
%         p_bar = - theta * p + (1 + theta) * p_tilde;
% 
%         denom = 1 + sigma * lambda;
%         u_tilde = (u - sigma * divop * p_bar + sigma * lambda * img(:)) / denom;
% 
%         u_new = u_tilde - theta*sigma*divop*(p-p_tilde);
%         tol_pred(step) = norm(u_new-u)/norm(u_new);
%         toll_pred(step) = norm(u_new-doubleI)/norm(doubleI);
%         p = p_tilde - tau*nabla*(u-u_tilde);
% 
%         u = u_new;
% 
%         criterion_pred(step) = Fval(u, clear_img, 1, 0) + Gval(u, img, 0);

end
out_img_pred = u;

% primal-dual
% L = sqrt(8);
% tau = 0.01;
% sigma = 1/(tau*L^2);
% alpha = 1;%0.5
% lambda = 16;
% u = img(:);
% p = zeros(N * 2, 1);
% for step = 1:num_steps
%     
%     pp = p + tau * nabla * u;
%     p_len = sqrt(pp(1:N).^2 + pp(N+1:end).^2);
%     p_len = max(1, p_len);
%     p_tilde = pp ./ repmat(p_len, 2, 1);
%     p_bar = - alpha * p + (1 + alpha) * p_tilde;
%     p = p_tilde;
% 
%     denom = 1 + sigma * lambda;
%     u_new = (u - sigma * divop * p_bar + sigma * lambda * img(:)) / denom;
%     
%     tol_prim(step) = norm(u_new-u)/norm(u_new);
%     toll_prim(step) = norm(u_new-doubleI)/norm(doubleI);
%     u = u_new;
%     % compute the criterion function value
%     criterion_prim(step) = Fval(u, clear_img, 1, 0) + Gval(u, img, 0);
% 
% end
% out_img_prim = u;

% Generalized primal-dual
L = sqrt(8);

% tau = 0.01;
% sigma = 1/(tau*L^2);
% alpha = 0.5;%0.5

alpha = 0.5; theta = 1; v = 1;
% r = sqrt(L*(1-alpha+alpha^2))/v; 
% s = sqrt(L*(1-alpha+alpha^2))*v;
r = 0.01;
s = 1/(L^2*r);

lambda = 16;
u = img(:);
p = zeros(N * 2, 1);
u_bar = zeros(N,1);
for step = 1:num_steps
    
%     pp = p + tau * nabla * u;
%     p_len = sqrt(pp(1:N).^2 + pp(N+1:end).^2);
%     p_len = max(1, p_len);
%     p_tilde = pp ./ repmat(p_len, 2, 1);
%     p_bar = - alpha * p + (1 + alpha) * p_tilde;
% 
%     denom = 1 + sigma * lambda;
%     u_new = (u - sigma * divop * p_bar + sigma * lambda * img(:)) / denom;
%     u_new = u_new - (1-alpha)*sigma*divop*(p_tilde-p);
% 
%     p = p_tilde;
% 
%     tol_gener(step) = norm(u_new-u)/norm(u_new);
%     toll_gener(step) = norm(u_new-doubleI)/norm(doubleI);
%     u = u_new;
    u_new = (r*u - divop*p + lambda*img(:))/(lambda+r);
    u_bar = (1+alpha)*u_new-alpha*u;

    pp = p + nabla * u_bar/s;
    p_len = sqrt(pp(1:N).^2 + pp(N+1:end).^2);
    p_len = max(1, p_len);
    p_bar = pp ./ repmat(p_len, 2, 1);

    p = p_bar - (1-alpha)*nabla*(u_new-u)/s;
    u = u_new;

    toll_gener(step) = norm(u-clear_img(:))/norm(clear_img(:));

end

out_img_gener = u;


subplot(2,2, 1);
imshow(clear_img);
title(' ');
set(gca,'linewidth',3,'fontsize',15,'fontname','Times');

subplot(2,2, 2);
PSNR = psnr(img,reshape(clear_img,H,W));
imshow(img);
title([sprintf('PSNR = %.4f',PSNR)]);
set(gca,'linewidth',3,'fontsize',15,'fontname','Times');

% subplot(2,2, 2);
% PSNR = psnr(reshape(median(out_img_prim,2),H,W),reshape(clear_img,H,W));
% imshow(reshape(out_img_prim, size(clear_img)));
% title('primal-dual');
% title([sprintf('PSNR = %.4f',PSNR)]);
% set(gca,'linewidth',3,'fontsize',15,'fontname','Times');

subplot(2,2, 3);
PSNR = psnr(reshape(median(out_img_gener,2),H,W),reshape(clear_img,H,W));
imshow(reshape(out_img_gener, size(clear_img)));
title('primal-dual');
title([sprintf('PSNR = %.4f',PSNR)]);
set(gca,'linewidth',3,'fontsize',15,'fontname','Times');

subplot(2,2, 4);
PSNR = psnr(reshape(median(out_img_pred,2),H,W),reshape(clear_img,H,W));
imshow(reshape(out_img_pred, size(clear_img)));
title('predict-correct');
title([sprintf('PSNR = %.4f',PSNR)]);
set(gca,'linewidth',3,'fontsize',15,'fontname','Times');

%%
% subplot(1,1,1)
% mediumgreen = [0.4660 0.6740 0.1880];
% mediumblue = [0 0.4470 0.7410];
% semilogy(1:num_steps,toll_gener,'Color',mediumblue,'LineWidth',2); hold on;
% semilogy(1:num_steps,toll_pred,'Color',mediumgreen,'LineWidth',2);
% legend('Generalized CP-PPA','Non-PPP');
% xlabel('Iteration','FontName','Times New Roman','FontSize',25)
% ylabel('Error','FontName','Times New Roman','FontSize',25)
% set(gca,'linewidth',3,'fontsize',18,'fontname','Times');