clc
close all
clear variables

N = 1000;                   % resolution
lambda = 0.633e-6;          % wavelength
SLM_pixel = 8e-6;           % SLM pixel resolution
del = SLM_pixel;            % SLM pixel resolution - AGAIN!
x = -N/2:N/2-1;             % x axis span
y = -N/2:N/2-1;             % y axis span
[X,Y] = meshgrid(x*del,y*del);  %create x,y meshgrid
%% Object construction
A = rgb2gray(imread('..\common\1951usaf_test_target.jpg'));    % load iamge
figure;
imagesc(A)                  % show image
A = A(200:300,40:200);      % select specified rows and columns from imgae
A(N,N) = 0;                 % increase image size to 1000x1000
A = circshift(A,[450 450]); % center image elements

P = exp(1i*double(A));      % convert A to double and vreate phase object
figure;
imagesc(abs(P))                  % show phase
%%  Phase functions
f1 = 100e-3;                % focus length of first lens
f2 = 100e-3;                % focus length of second lens
z_o = 100e-3;               % distance from object plane to first lens
z_i = 100e-3;               % distance from second lebs to image plane
Qlens_f1 = exp(-1i*(pi/(lambda*f1))*(X.^2+Y.^2));   % first lens transfer function
Qlens_f2 = exp(-1i*(pi/(lambda*f2))*(X.^2+Y.^2));   % second lens transfer function
prop_z_o = exp(1i*(pi/(lambda*z_o))*(X.^2+Y.^2));   % propagation transfer function over z_o
prop_z_i = exp(1i*(pi/(lambda*z_i))*(X.^2+Y.^2));   % propagation transfer function over z_i
%% Complex Field Propagation of A
A_before_lens1 = ift2(ft2(A).*ft2(prop_z_o));           % image after propagation over z_o
A_after_lens1 = A_before_lens1.*Qlens_f1;               % image after propagation trough lens 1
A_before_lens2 = ift2(ft2(A_after_lens1).*ft2(exp(1i*(pi/(lambda*(f1+f2)))*(X.^2+Y.^2))));  %image after propagating between lenses
A_after_lens2 = A_before_lens2.*Qlens_f2;               % image after propagating through lens 2
A_image_plane = ift2(ft2(A_after_lens2).*ft2(prop_z_i));% image after propagationg over z_i
%% Complex Field Propagation of P
P_before_lens1 = ift2(ft2(P).*ft2(prop_z_o));           % image after propagation over z_o
P_after_lens1 = P_before_lens1.*Qlens_f1;               % image after propagation trough lens 1
P_before_lens2 = ift2(ft2(P_after_lens1).*ft2(exp(1i*(pi/(lambda*(f1+f2)))*(X.^2+Y.^2))));  %image after propagating between lenses
P_after_lens2 = P_before_lens2.*Qlens_f2;               % image after propagating through lens 2
P_image_plane = ift2(ft2(P_after_lens2).*ft2(prop_z_i));% image after propagationg over z_i
%% Camera
I_amp = A_image_plane.*conj(A_image_plane);         % intensity of the image at imaging plane I = u*(u*)
I_phase = P_image_plane.*conj(P_image_plane);
figure;
imagesc(I_amp)                               % show intensity of normal object
title('Intensity of normal object')
figure;
imagesc(I_phase)                               % show intensity of phase object1
title('Intensity of phase object')
