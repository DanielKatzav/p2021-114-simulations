clc
close all
clear variables

N = 1000;                   % resolution
lambda = 0.633e-6;          % wavelength
SLM_pixel = 8e-6;           % SLM pixel resolution
del = SLM_pixel;            % SLM pixel resolution - AGAIN!
delta_z = 10e-3;            % distance to move the image plane out of focus
x = -N/2:N/2-1;             % x axis span
y = -N/2:N/2-1;             % y axis span
[X,Y] = meshgrid(x*del,y*del);  %create x,y meshgrid
%% Object construction
A = rgb2gray(imread('..\common\1951usaf_test_target.jpg'));    % load iamge
figure;
title('USAF Resolution chart as regular object')
imagesc(A)                  % show image
A = A(200:300,40:200);      % select specified rows and columns from imgae
A(N,N) = 0;                 % increase image size to 1000x1000
A = circshift(A,[450 450]); % center image elements

P = exp(1i*double(A));      % convert A to double and vreate phase object
figure;
title('USAF Resolution chart as phase object')
imagesc(abs(P))                  % show phase object
%%  Phase functions
f1 = 100e-3;                % focus length of first lens
f2 = 100e-3;                % focus length of second lens
z_o = 100e-3;               % distance from object plane to first lens
z_i = 100e-3;               % distance from second lebs to image plane
%% Complex Field Propagation of P
P_image_plane = propagation4f(P,[z_o z_i],[f1 f2], lambda, X,Y);
P_after_image_plane = propagation4f(P,[z_o z_i+delta_z],[f1 f2], lambda, X,Y);
P_before_image_plane = propagation4f(P,[z_o z_i-delta_z],[f1 f2], lambda, X,Y);
%% Camera        
I_phase = P_image_plane.*conj(P_image_plane);   % intensity of the image at imaging plane I = u*(u*)
figure;
imagesc(I_phase)                                % show intensity of phase object1
title('Intensity of phase object at image plane')

I_phase_after = P_after_image_plane.*conj(P_after_image_plane);   % intensity of the image at imaging plane I = u*(u*)
figure;
imagesc(I_phase_after)                                % show intensity of phase object1
title('Intensity of phase object after image plane')

I_phase_before = P_before_image_plane.*conj(P_before_image_plane);   % intensity of the image at imaging plane I = u*(u*)
figure;
imagesc(I_phase_before)                                % show intensity of phase object1
title('Intensity of phase object before image plane')
