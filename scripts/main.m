clc
close all
clear variables

%% Load image and set parameters
image = rgb2gray(imread('..\common\1951usaf_test_target.jpg'));    % load iamge
lambda = 0.633e-6;              % wavelength
k_0 = 2*pi/lambda;              % wave number
res = 1000;                     % resolution
SLM_pixel = 8e-6;               % SLM pixel resolution
delta_z = 10e-3;                % distance to move the image plane out of focus
distances = [100e-3 100e-3];    % distances z_o and z_i
focus = [100e-3 100e-3];        % focus lengths f1 and f2
graphs = false;                  % true - draw graphs, false - dont
SLM_type = 1;                   % to be updated - determine type of mask to put on SLM

%% Simulate 3 focus planes in a 4f system
[I_before, I_image, I_after] = simulation_4f_system(image,lambda,distances,focus,res,SLM_pixel,delta_z,graphs);
%results are the intensities of the image before, at and after the focus
%plane, respectively.

%% Calculate phase using TIE
[k_x, k_y] = calc_spatial_freqs(k_0, SLM_type, I_image);
spat_freqs = [k_x, k_y];
phase = TIE(I_before,I_image,I_after,delta_z,k_0,spat_freqs,graphs);
