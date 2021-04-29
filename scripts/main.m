clc
close all
clear variables

%% Load image and set parameters
image = rgb2gray(imread('..\common\USAF-1951.png'));    % load iamge
image = 255 - image;            % invert colors (use on new USAF-1951.pgn)
lambda = 0.633e-6;              % wavelength
k_0 = 2*pi/lambda;              % wave number
resolution = 1000;              % resolution
SLM_pixel = 8e-6;               % SLM pixel resolution
delta_z = 500e-6;                 % distance to move the image plane out of focus
focus_dist = 100e-3;
distances = [focus_dist focus_dist];    % distances z_o and z_i
focus = [focus_dist focus_dist];        % focus lengths f1 and f2
graphs = false;                  % true - draw graphs, false - dont
SLM_type = 2;                   % to be updated - determine type of mask to put on SLM
x = -resolution/2:resolution/2-1;                         % x axis span
y = -resolution/2:resolution/2-1;                         % y axis span
[X,Y] = meshgrid(x*SLM_pixel,y*SLM_pixel);  % create x,y meshgrid
%% Simulate 3 focus planes in a 4f system
[I_before, I_image, I_after, lapl] = simulation_4f_system(image,lambda,distances,focus,resolution,X,Y,delta_z,graphs, SLM_type);
%results are the intensities of the image before, at and after the focus
%plane, respectively.

%% Calculate phase using TIE
phase = TIE(I_before,I_image,I_after,delta_z,k_0,graphs, lapl, SLM_pixel, SLM_type);
