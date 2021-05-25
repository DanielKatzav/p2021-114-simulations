function [I_before_image_plane,I_image_plane, I_after_image_plane, lapl] = simulation_4f_system(image,lambda,distances,focus,resolution,X, Y,delta_z,graphs,SLM_type)
%simulation_4f_system will simulate a 4f system with the following
%parameters:
%A - matrix of a loaded image in grayscale.
%lambda - simulated laser wavelength in meters
%distances - 1x2 array of distances from object to first lens and from
%focus - 1x2 arrat of focus lengths of lens1 and lens2, respectively.
%res - image resolution, integer that has to be greater that loaded image
% resolution. Final resolution will be res x res size.
%SLM_pixel - SLM pixel resolution in meters
%delta_z - distance to move the image plane out of focus in meters
% second lens to image, respectively.
%graphs - Boolean to determine drawing graphs. True - draw graphs. False -
% dont.


%% Object construction
h_start = 760;                      % height starting point for cropping
h_end = 860;                        % height end point for cropping
w_start = 30;                       % width starting point for cropping
w_end = 100;                         % width end point for cropping
resize_factor = 1;                 % resize factor for scaling image


threshold = 240;                                    % threshold for binary phase usage
cropped_img = image(h_start:h_end,w_start:w_end);   % select specified rows and columns from image
cropped_img = imresize(cropped_img, resize_factor);
cropped_img_size = size(cropped_img);
cropped_img(resolution,resolution) = 0;             % increase image size to 1000x1000
cropped_img = centerImage(cropped_img, cropped_img_size);
binary_img = uint8((cropped_img >= threshold));     % create binary values depending on threshold

saveFigure(binary_img,'USAF Resolution chart as object','',graphs);

phase_const = pi/3;                                  % constant to multiply binary img, s.t. exp doesnt zero
phase_obj = complex(exp(1i*double(binary_img)*phase_const));      % convert A to double and create phase object

mask = ones(126,113);                               % create a mask to remove noise in TIE
mask_size = size(mask);
mask(resolution,resolution) = 0;
mask = centerImage(mask, mask_size);

phase_obj = phase_obj.*mask;                        % mask phase object to clear shift noise

saveFigure(angle(phase_obj),'Phase of USAF Resolution chart as phase object','',graphs);
saveFigure(abs(phase_obj),'Amplitude of USAF Resolution chart as phase object','',graphs);


%% Get laplacian using del2 to compare with 4f system
lapl = del2(double(binary_img)*phase_const);
saveFigure(abs(lapl),'Phase object Laplacian using del2() function','',graphs);


%%  Phase functions
f1 = focus(1);                % focus length of first lens
f2 = focus(2);                % focus length of second lens
z_o = distances(1);               % distance from object plane to first lens
z_i = distances(2);               % distance from second lebs to image plane

f_SLM = focus(2)^2/delta_z;
%% Complex Field Propagation of P

P_image_plane = propagation4f(phase_obj,[z_o z_i],[f1 f2],f_SLM, lambda, X,Y,graphs,"at Image Plane", SLM_type);
if SLM_type == 2
    graphs = false;         %disable graph printing of defocused planes while muxing
end
P_after_image_plane = propagation4f(phase_obj,[z_o z_i+delta_z],[f1 f2],f_SLM, lambda, X,Y,graphs, "After Image Plane", SLM_type);
P_before_image_plane = propagation4f(phase_obj,[z_o z_i-delta_z],[f1 f2],f_SLM, lambda, X,Y,graphs, "Before Image Plane ", SLM_type);
graphs = true;
%% Camera        
I_image_plane = P_image_plane.*conj(P_image_plane);   % intensity of the image at imaging plane I = u*(u*)
I_after_image_plane = P_after_image_plane.*conj(P_after_image_plane);   % intensity of the image at imaging plane I = u*(u*)
I_before_image_plane = P_before_image_plane.*conj(P_before_image_plane);   % intensity of the image at imaging plane I = u*(u*)

saveFigure(I_image_plane,'Intensity of phase object at image plane ','',graphs);
if SLM_type == 1
    saveFigure(I_after_image_plane,'Intensity of phase object after image plane ','',graphs);
    saveFigure(I_before_image_plane,'Intensity of phase object before image plane ','',graphs);
end
end