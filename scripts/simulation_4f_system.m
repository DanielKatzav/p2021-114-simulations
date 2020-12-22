function [I_before_image_plane,I_image_plane, I_after_image_plane] = simulation_4f_system(image,lambda,distances,focus,res,SLM_pixel,delta_z,graphs)
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

x = -res/2:res/2-1;             % x axis span
y = -res/2:res/2-1;             % y axis span
[X,Y] = meshgrid(x*SLM_pixel,y*SLM_pixel);  %create x,y meshgrid
%% Object construction
% A = rgb2gray(imread('..\common\1951usaf_test_target.jpg'));    % load iamge

cropped_img = image(200:300,40:200);      % select specified rows and columns from imgae
cropped_img(res,res) = 0;                 % increase image size to 1000x1000
cropped_img = circshift(cropped_img,[450 450]); % center image elements
if graphs
    figureToSave = figure;
    imagesc(cropped_img)                  % show image
    title('USAF Resolution chart as object')
    figFileName = strcat("../Docs/images/", get(get(gca,'title'),'string'), ".jpg");
    saveas(figureToSave, figFileName)
end
phase_div = 2*pi/256;                     % divide 2pi to accomodate grayscale
phase = exp(1i*double(cropped_img)*phase_div);      % convert A to double and vreate phase object
if graphs
    figureToSave = figure;
    imagesc(abs(phase).^2)          % show phase object
    title('USAF Resolution chart as phase object')
    figFileName = strcat("../Docs/images/", get(get(gca,'title'),'string'), ".jpg");
    saveas(figureToSave, figFileName)
    
end
%%  Phase functions
f1 = focus(1);                % focus length of first lens
f2 = focus(2);                % focus length of second lens
z_o = distances(1);               % distance from object plane to first lens
z_i = distances(2);               % distance from second lebs to image plane
%% Complex Field Propagation of P
P_image_plane = propagation4f(phase,[z_o z_i],[f1 f2], lambda, X,Y,graphs,"at image plane");
P_after_image_plane = propagation4f(phase,[z_o z_i+delta_z],[f1 f2], lambda, X,Y,graphs, "after image plane");
P_before_image_plane = propagation4f(phase,[z_o z_i-delta_z],[f1 f2], lambda, X,Y,graphs, "before image plane ");
%% Camera        
I_image_plane = P_image_plane.*conj(P_image_plane);   % intensity of the image at imaging plane I = u*(u*)
I_after_image_plane = P_after_image_plane.*conj(P_after_image_plane);   % intensity of the image at imaging plane I = u*(u*)
I_before_image_plane = P_before_image_plane.*conj(P_before_image_plane);   % intensity of the image at imaging plane I = u*(u*)

if graphs
    figureToSave = figure;
    imagesc(I_image_plane)                                % show intensity of phase object1
    title('Intensity of phase object at image plane')
    figFileName = strcat("../Docs/images/", get(get(gca,'title'),'string'), ".jpg");
    saveas(figureToSave, figFileName)
    
    figureToSave = figure;
    imagesc(I_after_image_plane)                                % show intensity of phase object1
    title('Intensity of phase object after image plane')
    figFileName = strcat("../Docs/images/", get(get(gca,'title'),'string'), ".jpg");
    saveas(figureToSave, figFileName)
    
    figureToSave = figure;
    imagesc(I_before_image_plane)                                % show intensity of phase object1
    title('Intensity of phase object before image plane')
    figFileName = strcat("../Docs/images/", get(get(gca,'title'),'string'), ".jpg");
    saveas(figureToSave, figFileName)

end
end