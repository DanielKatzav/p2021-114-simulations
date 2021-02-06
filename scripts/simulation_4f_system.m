function [I_before_image_plane,I_image_plane, I_after_image_plane, lapl] = simulation_4f_system(image,lambda,distances,focus,resolution,X, Y,delta_z,graphs)
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
threshold = 250;                                    % threshold for binary phase usage
cropped_img = image(230:270,60:90);                % select specified rows and columns from image
cropped_img = imresize(cropped_img, 20);
cropped_img(resolution,resolution) = 0;             % increase image size to 1000x1000
cropped_img = circshift(cropped_img,[85 205]);     % center image elements
binary_img = uint8((cropped_img >= threshold));     % create binary values depending on threshold

if graphs
    figureToSave = figure;
    imagesc(binary_img)                  % show image
    colorbar();
    title('USAF Resolution chart as object')
    figFileName = char(strcat("../Docs/images/", get(get(gca,'title'),'string'), ".jpg"));
    saveas(figureToSave, figFileName)
end

phase_const = pi/3;                                  % constant to multiply binary img, s.t. exp doesnt zero
phase = complex(exp(1i*double(binary_img)*phase_const));      % convert A to double and create phase object
if graphs
    figureToSave = figure;
    imagesc(angle(phase))          % show phase object
    colorbar();
    title('USAF Resolution chart as phase object')
    figFileName = char(strcat("../Docs/images/", get(get(gca,'title'),'string'), ".jpg"));
    saveas(figureToSave, figFileName)
    
end

%% Get laplacian using del2 to compare with 4f system
lapl = del2(phase);

if ~graphs
    figureToSave = figure;
    imagesc(abs(lapl))          % show laplacian object
    colorbar();
    title('Phase object Laplacian using del2')
    figFileName = char(strcat("../Docs/images/", get(get(gca,'title'),'string'), ".jpg"));
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
    colorbar();
    title('Intensity of phase object at image plane ')
    figFileName = char(strcat("../Docs/images/", get(get(gca,'title'),'string'), ".jpg"));
    saveas(figureToSave, figFileName)
    
    figureToSave = figure;
    imagesc(I_after_image_plane)                                % show intensity of phase object1
    colorbar();
    title('Intensity of phase object after image plane ')
    figFileName = char(strcat("../Docs/images/", get(get(gca,'title'),'string'), ".jpg"));
    saveas(figureToSave, figFileName)
    
    figureToSave = figure;
    imagesc(I_before_image_plane)                                % show intensity of phase object1
    colorbar();
    title('Intensity of phase object before image plane ')
    figFileName = char(strcat("../Docs/images/", get(get(gca,'title'),'string'), ".jpg"));
    saveas(figureToSave, figFileName)

end
end