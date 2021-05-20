function [rec_phase] = TIE(I_before,I_image,I_after,delta_z,k_0,graphs, compare, SLM_pixel, SLM_type)
%TIE will calculate the phase of an object using the Transport of Intensity
%equation. It requires 3 images at 3 focal planes, distanced delta_z from
%one another.
%k_0 is the wave number, and is equal to 2*pi/lambda
%spat_freqs is an array of the x and y spatial frequencies, k_x, k_y
%respectively.
%compare is the data of the laplacian received using the del2 function
%the graphs parameter will determine whether to draw graphs of not. 
if SLM_type == 1
    I_after = I_after(100:900,100:900);
    I_before = I_before(100:900,100:900);
    dIdz = (I_after - I_before)./(2*delta_z);           % approximate the derivative with respect to z axis
    I_avg = mean(mean(I_image(100:900,100:900)));       % avarage value of intensity at image plane
    I = -k_0*dIdz./I_avg;                              % Laplacian
else
    resolution = length(I_image);
    margin = 80;
    I_left = I_image(479 - margin:541 + margin, 394 - margin:444 + margin);                        % left part of image
    I_right = I_image(479 - margin:541 + margin, 569 - margin:619 + margin);                  %right part of image
    figure;
    imagesc(I_left)
    title('left')
    figure;
    imagesc(I_right);
    title('right')
    dIdz = (I_left - I_right)./(2*delta_z);                     % approximate the derivative with respect to z axis         
    I_avg = mean(mean(I_image(100:900,100:900)));       % avarage value of intensity at image plane
    I = -k_0*dIdz./I_avg;                              % Laplacian
    I(resolution,resolution) = 0;
    I = centerImage(I,size(dIdz));
    figure;
    imagesc(I)
    title('dIdz with SLM')
end
% set Dirichlet conditions on Laplacian, zeroing the edges
% dirich_border = 1;
% I(1:end,1:dirich_border) = 0;
% I(1:end,end-dirich_border:end) = 0;
% I(1:dirich_border,1:end) = 0;
% I(end-dirich_border:end,1:end) = 0; 
% expanded Dirichlet condition from only image edges to 100 first
% rows/columns. didnt help so much

% for some reason we get an image mirrored on the 45 deg axis. these bits
% of code are used to correct this when necessary
I = flip(I);
I = flip(I,2);

% normalize both results to better compare them
I_abs = abs(I);
comp_abs = abs(compare);
I_norm = I./max(I_abs(:));
comp_norm = compare./max(comp_abs(:));

if graphs
   
   figure;
   imagesc(I) 
   title('\nabla^2 using approximation')
   colorbar
   
   
   figure;
   imagesc(I_thresh) 
   title('\nabla^2 using approximation - threshold')
   colorbar
   

   figure;
   imagesc(compare) 
   title('\nabla^2 using del2')
   colorbar

end

if graphs
       
   figure;
   imagesc(comp_norm) 
   title('\nabla^2 using del2 - Normalized')
   colorbar
   
   figure;
   imagesc(I_norm) 
   title('\nabla^2 using approximation - Normalized')
   colorbar
   
   figure;
   imagesc(I_thresh_norm) 
   title('\nabla^2 using approximation - threshold + normalize')
   colorbar

end


% recipcoral sum of spatial freq's using nati's code:

N = length(I_image);
fsqr=repmat((1i*2*pi/(N*SLM_pixel))*(-(N-1)/2:(N-1)/2),N,1).^2+repmat((1i*2*pi/(N*SLM_pixel))*(-(N-1)/2:(N-1)/2)',1,N).^2; 
k_recip = 1./fsqr; 
k_recip(~isfinite(k_recip))=0; 



if graphs                           % graph of spatial frequencies
    figure;
    imagesc(k_recip)
    colorbar
    title('2D Spatial Frequencies')
end

rec_phase = ift2(k_recip .* ft2(I));                 % reconstructing phase
image_data = real(rec_phase);                       % get real part of reconstructed phase


slice = image_data(360:440, 400);


if ~graphs                           % graphs of reconstructed image
    
    figure;
    plot(slice)


    figureToSave = figure;
    imagesc(image_data);
    colorbar();
    title("Reconstructed data of image")
    figFileName = char(strcat("../Docs/images/", get(get(gca,'title'),'string'), ".jpg"));
    saveas(figureToSave, figFileName)
%     
%     figureToSave = figure;
%     imagesc(image_data_comp);
%     colorbar();
%     title("Reconstructed data of image for comparison")
%     figFileName = char(strcat("../Docs/images/", get(get(gca,'title'),'string'), ".jpg"));
%     saveas(figureToSave, figFileName)
%     
%     figureToSave = figure;
%     imagesc(image_data_thresh);
%     colorbar();
%     title("Reconstructed data of image with threshold")
%     figFileName = char(strcat("../Docs/images/", get(get(gca,'title'),'string'), ".jpg"));
%     saveas(figureToSave, figFileName)
%     
%     figureToSave = figure;
%     imagesc(image_data_thresh_norm);
%     colorbar();
%     title("Reconstructed data of image with threshold and normalization")
%     figFileName = char(strcat("../Docs/images/", get(get(gca,'title'),'string'), ".jpg"));
%     saveas(figureToSave, figFileName)
%     
end
end