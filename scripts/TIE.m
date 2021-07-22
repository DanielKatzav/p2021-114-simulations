function [rec_phase] = TIE(I_before,I_image,I_after,delta_z,k_0,graphs, compare, SLM_pixel, SLM_type)
%TIE will calculate the phase of an object using the Transport of Intensity
%equation. It requires 3 images at 3 focal planes, distanced delta_z from
%one another.
%k_0 is the wave number, and is equal to 2*pi/lambda
%spat_freqs is an array of the x and y spatial frequencies, k_x, k_y
%respectively.
%compare is the data of the laplacian received using the del2 function
%the graphs parameter will determine whether to draw graphs of not. 
resolution = length(I_image);
if SLM_type == 1
    I_after = I_after(100:900,100:900);
    I_before = I_before(100:900,100:900);
    dIdz = (I_after - I_before)./(2*delta_z);           % approximate the derivative with respect to z axis
    I_avg = mean(mean(I_image(100:900,100:900)));       % avarage value of intensity at image plane
    I = -k_0*dIdz./I_avg;                              % Laplacian
    I(resolution,resolution) = 0;
else
    margin = 80;
    I_left = I_image(479 - margin:541 + margin, 394 - margin:444 + margin);                        % left part of image
    I_right = I_image(479 - margin:541 + margin, 569 - margin:619 + margin);                  %right part of image
    
    saveFigure(I_left,'Left Shifted Image at Image Plane','',graphs);
    saveFigure(I_right,'Right Shifted Image at Image Plane','',graphs);

    dIdz = (I_left - I_right)./(2*delta_z);                     % approximate the derivative with respect to z axis         
    I_avg = mean(mean(mean(I_left))+mean(mean(I_right)));       % avarage value of intensity at image plane
    I = -k_0*dIdz./I_avg;                              % Laplacian
    I(resolution,resolution) = 0;
    I = centerImage(I,size(dIdz));
end

I = flip(I);
I = flip(I,2);

   
saveFigure(I,'\nabla^2 Using TIE Method','',graphs);
saveFigure(compare,'\nabla^2 Using del2() Function','',graphs);

% recipcoral sum of spatial freq's using nati's code:

N = length(I_image);
fsqr=repmat((1i*2*pi/(N*SLM_pixel))*(-(N-1)/2:(N-1)/2),N,1).^2+repmat((1i*2*pi/(N*SLM_pixel))*(-(N-1)/2:(N-1)/2)',1,N).^2; 
k_recip = 1./fsqr; 
k_recip(~isfinite(k_recip))=0; 

saveFigure(k_recip,'2D Spatial Frequencies','',graphs);        % graph of spatial frequencies

rec_phase = ift2(k_recip .* ft2(I));                 % reconstructing phase
image_data = real(rec_phase);                       % get real part of reconstructed phase

% slice = image_data(400:600, 500);           % intersection slice of the reconstructed data
% 
% 
% saveFigure(slice,'Slice of Middle of Reconstruction','',graphs);

figureToSave = figure;
histogram(image_data(440:560,440:560));     % histogram of the reconstructed data
title('Histogram of reconstructed Data')
figFileName = char(strcat("../Docs/images/", get(get(gca,'title'),'string'), ".jpg"));
saveas(figureToSave, figFileName)

saveFigure(image_data,'Reconstructed data of image','',graphs);

end