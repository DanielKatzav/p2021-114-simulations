function [rec_phase] = TIE(I_before,I_image,I_after,delta_z,k_0,k_x, k_y,graphs)
%TIE will calculate the phase of an object using the Transport of Intensity
%equation. It requires 3 images at 3 focal planes, distanced delta_z from
%one another.
%k_0 is the wave number, and is equal to 2*pi/lambda
%spat_freqs is an array of the x and y spatial frequencies, k_x, k_y
%respectively.
%the graphs parameter will determine whether to draw graphs of not. 
dIdz = (I_after - I_before)./(2*delta_z);           % approximate the derivative with respect to z axis
I = k_0*dIdz./I_image;                              % Fourier transform argument
k_recip = 1./(k_x.^2 + k_y.^2);                     % reciprocal of sum of spatial freq's
rec_phase = ift2(k_recip * ft2(I));                 % reconstructing phase
image_data = (angle(rec_phase)+pi)*256/(2*pi);

if ~graphs
    figureToSave = figure;
    imagesc(image_data);
    colorbar();
    title("Reconstructed data of image")
    figFileName = char(strcat("../Docs/images/", get(get(gca,'title'),'string'), ".jpg"));
    saveas(figureToSave, figFileName)
end
end