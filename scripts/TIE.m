function [rec_phase] = TIE(I_before,I_image,I_after,delta_z,k_0,k_x, k_y,graphs, compare)
%TIE will calculate the phase of an object using the Transport of Intensity
%equation. It requires 3 images at 3 focal planes, distanced delta_z from
%one another.
%k_0 is the wave number, and is equal to 2*pi/lambda
%spat_freqs is an array of the x and y spatial frequencies, k_x, k_y
%respectively.
%compare is the data of the laplacian received using the del2 function
%the graphs parameter will determine whether to draw graphs of not. 
dIdz = (I_after - I_before)./(2*delta_z);           % approximate the derivative with respect to z axis
I = -k_0*dIdz./I_image;                              % Laplacian
% set Dirichlet conditions on Laplacian, zeroing the edges
I(1:end,1:100) = 0;
I(1:end,end-100:end) = 0;
I(1:100,1:end) = 0;
I(end-100:end,1:end) = 0; 
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

del_difference = comp_norm - I_norm;

if ~graphs
   figure;
   imagesc(I) 
   title('\nabla^2 using approximation')
   colorbar
   
   figure;
   imagesc(I_norm) 
   title('\nabla^2 using approximation - Normalized')
   colorbar
   
   figure;
   imagesc(compare) 
   title('\nabla^2 using del2')
   colorbar
   
   figure;
   imagesc(comp_norm) 
   title('\nabla^2 using del2 - Normalized')
   colorbar
   
   figure;
   imagesc(del_difference) 
   title('Difference between approx. and del2 \nabla^2 (normalized)')
   colorbar


end

% k_recip = 1./(k_x.^2 + k_y.^2);% reciprocal of sum of spatial freq's

% recipcoral sum of spatial freq's using nati's code:

N = length(I);
fsqr=repmat((1i*2*pi/(2*N))*(-(N-1)/2:(N-1)/2),N,1).^2+repmat((1i*2*pi/(2*N))*(-(N-1)/2:(N-1)/2)',1,N).^2; 
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

rec_phase_comp = ift2(k_recip .* ft2(compare));
image_data_comp = real(rec_phase_comp);



if ~graphs                           % graphs of reconstructed image
    figureToSave = figure;
    imagesc(image_data);
    colorbar();
    title("Reconstructed data of image")
    figFileName = char(strcat("../Docs/images/", get(get(gca,'title'),'string'), ".jpg"));
    saveas(figureToSave, figFileName)
    
    figureToSave = figure;
    imagesc(image_data_comp);
    colorbar();
    title("Reconstructed data of image for comparison")
    figFileName = char(strcat("../Docs/images/", get(get(gca,'title'),'string'), ".jpg"));
    saveas(figureToSave, figFileName)
    
end
end