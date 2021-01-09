function [rec_phase] = TIE(I_before,I_image,I_after,delta_z,k_0,k_x, k_y,graphs)
%TIE will calculate the phase of an object using the Transport of Intensity
%equation. It requires 3 images at 3 focal planes, distanced delta_z from
%one another.
%k_0 is the wave number, and is equal to 2*pi/lambda
%spat_freqs is an array of the x and y spatial frequencies, k_x, k_y
%respectively.
%the graphs parameter will determine whether to draw graphs of not. 
dIdz = (I_after - I_before)./(2*delta_z);           % approximate the derivative with respect to z axis
I = k_0*dIdz./I_image;                              % Laplacian
% set Dirichlet conditions on Laplacian, zeroing the edges
I(1:end,1) = 0;
I(1:end,end) = 0;
I(1,1:end) = 0;
I(end,1:end) = 0; 


I_del = del2(k_0.*I_image);
I_del(1:end,1) = 0;
I_del(1:end,end) = 0;
I_del(1,1:end) = 0;
I_del(end,1:end) = 0; 
I_del = I_del./1e15;


if ~graphs
   figure;
   imagesc(I) 
   title('\nabla^2 using approximation')
   colorbar
   figure;
   imagesc(I_del);
   title('\nabla^2 using del2 function')
   colorbar
end

% k_recip = 1./(k_x.^2 + k_y.^2);                     % reciprocal of sum of spatial freq's
N = length(I);
fsqr=repmat((1i*2*pi/(2*N))*(-(N-1)/2:(N-1)/2),N,1).^2+repmat((1i*2*pi/(2*N))*(-(N-1)/2:(N-1)/2)',1,N).^2;% 
k_recip = 1./fsqr; 
k_recip(~isfinite(k_recip))=0; 

rec_phase = ift2(k_recip * ft2(I));                 % reconstructing phase
image_data = real(rec_phase);
rec_phase_del = ift2(k_recip * ft2(I_del));
image_data_del = real(rec_phase_del);
if graphs
    figure;
    imagesc(k_recip)
    colorbar
    title('2D Spatial Frequencies')
end

if ~graphs
    figureToSave = figure;
    imagesc(image_data);
    colorbar();
    title("Reconstructed data of image")
    figFileName = char(strcat("../Docs/images/", get(get(gca,'title'),'string'), ".jpg"));
    saveas(figureToSave, figFileName)
    
    figureToSave = figure;
    imagesc(image_data_del);
    colorbar();
    title("Reconstructed data of image from del2")
    figFileName = char(strcat("../Docs/images/", get(get(gca,'title'),'string'), ".jpg"));
    saveas(figureToSave, figFileName)
end
end