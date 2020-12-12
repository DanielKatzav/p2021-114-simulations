function [phase] = TIE(I_before,I_image,I_after,delta_z,k_0,spat_freqs,graphs)
%TIE will calculate the phase of an object using the Transport of Intensity
%equation. It requires 3 images at 3 focal planes, distanced delta_z from
%one another.
%k_0 is the wave number, and is equal to 2*pi/lambda
%spat_freqs is an array of the x and y spatial frequencies, k_x, k_y
%respectively.
%the graphs parameter will determine whether to draw graphs of not. 
dIdz = (I_after - I_before)./(2*delta_z);           % approximate the derivative with respect to z axis
I = k_0*dIdz./I_image;                            % Fourier transform argument
k_recip = 1/(spat_freqs(1)^2 + spat_freqs(2)^2);    % reciprocal of sum of spatial freq's
phase = ift2(k_recip * ft2(I));                     % reconstructing phase

if graphs
    rec_image = log(phase);
    abs_image = rec_image.*conj(rec_image);
    real_image = real(rec_image);
    imag_image = imag(rec_image);
    figure;
    imagesc(abs_image);
    title("Reconstructed abs of image")
    figure;
    imagesc(real_image);
    title("Reconstructed Real of image")
    figure;
    imagesc(imag_image);
    title("Reconstructed Imaginary of image")
    figure;
    imagesc(phase.*conj(phase));
    title("Reconstructed intensity of image")
    
end
end