function [propagated] = lensPropagate(A,f,lambda,X,Y,graphs, nameOfPlane)
% lensPropagate calculates the propagation of an EM wave through lens with
% focus length f. A is a 2D wave front complex values matrix, lambda is the
% wavelength, X,Y are the meshgrid values. The function uses fft and ifft
% to calculate propagation using Fourier optics transfer functions.

Qlens = exp(-1i*(pi/(lambda*f))*(X.^2+Y.^2));   % propagation transfer function of lens with focus f
propagated = A.*Qlens;                          % propagation through the lens
if graphs
   figureToSave = figure;
   imagesc(propagated.*conj(propagated))
   colorbar();
   title("intensity after lens propagation")
   figFileName = char(strcat("../Docs/images/", get(get(gca,'title'),'string'), nameOfPlane, ".jpg"));
   saveas(figureToSave, figFileName)
end
end