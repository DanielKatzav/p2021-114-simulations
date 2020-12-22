function [propagated] = distancePropagate(A,d,lambda,X,Y,graphs, nameOfPlane)
% distancePropagate calculates the propagation of an EM wave over distance
% d in air. A is a 2D wave front complex values matrix, lambda is the
% wavelength, X,Y are the meshgrid values. The function uses fft and ifft
% to calculate propagation using Fourier optics transfer functions.

prop_d = exp(1i*(pi/(lambda*d))*(X.^2+Y.^2));   % propagation transfer function over d
propagated = ift2(ft2(A).*ft2(prop_d));         % propagation through the distance
if graphs
   figureToSave = figure;
   imagesc(propagated.*conj(propagated))
   title("intensity after distance propagation")
    figFileName = strcat("../Docs/images/", get(get(gca,'title'),'string'), nameOfPlane, ".jpg");
    saveas(figureToSave, figFileName)
end

end