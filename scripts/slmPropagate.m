function [propagated] = slmPropagate(A,f,lambda,X,Y,graphs, nameOfPlane, SLM_type)
% slmPropagate calculates the propagation of an EM wave through a specified 
% SLM, using the parameters of such SLM.
x_shift = -100;
y_shift = 0;
% Qlens = exp(-1i*(pi/(lambda*f))*(X.^2+Y.^2));   % propagation transfer function of lens with focus f
Qshift = exp(1i*(pi/lambda)*(X*x_shift + Y*y_shift)); % transfer function of linear phase shift
Qlens = 1;
%TODO: check if such function is correct!
if SLM_type == 1
    propagated = A;
else
    propagated = A.*Qlens.*Qshift;
end

if graphs
   figureToSave = figure;
   imagesc(propagated.*conj(propagated))
   colorbar();
   title("intensity after SLM trnasfer ")
    figFileName = char(strcat("../Docs/images/", get(get(gca,'title'),'string'), nameOfPlane, ".jpg"));
    saveas(figureToSave, figFileName)
end
end
