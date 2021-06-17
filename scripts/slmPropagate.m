function [propagated] = slmPropagate(A,f,lambda,X,Y,graphs, nameOfPlane, SLM_type)
% slmPropagate calculates the propagation of an EM wave through a specified 
% SLM, using the parameters of such SLM.
N = length(X);
pattern = [1 -1; -1 1];
x_shift = 100;
y_shift = 0;
checkers = repmat(pattern, N/2, N/2); 
Qlens = exp(-1i*(pi/(lambda*f))*(X.^2+Y.^2)).^checkers;   % propagation transfer function of lens with focus f
Qshift = exp(1i*(pi/lambda)*(X*x_shift + Y*y_shift)).^checkers; % transfer function of linear phase shift
if SLM_type == 1
    propagated = A;
    graphs = false;     % disable graphs printing while using standard TIE
else
    propagated = A.*Qlens.*Qshift;
end

saveFigure(propagated.*conj(propagated),'Intensity After SLM Transfer ',nameOfPlane,graphs);
saveFigure(angle(Qlens.*Qshift),'Mask Applied to Fourier Plane','',graphs);
end
