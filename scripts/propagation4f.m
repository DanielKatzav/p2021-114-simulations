function [propagated] = propagation4f(A,distances,focus,lambda,X,Y,graphs, nameOfPlane)
% propagation4f will simulate the propagation of an EM wave through a 4f
% system. Input must include a 1x2 list of travel distances where the first
% element is the distance from object to first lens, and second element is
% distance from second lens to image [z_o z_i], and must include a 1x2 list
% of lens focus lengths [f1 f2]. lambda is the wavelength and X,Y are the
% meshgrid values for the 2D wavefront.
if length(focus) ~= 2 || length(distances) ~= 2
    disp('Must provide a 2 element array for distances and foci')
    propagated = -1;
    return
end

A_before_lens1 = distancePropagate(A,distances(1), lambda, X,Y,graphs, strcat(nameOfPlane, "before first lens"));          % image after propagation over z_o
A_after_lens1 = lensPropagate(A_before_lens1,focus(1),lambda,X,Y,graphs, strcat(nameOfPlane, "after first lens"));        % image after propagation trough lens 1
A_before_fourier = distancePropagate(A_after_lens1,focus(1),lambda,X,Y,graphs,strcat(nameOfPlane, "at fourier plane"));
A_before_lens2 = distancePropagate(A_before_fourier,focus(2),lambda,X,Y,graphs, strcat(nameOfPlane, "before second lens"));  % image after propagating between lenses
A_after_lens2 = lensPropagate(A_before_lens2,focus(2),lambda,X,Y,graphs, strcat(nameOfPlane, "after second lens"));        % image after propagating through lens 2
propagated = distancePropagate(A_after_lens2,distances(2), lambda, X,Y,graphs, strcat(nameOfPlane, "propagated"));  % image after propagationg over z_i

end