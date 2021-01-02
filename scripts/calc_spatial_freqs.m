function [k_x,k_y] = calc_spatial_freqs(SLM_type, X, Y)
if SLM_type == 1
    % take a meshgrid of the plane and calculate the spatial frequencies of
    % that plane

    X(X==0) = sqrt(realmin);            %replace zeros with sqrt(realmin)
    Y(Y==0) = sqrt(realmin);            %replace zeros with sqrt(realmin)
    %replacing the zeros with the sqrt of the minimum values so we can
    %divide by the elements with zero value. the swaure root is to ensure
    %that we can later on sqaure those numbers without geting Inf value
    k_x = 2*pi./X;          % grid of x axis spatial frequencies
    k_y = 2*pi./Y;          % grid of y axis spatial frequencies
  
end
end