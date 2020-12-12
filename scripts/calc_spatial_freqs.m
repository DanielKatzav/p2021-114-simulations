function [k_x,k_y] = calc_spatial_freqs(lambda,SLM_type)
k_0 = 2*pi/lambda;
if SLM_type == 1
    k_x = k_0/sqrt(2);
    k_y = k_x;
end
end