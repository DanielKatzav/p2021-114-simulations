function [k_x,k_y] = calc_spatial_freqs(k_0, SLM_type, image_intensity)
if SLM_type == 1
    
    [width,height] = size(image_intensity);
    center_x = round(width/2);
    center_y =  round(height/2);
    
    k_x = ones(width);
    k_y = ones(height);
    k_x(center_x, :) = realmax;
    k_y(:, center_y) = realmax;
    
    for x = 1:center_x-1 
        k_x(center_x+x,:) = 2*pi/x;
        k_x(center_x-x,:) = -2*pi/x;
    end
    
    for y = 1:center_y-1 
        k_y(:, center_y+y) = 2*pi/y;
        k_y(:, center_y-y) = -2*pi/y;
    end   
    
end
end