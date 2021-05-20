function [centered] = centerImage(full_img,partial_img)
% center a small image inside a larger image. 
h_move = fix((size(full_img,1) - partial_img(1) + 1)/2);  % pixels to shift height for centering
w_move = fix((size(full_img,2) - partial_img(2) + 1)/2);  % pixels to shift width for centering
centered = circshift(full_img,[h_move w_move]);     % center image elements
end
