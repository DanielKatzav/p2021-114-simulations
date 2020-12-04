function  y=ift2(mt)
y = fftshift(ifft2(mt));
[m,n] = size(mt);
M = mod(m,2);
N = mod(n,2);
Xind = ceil(-m/2) : 1 : ceil(m/2)-1;
Yind = ceil(-n/2) : 1 : ceil(n/2)-1;
if (M==0) && (N==0)
    correctionX = exp(-1i*pi*(Xind));
    correctionY = exp(-1i*pi*(Yind));
    correcrionMatrix = correctionX.'*correctionY;
elseif (M==1) && (N==1)
    correctionX = exp(-1i*pi*(m-1)*(Xind/m));
    correctionY = exp(-1i*pi*(n-1)*(Yind/n));
    correcrionMatrix = correctionX.'*correctionY;
else
    warning('Your matrix size is not compatible')
end
y = y.*correcrionMatrix;
