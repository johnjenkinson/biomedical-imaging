clear;
%Low-pass Gaussian Filter (LPGF);
f = imread('monarchg.png');
[N M] = size(f);
f = double(f);
% Compute padded size to prepare for FFT-based filtering.
NM = [N M];
PQ = 2*NM;
% Set up range of variables.
u = single(0:(PQ(1)-1));
v = single(0:(PQ(2)-1));
% Compute the indices for use in meshgrid.
idx = find(u > PQ(1)/2);
u(idx) = u(idx) - PQ(1);
idy = find(v > PQ(2)/2);
v(idy) = v(idy) - PQ(2);
% Compute the meshgrid arrays.
[V, U] = meshgrid(v, u);
% C = hypot(A,B) returns SQRT(ABS(A).^2+ABS(B).^2)
D = hypot(V, U);
% Obtain the FFT of the padded input.
F = fft2(f, PQ(1), PQ(2));
% Gaussian filter.
%D = D.^2;
D0 = 0.05*PQ(2);
% D0 = D0^2;
H = exp(-(D.^2)/(2*D0^2));
% Perfrom filtering.
g = ifft2(H.*F);
% Crop to original size.
g = g(1:size(f, 1), 1:size(f, 2)); % g is of class single here.
% Convert output to original class type.
g = uint8(g);
% Plot results.
figure
subplot(2, 2, 1)
imshow(f, [])
title('Original image, f');
subplot(2, 2, 2)
imshow(abs(fftshift(F)), [0 10000])
title('F = DFT(f)');
subplot(2, 2, 3)
imshow(fftshift(H))
title('Gaussian frequency filter function, H');
subplot(2, 2, 4)
imshow(g, [])
title('Filtered image,g');
