% for color monarch image
% f = imread('monarch.png');
% f = f(:,:,1); to convert to grayscale
clear; close; clc;
tic
f = imread('monarchg.png');
f = double(f);
[N M] = size(f);

%% attempt 2
figure;
subplot(2, 2, 1)
imshow(f, [])

wx = linspace(-pi, pi, N);
wy = linspace(-pi, pi, M);
[U V] = meshgrid(wy, wx);
D = hypot(U, V);

H = zeros(N, M);
D0 = 10;
D0 = D0^2;
for i = 1:N
    for j = 1:M
        arg = -D(i,j)/(2*D0);
        H(i, j) = exp(arg);
    end
end

H = H / sum(H(:));
% H = fftshift(fft2(H));
F = fftshift(fft2(f));

subplot(2, 2, 2)
imshow(abs(F)/N/20)

subplot(2, 2, 3)
imshow(H*500)

% Circular convolution.
G = conv2(F,H);

% Inverse 2-D DFT.
g = real(ifft2(G));
g = uint8(g);

subplot(2, 2, 4)
imshow(g)
toc
%% attempt 1
% 
% % Define the Gaussian Low-pass filter.
% H = zeros(N, M);
% %cu=0; cv=0;
% for u = 1:N
%     %u1=(u-N/2)*(2*pi/N);
%     for v = 1:M
%         %v1=(v-M/2)*(2*pi/M);
% %         D = (u1 - cu).^2 + (v1 - cv).^2;
%         %H(u,v) = exp(-D/(2*D0)); % n = 1.
%         %H(u,v) = sqrt(1/(1 + (D0/D)^6)); 
%     end
% end
% 
% subplot(2, 2, 3)
% imshow(H)
% 
% % Circular convolution
% G = F.*H; % image filtration
% g1 = real(ifft2(G));
% g=uint8(g1);
% 
% subplot(2, 2, 4)
% imshow(g1*10)
