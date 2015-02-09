%% EE 5353: Biomedical Imaging, UTSA, April 2014
%john jenkinson
clear all; close all; clc;
f = imread('capLAN_BL1.TIF');
% 1. Compose the gray-scale image.
f = (f(:,:,1) + f(:,:,2) + f(:,:,3)) / 3;
% 2,3. Calculate dilation and erosion
% opening and closing.
B1 = strel('disk', 3, 8);
fd1 = imdilate(f, B1);
fe1 = imerode(f, B1);
fo1 = imopen(f, B1);
fo1x = imdilate(fe1, B1);
fc1 = imclose(f, B1);
fc1x = imerode(fd1, B1);
if (fo1 == fo1x)>0;
    disp('matlab imopen equals matlab erosion followed by matlab dilation');
else
    disp('matlab imopen does not equal matlab erosion followed by matlab dilation');
end
if (fc1 == fc1x)>0;
    disp('matlab imclose equals matlab dilation followed by matlab erosion');
else
    disp('matlab imclose does not equal matlab dilation followed by matlab erosion');
end
cross = [0 1 0; 1 1 1; 0 1 0];
B2 = strel(cross);
fd2 = imdilate(f, B2);
fe2 = imerode(f, B2);
fo2 = imopen(f, B2);
fc2 = imclose(f, B2);
xcross = [1 0 1; 0 1 0; 1 0 1];
B3 = strel(xcross);
fd3 = imdilate(f, B3);
fe3 = imerode(f, B3);
fo3 = imopen(f, B3);
fc3 = imclose(f, B3);
dil_fig=figure;
set(dil_fig,'Name','Dilation','Menubar','None')
subplot(2,2,1)
imshow(f,[])
stitle=title('Original image');
set(stitle,'FontName','Times','Fontsize',12);
subplot(2,2,2)
imshow(fd1,[])
stitle=title('Dilation by disk 3x3');
set(stitle,'FontName','Times','Fontsize',12);
subplot(2,2,3)
imshow(fd2,[])
stitle=title('Dilation by cross 3x3');
set(stitle,'FontName','Times','Fontsize',12);
subplot(2,2,4)
imshow(fd3,[])
stitle=title('Dilation by xcross 3x3');
set(stitle,'FontName','Times','Fontsize',12);
%
ero_fig=figure;
set(ero_fig,'Name','Erosion','Menubar','None')
subplot(2,2,1)
imshow(f,[])
stitle=title('Original image');
set(stitle,'FontName','Times','Fontsize',12);
subplot(2,2,2)
imshow(fe1,[])
stitle=title('Erosion by disk 3x3');
set(stitle,'FontName','Times','Fontsize',12);
subplot(2,2,3)
imshow(fe2,[])
stitle=title('Erosion by cross 3x3');
set(stitle,'FontName','Times','Fontsize',12);
subplot(2,2,4)
imshow(fe3,[])
stitle=title('Erosion by xcross 3x3');
set(stitle,'FontName','Times','Fontsize',12);
%
open_fig=figure;
set(open_fig,'Name','Opening','Menubar','None')
subplot(2,2,1)
imshow(f,[])
stitle=title('Original image');
set(stitle,'FontName','Times','Fontsize',12);
subplot(2,2,2)
imshow(fo1,[])
stitle=title('Opening by disk 3x3');
set(stitle,'FontName','Times','Fontsize',12);
subplot(2,2,3)
imshow(fo2,[])
stitle=title('Opening by cross 3x3');
set(stitle,'FontName','Times','Fontsize',12);
subplot(2,2,4)
imshow(fo3,[])
stitle=title('Opening by xcross 3x3');
set(stitle,'FontName','Times','Fontsize',12);
%
open_fig=figure;
set(open_fig,'Name','Closing','Menubar','None')
subplot(2,2,1)
imshow(f,[])
stitle=title('Original image');
set(stitle,'FontName','Times','Fontsize',12);
subplot(2,2,2)
imshow(fc1,[])
stitle=title('Closing by disk 3x3');
set(stitle,'FontName','Times','Fontsize',12);
subplot(2,2,3)
imshow(fc2,[])
stitle=title('Closing by cross 3x3');
set(stitle,'FontName','Times','Fontsize',12);
subplot(2,2,4)
imshow(fc3,[])
stitle=title('Closing by xcross 3x3');
set(stitle,'FontName','Times','Fontsize',12);
% 4. Caculate gradients by morphological operators.
G1=f-fe1;%f-(f erode B) !!best gradient operation for the image
G2=fd1-f;%(f dilate B)-f
G3=fd1-fe1;%(f dilate B)-(f erode B)
G4=fc1-fo1;%(f open B)-(f close B)
grad_fig=figure;
set(grad_fig,'Name','Gradients by Opening & Closing','Menubar','None')
subplot(2,2,1)
imshow(4*G1,[])
stitle=title('f-(f erode B) !!best gradient operation for the image');
set(stitle,'FontName','Times','Fontsize',12);
subplot(2,2,2)
imshow(4*G2,[])
stitle=title('(f dilate B)-f');
set(stitle,'FontName','Times','Fontsize',12);
subplot(2,2,3)
imshow(4*G3,[])
stitle=title('(f dilate B)-(f erode B)');
set(stitle,'FontName','Times','Fontsize',12);
subplot(2,2,4)
imshow(4*G4,[])
stitle=title('(f open B)-(f close B)');
set(stitle,'FontName','Times','Fontsize',12);
% 5. Calculate global threshold and threshold by Otsu's method.
count=0;
T=mean2(f);
done=false;
while ~done
    count=count+1;
    g=f>T;
    Tnext=0.5*(mean(f(g))+mean(f(~g)));
    done=abs(T-Tnext)<0.5;
    T=Tnext;
end
T
g=im2bw(f,T/255);
% Otsu's method.
[level em]=graythresh(f);
255*level
Terror=norm(T-255*level)
gt_fig=figure;
set(gt_fig,'Name','Global Thresholding','Menubar','None')
subplot(2,3,1)
imshow(f,[])
stitle=title('Original image, f');
set(stitle,'FontName','Times','Fontsize',12);
subplot(2,3,2)
imshow(g,[])
stitle=title('Threshold image, g');
set(stitle,'FontName','Times','Fontsize',12);
subplot(2,3,3)
imshow(g,[])
stitle=title('Otsu threshold');
set(stitle,'FontName','Times','Fontsize',12);
subplot(2,3,[4,5,6])
imhist(f)
stitle=title('Threshold image histogram');
set(stitle,'FontName','Times','Fontsize',12);
% 6. Segmentation by watershed.
gc=~g;
D=bwdist(gc);
L=watershed(-D);
w=L==0;
g2=g&~w;
water_fig=figure;
set(water_fig,'Name','Watershed segmentation','Menubar','None')
subplot(3,2,1)
imshow(g,[])
stitle=title('Binary image');
set(stitle,'FontName','Times','Fontsize',12);
subplot(3,2,2)
imshow(gc,[])
stitle=title('Binary complement');
set(stitle,'FontName','Times','Fontsize',12);
subplot(3,2,3)
imshow(D,[])
stitle=title('Distance transform');
set(stitle,'FontName','Times','Fontsize',12);
subplot(3,2,4)
imshow(L,[])
stitle=title('Watershed ridge lines');
set(stitle,'FontName','Times','Fontsize',12);
subplot(3,2,[5,6])
imshow(g2,[])
stitle=title('Segmented image');
set(stitle,'FontName','Times','Fontsize',12);
%
disp('Gradient operator: f-(f erode B) !!best method for segmentation');



