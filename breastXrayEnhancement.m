% EE 5353 - Biomedical Imaging, UTSA
% john jenkinson 2014
% project 1... X-ray Image Enhancement
% Galaxy Image Enhancment (testing)

close all; clear all; clc;
%% A. Read the image from file
% Breast Xray Image
Y=imread('breast_Xray.tif');
[N M]=size(Y);
Y_double=double(Y);
Y_neg=255-Y; 
%Y_neg=imadjust(Y,[0 1],[1 0],1); negative image by function
%Y_neg=imcomplement(Y); negative image matlab function
Y_neg_double=double(Y_neg);

% Spiral Galaxy Image
G=imread('AC-2387spiral.jpg');
G=G(:,:,1);
[V D]=size(G);
G_double=double(G);
G_neg=255-G;
G_neg_double=double(G_neg);

%% B. Process image 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Grayscale transformations (GST)

% Logarithmic GST 
% Xray
c=255/log(256);
r=0:255;
T1=c*log(1+r);
Y_log=zeros(N,M);
Y_neg_log=zeros(N,M);
for n=1:N
    for m=1:M
        r=Y_double(n,m)+1;
        s=T1(r);
        Y_log(n,m)=s;
        %
        u=Y_neg_double(n,m)+1;
        i=T1(u);
        Y_neg_log(n,m)=i;
    end
end
% Galaxy
G_log=zeros(V,D);
G_neg_log=zeros(V,D);
for v=1:V
    for d=1:D
        o=G_double(v,d)+1;
        p=T1(o);
        G_log(v,d)=p;
        %
        x=G_neg_double(v,d)+1;
        a=T1(x);
        G_neg_log(v,d)=a;
    end
end

% Square-root GST
% Xray
c=255/sqrt(256);
r=0:255;
T2=c*sqrt(1+r);
Y_sqrt=zeros(N,M);
Y_neg_sqrt=zeros(N,M);
for n=1:N
    for m=1:M
        r=Y(n,m)+1;
        s=T2(r);
        Y_sqrt(n,m)=s;
        %
        u=Y_neg(n,m)+1;
        i=T2(u);
        Y_neg_sqrt(n,m)=i;
    end
end
% Galaxy
G_sqrt=zeros(V,D);
G_neg_sqrt=zeros(V,D);
for v=1:V
    for d=1:D
        o=G_double(v,d)+1;
        p=T1(o);
        G_sqrt(v,d)=p;
        %
        x=G_neg_double(v,d)+1;
        a=T1(x);
        G_neg_sqrt(v,d)=a;
    end
end

% Power transform
% Xray
r=0:255;
g=1;
T3=(1+r).^g;
mm=max(T3);
c=255/(mm*256);
T3=c*(1+r).^g;
Y_pwr=zeros(N,M);
Y_neg_pwr=zeros(N,M);
for n=1:N
    for m=1:M
        r=Y(n,m)+1;
        s=T3(r);
        Y_pwr(n,m)=s;
        %
        u=Y_neg(n,m)+1;
        i=T3(u);
        Y_neg_pwr(n,m)=i;
    end
end
% Galaxy
G_pwr=zeros(V,D);
G_neg_pwr=zeros(V,D);
for v=1:V
    for d=1:D
        o=G_double(v,d)+1;
        p=T1(o);
        G_pwr(v,d)=p;
        %
        x=G_neg_double(v,d)+1;
        a=T1(x);
        G_neg_pwr(v,d)=a;
    end
end

% Gray-level slicing GST
% this transform is used if the user wishes to highlight a certain range
% of intensities of the image.  Brightness of the image outside this range
% will be preserved.
% Xray
Y_slice=zeros(N,M);
Y_neg_slice=zeros(N,M);
for n=1:N
    for m=1:M
        if (Y(n,m)>=48) && (Y(n,m)<=96)
        Y_slice(n,m)=120;
        else
            Y_slice(n,m)=Y(n,m);
        end
        if (Y_neg(n,m)>=120) && (Y_neg(n,m)<=180)
        Y_neg_slice(n,m)=79;
        else
            Y_neg_slice(n,m)=Y_neg(n,m);
        end
    end
end
% Galaxy
G_slice=zeros(V,D);
G_neg_slice=zeros(V,D);
for v=1:V
    for d=1:D
        if (G(v,d)>=48) && (G(v,d)<=96)
        G_slice(v,d)=120;
        else
            G_slice(v,d)=G(v,d);
        end
        if (G_neg(v,d)>=120) && (G_neg(v,d)<=180)
        G_neg_slice(v,d)=79;
        else
            G_neg_slice(v,d)=G_neg(v,d);
        end
    end
end

% Constrast-stretching
r=0:255;
Y_stretch=zeros(N,M);
Y_neg_stretch=zeros(N,M);
E=2;
i=80;
T=1./(1+(i./r).^E);
for n=1:N
    for m=1:M
        r=Y(n,m)+1;
        s=T(r);
        Y_stretch(n,m)=s;
        %Y_stretch(n,m)=1./(1+(i./Y_double(n,m)).^E);
        %
        u=Y_neg(n,m)+1;
        i=T(u);
        Y_neg_stretch(n,m)=i;
    end
end
% Galaxy
r=0:255;
G_stretch=zeros(V,D);
G_neg_stretch=zeros(V,D);
E=2;
i=80;
T=1./(1+(i./r).^E);
for v=1:V
    for d=1:D
        o=G(v,d)+1;
        p=T(o);
        G_stretch(v,d)=p;
        %G_stretch(v,d)=1./(1+(i./G_double(n,m)).^E);
        %
        x=G_neg(v,d)+1;
        a=T(x);
        G_neg_stretch(v,d)=a;
    end
end

% Histogram equalization
% Xray
Y_Heq=histeq(Y);
Y_neg_Heq=histeq(Y_neg);
% Galaxy
G_Heq=histeq(G);
G_neg_Heq=histeq(G_neg);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. Propose best method for image enhancement

% Method of log-transform then power-transform (power g=1:0.5:6)
for g=1:0.5:5
c=255/log(256);
Y_double=double(Y);
r=0:255;
CT1=c*log(1+r);
%g=1;
CT2=(1+r).^g;
mm=max(CT2);
c=255/(mm*256);
CT2=c*(1+r).^g;
Y_log_pwr=zeros(N,M);
for n=1:N
    for m=1:M
        r=Y_double(n,m);
        s=CT1(r);
        s=round(s);
        p=CT2(s);
        Y_log_pwr(n,m)=p;
    end
end
imagesc(Y_log_pwr)
title('Method of log-transform then power-transform (power g=1:0.5:6)');
axis image
colormap(gray)
pause(1)
clear CT1;clear CT2;
end

% Method of log-transform then square-root transfrom
c=255/log(256);
Y_double=double(Y);
r=0:255;
CT3=c*log(1+r);
c=255/sqrt(256);
CT4=c*sqrt(1+r);
Y_log_sqrt=zeros(N,M);
for n=1:N
    for m=1:M
        r=Y_double(n,m);
        s=CT3(r);
        s=round(s);
        p=CT4(s);
        Y_log_sqrt(n,m)=p;
    end
end

% Method of square-root transfrom then power-transfrom
r=0:255;
c=255/sqrt(256);
CT5=c*sqrt(1+r);
g=1;
CT6=(1+r).^g;
mm=max(CT6);
c=255/mm;
CT6=c*(1+r).^g;
Y_sqrt_pwr=zeros(N,M);
for n=1:N
    for m=1:M
        r=Y(n,m);
        s=CT5(r);
        s=round(s);
        p=CT6(s);
        Y_sqrt_pwr(n,m)=p;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. Direct method of enhancement
% (1)fn,m --> [GST] --> \hat{f}n,m
% (2)fn,m --> gn,m=255-fn,m --> [GST] --> \hat{g}n,m --> \hat{f}n,m =
%   = (255-\hat{g}n,m) 
% Both(1) and (2) are performed for each GST in section B.1.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4. EME Calculations

% breast X-ray
EME_Y=eme2(Y_double,N,M,7)
EME_Y_log=eme2(Y_log,N,M,7)
EME_Y_sqrt=eme2(Y_sqrt,N,M,7)
EME_Y_pwr=eme2(Y_pwr,N,M,7)
EME_Y_slice=eme2(Y_slice,N,M,7)
EME_Y_Heq=eme2(Y_Heq,N,M,7)

% Spiral Galaxy
EME_G=eme2(G_double,V,D,7)
EME_G_log=eme2(G_log,V,D,7)
EME_G_sqrt=eme2(G_sqrt,V,D,7)
EME_G_pwr=eme2(G_pwr,V,D,7)
EME_G_slice=eme2(G_slice,V,D,7)
EME_G_Heq=eme2(G_Heq,V,D,7)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5. Display results

% Direct Method Enhancement X-ray
h_fig1=figure;
set(h_fig1,'Name','Direct Method GSTs');
colormap(gray);
subplot(3,2,1);
imagesc(Y); title('breast Xray');
axis image; axis off;
subplot(3,2,2);
imagesc(Y_log); title('logarithmic GST');
axis image; axis off;
subplot(3,2,3);
imagesc(Y_sqrt); title('square root GST');
axis image; axis off;
subplot(3,2,4);
imagesc(Y_pwr); title('power GST g=1');
axis image; axis off;
subplot(3,2,5);
imagesc(Y_slice); title('grey-level slicing, [48:96] --> 120');
axis image; axis off;
subplot(3,2,6);
imagesc(Y_stretch); title('constrast stretching');
axis image; axis off;

% Histogram Equalization Xray
h_fig2=figure;
set(h_fig2,'Name','Histogram Equalization');
colormap(gray)
subplot(2,2,1)
imagesc(Y); axis image
title('Original Image')
subplot(2,2,2)
imagesc(Y_Heq); axis image 
title('Histogram equalization');
subplot(2,2,3)
imagesc(Y_neg); axis image
title('Original Negative Image')
subplot(2,2,4)
imagesc(Y_neg_Heq); axis image 
title('Negative Histogram equalization');


% Best GST Method Enhancement X-ray
h_fig3=figure;
set(h_fig3,'Name','Best Method GST');
colormap(gray);
subplot(1,3,1);
imagesc(Y_log_pwr); title('log --> power');
axis image; axis off;
subplot(1,3,2);
imagesc(Y_log_sqrt); title('log --> sqrt');
axis image; axis off;
subplot(1,3,3);
imagesc(Y_sqrt_pwr); title('sqrt --> power');
axis image; axis off;

% Negative Image Enhancment X-ray
h_fig4=figure;
set(h_fig4,'Name','Negative Image GSTs');
colormap(gray);
subplot(3,2,1);
imagesc(Y_neg); title('negative breast Xray');
axis image; axis off;
subplot(3,2,2);
imagesc(Y_neg_log); title('logarithmic GST');
axis image; axis off;
subplot(3,2,3);
imagesc(Y_neg_sqrt); title('square root GST');
axis image; axis off;
subplot(3,2,4);
imagesc(Y_neg_pwr); title('power GST g=1');
axis image; axis off;
subplot(3,2,5);
imagesc(Y_neg_slice); title('grey-level slicing, [48:96] --> 120');
axis image; axis off;
subplot(3,2,6);
subplot(3,2,6);
imagesc(Y_neg_stretch); title('constrast stretching');
axis image; axis off;


% Direct Method Enhancement Galaxy
h_fig5=figure;
set(h_fig5,'Name','Direct Method GSTs');
colormap(gray);
subplot(3,2,1);
imagesc(G); title('Spiral Galaxy');
axis image; axis off;
subplot(3,2,2);
imagesc(G_log); title('logarithmic GST');
axis image; axis off;
subplot(3,2,3);
imagesc(G_sqrt); title('square root GST');
axis image; axis off;
subplot(3,2,4);
imagesc(G_pwr); title('power GST g=1');
axis image; axis off;
subplot(3,2,5);
imagesc(G_slice); title('grey-level slicing, [48:96] --> 120');
axis image; axis off;
subplot(3,2,6);
imagesc(G_Heq); title('histogram equalization');
axis image; axis off;

% Histogram Equalization Xray
h_fig6=figure;
set(h_fig6,'Name','Histogram Equalization');
colormap(gray)
subplot(2,2,1)
imagesc(G); axis image
title('Original Image')
subplot(2,2,2)
imagesc(G_Heq); axis image 
title('Histogram equalization');
subplot(2,2,3)
imagesc(G_neg); axis image
title('Original Negative Image')
subplot(2,2,4)
imagesc(G_neg_Heq); axis image 
title('Negative Histogram equalization');

% Negative Image Enhancment Galaxy
h_fig7=figure;
set(h_fig7,'Name','Negative Image GSTs');
colormap(gray);
subplot(3,2,1);
imagesc(G_neg); title('negative Galaxy');
axis image; axis off;
subplot(3,2,2);
imagesc(G_neg_log); title('logarithmic GST');
axis image; axis off;
subplot(3,2,3);
imagesc(G_neg_sqrt); title('square root GST');
axis image; axis off;
subplot(3,2,4);
imagesc(G_neg_pwr); title('power GST g=1');
axis image; axis off;
subplot(3,2,5);
imagesc(G_neg_slice); title('grey-level slicing, [48:96] --> 120');
axis image; axis off;
subplot(3,2,6);
imagesc(G_neg_stretch); title('contrast stretching');
axis image; axis off;









