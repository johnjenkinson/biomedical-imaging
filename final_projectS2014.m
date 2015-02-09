% EE 5353 - Biomedical Imaging
% Final Project - Image Reconstruction by Projections
% john jenkinson, UTSA Spring 2014
%
clear all; close all; clc;
%method 1.
tic
X=linspace(-3,3,256);
Y=linspace(-3,3,256);
phi=pi/9;
x0=-0.2;
y0=0.3;
a=2; a=a^2;
b=0.5; b=b^2;
c1=cos(phi);
c2=sin(phi);
u=zeros(256);
for y=1:length(Y)
    for x=1:length(X)    
                if(((X(x)*c1+Y(y)*c2-x0).^2/a+...
                (Y(y)*c1-X(x)*c2-y0).^2/b)<1)
        u(x,y)=20;
                else if(((X(x)*c1+Y(y)*c2-x0).^2/a+...
                (Y(y)*c1-X(x)*c2-y0).^2/b)==1)
        u(x,y)=100;
            end
        end
    end
end
toc
            %imshow(u,[])
%method 2.
tic
    [XX,YY]=meshgrid(X,Y);
    v = zeros(256,256);
    v( ((XX.*cos(phi) + YY.*sin(phi)-x0).^2./a+ ...
        (YY.*cos(phi)-XX.*sin(phi)-y0).^2./b)<=1) = 20;
toc
            % imshow(v)
tic
theta=0:179;
R=radon(u,theta);
I=iradon(R,theta);
            % imshow(I,[])
toc
E=norm(I(1:256,1:256)-u);
fprintf('\n error between sampled and reconstructed images is %6.4f \n',E)
hfig=figure;
set(hfig,'Menubar','None','Name','Image Reconstruction by Projections')
subplot(2,2,1); 
imshow(u,[])
stitle=title('Sampled image 256x256');
set(stitle,'FontName','Times','FontSize',12);
subplot(2,2,2)
imshow(I,[])
stitle=title('Reconstructed image 258x258');
set(stitle,'FontName','Times','FontSize',12);
subplot(2,2,[3,4]);
imshow(I(1:256,1:256)-u);
stitle=title('Reconstruction error I-u');
set(stitle,'FontName','Times','FontSize',12);


