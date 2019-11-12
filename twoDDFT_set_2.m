clear
clc
close all
% Load the image set_1
set_2 = dir('set_2');
len = length(set_2);
k = 0;
for i = 3:len
    filename = set_2(i).name;
    img = imread(filename);
    k = k + 1;
    x(:,k) = img(:);
end

nImages = k;                     %total number of images
imsize = size(img);       %size of image (they all should have the same size) 
nPixels = imsize(1)*imsize(2);   %number of pixels in image
x = double(x)/255;               %convert to double and normalize
%Calculate the average
avrgx = mean(x')';
for i=1:nImages
    x(:,i) = x(:,i) - avrgx; % substruct the average
end

image_index = 20; % Img

image_2d = x(:,image_index);
image_2d = reshape(image_2d+avrgx, imsize);

A=fft2(image_2d);
figure(1);
imshow(abs(ifft2(A)));title('2D-DFT Reconst');

%top 10000 coefficients
list100 = A(:);
list100 = list100';
sortl = sort(list100);
gate = sortl(120*160-9999);
A100 = A;

for i = 1:120
    for j = 1:160
        if A100(i,j) < gate
            A100(i,j) = 0;
        end
    end
end

figure(2);
subplot(2,2,1); imshow(abs(ifft2(A100))); title('Reconstructed Top 10000 Co');

%top 7000 coefficients
list20 = A(:);
list20 = list20';
sortl = sort(list20);
gate = sortl(120*160-6999);
A20 = A;

for i = 1:120
    for j = 1:160
        if A20(i,j) < gate
            A20(i,j) = 0;
        end
    end
end
subplot(2,2,2); imshow(abs(ifft2(A20))); title('Reconstructed Top 7000 Co');

%top 4000 coefficients
list10 = A(:);
list10 = list10';
sortl = sort(list10);
gate = sortl(120*160-3999);
A10 = A;

for i = 1:120
    for j = 1:160
        if A10(i,j) < gate
            A10(i,j) = 0;
        end
    end
end
subplot(2,2,3); imshow(abs(ifft2(A10))); title('Reconstructed Top 4000 Co');

%top 1000 coefficients
list5 = A(:);
list5 = list5';
sortl = sort(list5);
gate = sortl(120*160-999);
A5 = A;

for i = 1:120
    for j = 1:160
        if A5(i,j) < gate
            A5(i,j) = 0;
        end
    end
end
subplot(2,2,4); imshow(abs(ifft2(A5))); title('Reconstructed Top 1000 Co');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RMSE vs Top 10000 DFTCoef, Img 20
diff = abs(double(abs(ifft2(A100)) - image_2d));
rmse100 = sqrt(sum(diff(:).^2)/k);

% RMSE vs Top 7000 DFTCoef, Img 20
diff = abs(double(abs(ifft2(A20)) - image_2d));
rmse20 = sqrt(sum(diff(:).^2)/k);

% RMSE vs Top 4000 DFTCoef, Img 20
diff = abs(double(abs(ifft2(A10)) - image_2d));
rmse10 = sqrt(sum(diff(:).^2)/k);

% RMSE vs Top 1000 DFTCoef, Img 20
diff = abs(double(abs(ifft2(A5)) - image_2d));
rmse5 = sqrt(sum(diff(:).^2)/k);

%Plot RMSE vs DFTCoef Num
rmse_set_2 = [rmse100, rmse20, rmse10, rmse5];
DFTCoefnum = [10000, 7000, 4000, 1000];
figure(3);
plot(DFTCoefnum, rmse_set_2);
title('Set 2 RMSE vs DFTCoef Num');

