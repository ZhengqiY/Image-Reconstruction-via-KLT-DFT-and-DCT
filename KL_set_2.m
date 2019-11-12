clear
clc
close all

% Load the image set_2

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
%compute covariance matrix
cov_mat = x'*x;
[V,D] = eig(cov_mat);         %eigen values of cov matrix
V = x*V*(abs(D))^-0.5;               
 
%image decomposition coefficients
KLCoef =  x'*V;

%top 100 coefficients
list100 = KLCoef(:);
list100 = list100';
sortl = sort(list100);
gate = sortl(301);
KLCoef100 = KLCoef;

for i = 1:20
    for j = 1:20
        if KLCoef100(i,j) < gate
            KLCoef100(i,j) = 0;
        end
    end
end

%top 20 coefficients
list20 = KLCoef(:);
list20 = list20';
sortl = sort(list20);
gate = sortl(381);
KLCoef20 = KLCoef;

for i = 1:20
    for j = 1:20
        if KLCoef20(i,j) < gate
            KLCoef20(i,j) = 0;
        end
    end
end

%top 10 coefficients
list10 = KLCoef(:);
list10 = list10';
sortl = sort(list10);
gate = sortl(391);
KLCoef10 = KLCoef;

for i = 1:20
    for j = 1:20
        if KLCoef10(i,j) < gate
            KLCoef10(i,j) = 0;
        end
    end
end

%top 5 coefficients
list5 = KLCoef(:);
list5 = list5';
sortl = sort(list5);
gate = sortl(396);
KLCoef5 = KLCoef;

for i = 1:20
    for j = 1:20
        if KLCoef5(i,j) < gate
            KLCoef5(i,j) = 0;
        end
    end
end

%reconstruction of Image
i = 1:20;  %index of image to be recontructed
reconst = V*KLCoef';

% % % for image_index = i
% % %     
% % %     diff = abs(reconst(:,image_index) - x(:,image_index));
% % %     figure;
% % %     subplot(1,2,1); imshow((reshape(avrgx+reconst(:,image_index), imsize))); title('Reconstructed');
% % %     subplot(1,2,2); imshow((reshape(avrgx+x(:,image_index), imsize)));title('original');
% % %     saveas(gcf, sprintf('reconst set 2 %d.jpg', image_index));
% % %     close all
% % % end

% RMSE vs Top 100 KLCoef, Img 20
reconst100 = V*KLCoef100';
image_index = 20;
diff = abs(reconst100(:,image_index) - x(:,image_index));
rmse100 = sqrt(sum(diff.^2)/k);

% RMSE vs Top 20 KLCoef, Img 20
reconst20 = V*KLCoef20';
image_index = 20;
diff = abs(reconst20(:,image_index) - x(:,image_index));
rmse20 = sqrt(sum(diff.^2)/k);

% RMSE vs Top 10 KLCoef, Img 20
reconst10 = V*KLCoef10';
image_index = 20;
diff = abs(reconst10(:,image_index) - x(:,image_index));
rmse10 = sqrt(sum(diff.^2)/k);

% RMSE vs Top 5 KLCoef, Img 20
reconst5 = V*KLCoef5';
image_index = 20;
diff = abs(reconst5(:,image_index) - x(:,image_index));
rmse5 = sqrt(sum(diff.^2)/k);

figure
subplot(2,2,1); imshow((reshape(avrgx+reconst100(:,image_index), imsize))); title('Reconstructed Top 100 Co');
subplot(2,2,2); imshow((reshape(avrgx+reconst20(:,image_index), imsize))); title('Reconstructed Top 20 Co');
subplot(2,2,3); imshow((reshape(avrgx+reconst10(:,image_index), imsize))); title('Reconstructed Top 10 Co');
subplot(2,2,4); imshow((reshape(avrgx+reconst5(:,image_index), imsize))); title('Reconstructed Top 5 Co');

%Plot RMSE vs KLCoef Num
rmse_set_1 = [rmse100, rmse20, rmse10, rmse5];
KLCoefnum = [100, 20, 10, 5];
figure
plot(KLCoefnum, rmse_set_1);
title('Set 2 RMSE vs KLCoef Num');