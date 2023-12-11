% This script calbrates the two channels (561 and 647) used for imaging DNA
% loci in sequential FISH for the strain-specific marking specifially to
% adjust for Z-drift between aquisitions of channels
% ------------------------------------------------------------------------------------------------------------------------------------------
% Input:  multidimensional matrices (3D image stacks) of nuclear signal of
% the primary imaging round and the strain imaging round
%  
% 
% ------------------------------------------------------------------------------------------------------------------------------------------
% Output: DeltaZ560 and DeltaZ647 to adjust for Zdrift between primary imaging round and strain imaging 
% ------------------------------------------------------------------------------------------------------------------------------------------
% Silvia Gutnik
% 10.03.2021
% Version 1.0
% 
% using fitFoci and selectROI functions from the Zhuang lab
%%
clear all
close all

SampleNum = 18;
sz = 1608;

FileName = ['sequential/560_denoised0_' num2str(SampleNum) '.mat'];
load(FileName)
ImagePrim = ImageStack560;


ImageMean = max(ImagePrim,[],3);
figure(1)
imagesc(ImageMean)
colormap gray
caxis([100 300])
axis equal

selectROI    %to select the area you want to do the bead calibartion on 

%ouput of selectROI is roiList containing rect. rect contains x/y and
%distance to next x/y
%coordiantes of the position of the ROI you choose manually whereas x1/y1 is
%the x7y-coordinate of the upper left corner, x2/y2 of the upper right corner,


x1 = round(roiList.rect(1));   %rect(1) is the distance from 0 to first x ROI pixel in x direction
x2 = round(roiList.rect(1)+roiList.rect(3)); %rect(3) is the length of the ROI in x-direction
y1 = round(roiList.rect(2));    %rect(2) is the distance from 0 to first y ROI pixel in y direction
y2 = round(roiList.rect(2)+roiList.rect(4)); %rect(4) is the length of the ROI in y-direction
sz1 = y2 - y1 + 1;          %sz1 resulting ROI in y 
sz2 = x2 - x1 + 1;          %sz2 resulting ROI in x 






%%
FileName = ['sequential/560_denoised_strain_' num2str(SampleNum) '.mat'];
load(FileName)


%%
FileName = ['sequential/647_denoised_strain_' num2str(SampleNum) '.mat'];
load(FileName)



ImageMean560 = max(ImageStack560,[],3);
ImageMean647 = max(ImageStack647,[],3);
ImageMeanPrim = max(ImagePrim,[],3);

figure 
subplot(1,3,1)
imagesc(ImageMeanPrim)
colormap gray
axis square
subplot(1,3,2)
imagesc(ImageMean560)
colormap gray
axis square
subplot(1,3,3)
imagesc(ImageMean647)
colormap gray
axis square

RGB(:,:,1) = ImageMean560/max(max(ImageMean560))*2;
RGB(:,:,2) = ImageMean647/max(max(ImageMean647))*2;
RGB(:,:,3) = ImageMeanPrim/max(max(ImageMeanPrim))*2;
%RGB(find(RGB>1)) = 1;
figure
imagesc(RGB)
axis equal

%%


ImageStack560 = ImageStack560(ceil(y1):ceil(y3),ceil(x1):ceil(x2),1:100);

for i = 1:size(ImageStack560, 3)
    Std560(i) = std2(ImageStack560(:,:,i));
end

ImageStack647 = ImageStack647(ceil(y1):ceil(y3),ceil(x1):ceil(x2),1:100);

for i = 1:size(ImageStack647, 3)
    Std647(i) = std2(ImageStack647(:,:,i));
end


ImagePrim = ImagePrim(ceil(y1):ceil(y3),ceil(x1):ceil(x2),1:100);
for i = 1:size(ImagePrim, 3)
    StdPrim(i) = std2(ImagePrim(:,:,i));
end

GaussEqu = 'a*exp(-(x-b)^2/2/c^2)+d';


StartPoint = [max(StdPrim) 10 5 0];
f_left = fit([1:length(StdPrim)]', StdPrim', GaussEqu, 'Start', StartPoint);

StartPoint = [max(Std560) 5 5 0];
f_right1 = fit([1:length(Std560)]', Std560', GaussEqu, 'Start', StartPoint);


StartPoint = [max(Std647) 5 5 0];
f_right2 = fit([1:length(Std647)]', Std647', GaussEqu, 'Start', StartPoint);



figure(300)
subplot(1,3,1)
plot(f_left, 1:length(StdPrim), StdPrim);
legend('off');
subplot(1,3,2)
plot(f_right1, 1:length(Std560), Std560);
legend('off');
subplot(1,3,3)
plot(f_right2, 1:length(Std647), Std647);
legend('off');



DeltaZ560 = f_right1.b-f_left.b;         %f_left.b is fixed one 
DeltaZ647 = f_right2.b-f_left.b;

save(['DeltaZ_strain_' num2str(SampleNum) '.mat'],'DeltaZ560', 'DeltaZ647');





