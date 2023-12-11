% This script calbrates the two channels (561 and 647) used for imaging DNA
% loci in sequential FISH
% ------------------------------------------------------------------------------------------------------------------------------------------
% Input: multidimensional matrices (3D image stacks) of Tetraspeck beads
% run SG_beads_cal_import prior to the script
% ------------------------------------------------------------------------------------------------------------------------------------------
% Output: Warping matrix from the 647 channel to the 561 channel (tform + DeltaZ)
% ------------------------------------------------------------------------------------------------------------------------------------------
% Ahilya N. Sawh, PhD
% 23.08.2019
% Version 1.0
% 
% using fitFoci and selectROI functions from the Zhuang lab
% 
% updated by Silvia Gutnik 20.04.2023
% ------------------------------------------------------------


clear all
close all

sz = 1608;

FileName = ['beads/560_beads.mat'];
load(FileName)
for i = 1:size(ImageStack560, 3)
    Std560(i) = std2(ImageStack560(:,:,i));
end

%%
FileName = ['beads/647_beads.mat'];
load(FileName)
for i = 1:size(ImageStack647, 3)
    Std647(i) = std2(ImageStack647(:,:,i));
end
%%
FileName = ['beads/488_beads.mat'];
load(FileName)
for i = 1:size(ImageStack488, 3)
    Std488(i) = std2(ImageStack488(:,:,i));
end

%%
ImageMean560 = max(ImageStack560,[],3);
ImageMean647 = max(ImageStack647,[],3);
ImageMean488 = max(ImageStack488,[],3);

figure 
subplot(1,3,1)
imagesc(ImageMean488)
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
RGB(:,:,3) = zeros(sz,sz);
%RGB(find(RGB>1)) = 1;
figure(2)
imagesc(RGB)
axis equal
%%
WhetherROI = questdlg('Do you want to select ROIs ?'); %ask question
if(strcmp(WhetherROI, 'Yes'))
    selectROI;    
    lmol_match =[];
    rmol_match =[];
    for i = 1:length(roiList)
        [Xfit, Yfit, Zfit] = fitFoci(ImageStack560, roiList(i),2*i-1,1);      
        lmol_match = cat(1, lmol_match, [Xfit, Yfit, Zfit]);
        [Xfit, Yfit, Zfit] = fitFoci(ImageStack647, roiList(i),2*i,1);   
        rmol_match = cat(1, rmol_match, [Xfit, Yfit, Zfit]); %tform = cp2tform(movingPoints,fixedPoints,transformationType)
    end
    tform = cp2tform(rmol_match(:,1:2),lmol_match(:,1:2),'projective');
    save('tform.mat','tform');
end

%%
GaussEqu = 'a*exp(-(x-b)^2/2/c^2)+d';


StartPoint = [max(Std560) 20 5 0];
f_left = fit([1:length(Std560)]', Std560', GaussEqu, 'Start', StartPoint);


StartPoint = [max(Std647) 20 5 0];
f_right = fit([1:length(Std647)]', Std647', GaussEqu, 'Start', StartPoint);

figure(300)
subplot(1,2,1)
plot(f_left, 1:length(Std560), Std560);
legend('off');
subplot(1,2,2)
plot(f_right, 1:length(Std647), Std647);
legend('off');

DeltaZ = f_right.b-f_left.b;

save('DeltaZ.mat','DeltaZ');

TransImg = imtransform(RGB(:,:,2), tform, 'XData', [1 sz], 'Ydata', [1 sz]);
RGB(:,:,2) = TransImg;
figure(500)
imagesc(RGB)
axis equal



