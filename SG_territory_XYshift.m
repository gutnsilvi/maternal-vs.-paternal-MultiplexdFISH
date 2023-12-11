%% This script determines the xy shift between imaging rounds of the primary territory and strain specific marking
%
% Input: Hyb0 488 channel with beads and strain Hyb 488 channel with
% beads - shift is calculated by taking at least 2 beads
% saves a tform2 file which can be applied to shift the strain Hyb on top of Hyb0
% --------------------------------------------------------------
% Output: tform2 - Warping matrix from the 647 and 560 strain imaging round
% to 560 primary imaging round
%
% --------------------------------------------------------------
% Silvia Gutnik
% 10.03.2021
% Version 1.0
% 
% using fitFoci and selectROI functions from the Zhuang lab
% --------------------------------------------------------------



clear all 
close all


SampleNum = 18; 

lmol_match =[];
rmol_match =[];
FileName = ['sequential/488_strain_' num2str(SampleNum) '.mat'];
load(FileName)

ImageStackS = ImageStack488;

FileName = ['sequential/488_0_' num2str(SampleNum) '.mat'];
load(FileName)

ImageStack0 = ImageStack488;


ImageMean = max(ImageStack0,[],3);
figure(1)
imagesc(ImageMean)
colormap gray
caxis([100 300])
axis equal

WhetherROI = questdlg('Do you want to select ROIs ?'); %ask question
if(strcmp(WhetherROI, 'Yes'))
selectROI

figure(2)
imagesc(ImageMean)
colormap gray
axis equal
hold on

for j = 1: length(roiList)
    i= 2*j-1;
    [Xfit(j), Yfit(j), Zfit(j)] = fitFoci(ImageStack0, roiList(j), i ,1);
    
    figure(2)
    plot(Xfit(j), Yfit(j), 'x');
    lmol_match = cat(1, lmol_match, [Xfit(j), Yfit(j), Zfit(j)]);
    i = 2*j;   
    [Xfit(j), Yfit(j), Zfit(j)] = fitFoci(ImageStackS, roiList (j), i,1);
    
    figure(2)
    plot(Xfit(j), Yfit(j), 'x');
    rmol_match = cat(1, rmol_match, [Xfit(j), Yfit(j), Zfit(j)]);
end
 
tform2 = cp2tform(rmol_match(:,1:2),lmol_match(:,1:2),'nonreflective similarity');
save(['tform_strain_' num2str(SampleNum) '.mat'],'tform2');

end



%% Visualization


ImageMean488 = max(ImageStack0,[],3);
RGB(:,:,1) = ImageMean488/max(max(ImageMean488))*2;


ImageMean488 = max(ImageStackS,[],3);
RGB(:,:,2) = ImageMean488/max(max(ImageMean488))*2;
RGB(:,:,3) = zeros(1608,1608);
%RGB(find(RGB>1)) = 1;
figure(2)
imagesc(RGB)
axis equal
 
TransImg = imtransform(RGB(:,:,2), tform2, 'XData', [1 1608], 'Ydata', [1 1608]);
RGB(:,:,2) = TransImg;
figure(500)
imagesc(RGB)
axis equal



   