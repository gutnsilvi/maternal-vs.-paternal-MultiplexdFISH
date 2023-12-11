% This script opens and converts the bead calibration-file from nd2 format to multidemensional matrices for downstream
% processing
% ------------------------------------------------------------------------------------------------------------------------------------------
% Requirements: bfmatlab installed and in the path
% 
% Silvia Gutnik
% 15.01.2020
% Version 1.0 
% based on the script ANSrun_initialize_nd2files.m by Ahilya N. Sawh
% 
%% ----------------------------------------------------------------

clear all
close all

mkdir beads
sz = 1608;  %pixel size of images 
z = 41;     %number of z-sclices used for bead calibration
c = 4;      %number of channels imaged
zc = z*c;


%bead calibartion file import bioformates

data = bfopen(['beadCal054.nd2']);
ImageStack = zeros(sz,sz,zc);
for i = 1:zc %planes are organized by z position first then channel
    ImageStack(:,:,i) = data{1, 1}{i,1};
end

%channel 2 = 488
ImageStack488 = zeros(sz,sz,z);
ImageStack488 = ImageStack(:,:,2:4:zc);
ImageMax = max(ImageStack488,[],3);
figure
imagesc(ImageMax)
axis equal
colormap gray
title(['488_beads'])
saveas(gcf,['beads/488_beads'])
save( ['beads/488_beads.mat'],'ImageStack488');
clear ImageMax ImageStack488

%channel 3 = 560
ImageStack560 = zeros(sz,sz,z);
ImageStack560 = ImageStack(:,:,3:4:zc);
ImageMax = max(ImageStack560,[],3);
figure
imagesc(ImageMax)
axis equal
colormap gray
title(['560_beads'])
saveas(gcf,['beads/560_beads'])
save( ['beads/560_beads.mat'],'ImageStack560');
clear ImageMax ImageStack560

%channel 4 = 647
ImageStack647 = zeros(sz,sz,z);
ImageStack647 = ImageStack(:,:,4:4:zc);
ImageMax = max(ImageStack647,[],3);
figure
imagesc(ImageMax)
axis equal
colormap gray
title(['647_beads'])
saveas(gcf,['beads/647_beads'])
save( ['beads/647_beads.mat'],'ImageStack647');
clear ImageMax ImageStack647

%channel 1 = 405
ImageStack405 = zeros(sz,sz,z);
ImageStack405 = ImageStack(:,:,1:4:zc);
ImageMax = max(ImageStack405,[],3);
figure
imagesc(ImageMax)
axis equal
colormap gray
title(['405_beads'])
saveas(gcf,['beads/405_beads'])
save( ['beads/405_beads.mat'],'ImageStack405');
clear ImageStack ImageMax ImageStack405
close all

clear data
