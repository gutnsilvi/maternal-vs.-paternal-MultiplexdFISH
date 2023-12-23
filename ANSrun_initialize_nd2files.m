% This script opens and converts files from nd2 format to multidemensional matrices for downstream
% processing
% ------------------------------------------------------------------------------------------------------------------------------------------
% Requirements: bfmatlab installed and in the path
% ------------------------------------------------------------------------------------------------------------------------------------------
% Input: path to .nd2 file with multiple series(FOV/samples) and illumination channels
% 
% series: FOV/position
% 
% zstacks for sequential FISH hybs are numbered starting at 0 in increments
% of 3 because a snap image is taken before and after bleaching steps
% ie Hyb0 = Seq0000
%    Hyb1 = Seq0003
%    Hyb2 = Seq0006
%    Hyb3 = Seq0009
%    Hyb4 = Seq00012
%    Hyb5 = Seq00015
% ------------------------------------------------------------------------------------------------------------------------------------------
% Output: Multidemensional matrices for each channel, series, and hyb#
% ------------------------------------------------------------------------------------------------------------------------------------------
% Ahilya N. Sawh, PhD
% 06.03.2019
% Version 1.0
%% ------------------------------------------------------------------------------------------------------------------------------------------
clear all 
close all
tic

mkdir sequential

NumSeries = 20;
NumHybs = 11;
sz = 1608;

%% primary probe file

for Hyb = 0
    Seq = Hyb*3;
    
    data = bfopen(['Channel560,488,405_Seq000' num2str(Seq) '.nd2']);
    
    for series = 1:NumSeries
        ImageStack = zeros(sz,sz,453);
        for i = 1:453 %planes are organized by z position first then channel
            ImageStack(:,:,i) = data{series, 1}{i,1};
        end
        
        %channel 2 = 488
        ImageStack488 = zeros(sz,sz,151);
        ImageStack488 = ImageStack(:,:,2:3:453);
        ImageMax = max(ImageStack488,[],3);
        figure
        imagesc(ImageMax)
        axis equal
        colormap gray
        title(['488_Hyb' num2str(Hyb) '_FOV' num2str(series)])
        saveas(gcf,['sequential/488_Hyb' num2str(Hyb) '_FOV' num2str(series)])
        save( ['sequential/488_' num2str(Hyb) '_' num2str(series) '.mat'],'ImageStack488','-v7.3');
        clear ImageMax ImageStack488
        
        %channel 1 = 560
        ImageStack560 = zeros(sz,sz,151);
        ImageStack560 = ImageStack(:,:,1:3:453);
        ImageMax = max(ImageStack560,[],3);
        figure
        imagesc(ImageMax)
        axis equal
        colormap gray
        title(['560_Hyb' num2str(Hyb) '_FOV' num2str(series)])
        saveas(gcf,['sequential/560_Hyb' num2str(Hyb) '_FOV' num2str(series)])
        save( ['sequential/560_' num2str(Hyb) '_' num2str(series) '.mat'],'ImageStack560','-v7.3');
        clear ImageMax ImageStack560
        
   
        %channel 3 = 405
        ImageStack405 = zeros(sz,sz,151);
        ImageStack405 = ImageStack(:,:,3:3:453);
        ImageMax = max(ImageStack405,[],3);
        figure
        imagesc(ImageMax)
        axis equal
        colormap gray
        title(['405_Hyb' num2str(Hyb) '_FOV' num2str(series)])
        saveas(gcf,['sequential/405_Hyb' num2str(Hyb) '_FOV' num2str(series)])
        save( ['sequential/405_' num2str(Hyb) '_' num2str(series) '.mat'],'ImageStack405','-v7.3');
        clear ImageStack ImageMax ImageStack405
        close all
    end
    clear data
end

%% secondary hyb files
tic

for Hyb = 1:NumHybs 
    Seq = Hyb*3;
        if Seq <= 10
            data = bfopen(['Channel560,488,647_Seq000' num2str(Seq) '.nd2']);
        elseif Seq >= 10
            data = bfopen(['Channel560,488,647_Seq00' num2str(Seq) '.nd2']);
        end
        for series = 1:NumSeries
            ImageStack = zeros(sz,sz,453);
            for i = 1:453 %planes are organized by z position first then channel
                ImageStack(:,:,i) = data{series, 1}{i,1};
            end
            
            %channel 2 = 488
            ImageStack488 = zeros(sz,sz,151);
            ImageStack488 = ImageStack(:,:,2:3:453);
            ImageMax = max(ImageStack488,[],3);
            figure
            imagesc(ImageMax)
            axis equal
            colormap gray
            title(['488_Hyb' num2str(Hyb) '_FOV' num2str(series)])
            saveas(gcf,['sequential/488_Hyb' num2str(Hyb) '_FOV' num2str(series)])
            save( ['sequential/488_' num2str(Hyb) '_' num2str(series) '.mat'],'ImageStack488','-v7.3');
            clear ImageMax ImageStack488
            
            %channel 1 = 560
            ImageStack560 = zeros(sz,sz,151);
            ImageStack560 = ImageStack(:,:,1:3:453);
            ImageMax = max(ImageStack560,[],3);
            figure
            imagesc(ImageMax)
            axis equal
            colormap gray
            title(['560_Hyb' num2str(Hyb) '_FOV' num2str(series)])
            saveas(gcf,['sequential/560_Hyb' num2str(Hyb) '_FOV' num2str(series)])
            save( ['sequential/560_' num2str(Hyb) '_' num2str(series) '.mat'],'ImageStack560','-v7.3');
            clear ImageMax ImageStack560
            
            %channel 3 = 647
            ImageStack647 = zeros(sz,sz,151);
            ImageStack647 = ImageStack(:,:,3:3:453);
            ImageMax = max(ImageStack647,[],3);
            figure
            imagesc(ImageMax)
            axis equal
            colormap gray
            title(['647_Hyb' num2str(Hyb) '_FOV' num2str(series)])
            saveas(gcf,['sequential/647_Hyb' num2str(Hyb) '_FOV' num2str(series)])
            save( ['sequential/647_' num2str(Hyb) '_' num2str(series) '.mat'],'ImageStack647','-v7.3');
            clear ImageMax ImageStack647 ImageStack
            
            close all
        end
    clear data 
end

toc