% This script traces chromosomes by nearest neighbour method, using
% watershed segmentation of nuclei and primary probe signal
% ------------------------------------------------------------------------------------------------------------------------------------------
% Input: multidimensional matrices (3D image stacks) of nuclear signal,
% primary probe signal, DeltaZ, tform, and fitted foci result file for each sample
% 
% create folders before running: traces
% run section by section for each sample
%
% ------------------------------------------------------------------------------------------------------------------------------------------
% Output: Traces with coordinates of each FOV
% ------------------------------------------------------------------------------------------------------------------------------------------
% using traceChromosome_3D_L1 function 
% 
% Ahilya N. Sawh, PhD
% 24.09.2020
% adapted by Silvia Gutnik 28.10.2020
% 
%% ------------------------------------------------------------------------------------------------------------------------------------------
clear all
close all

% set parameters
global SampleNum
global flap
global x1
global x2
global y1
global y3

SampleNum = '20';

FOV = 1; %which embryo are we analyzing now (if > 1 in the sample)

save_traces = 1;

sz = 1608;

nuc = 1;    %1 = take into account nuclear Vol for segmentation 0 = only in rare cases of bad nuclear stain 0

save_age = 1;
man_age = 0;    

%%%%%%%%%%maxintensitythreshs%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Thresh405 = 3000;
Thresh560 = 300;

zeroThre560 = 90;
zeroThre405 = 1500;

TerrDil560 = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%for reconstruction
dilateThresh = 5; % 
erodeThresh = 5; % 

%for binarization and watershed strength%%%change terrmin if traces are too
%many/few per territory

TerrbinaryThresh = 10; 
TerrMin = 0.15; 
NucbinaryThresh = 300; 
NucMin= 0.45; 

TotalTADNum = 22;
TADsToExclude =[];

Color = colormap(jet(28));

%% find and segment nuclei
FileName = ['sequential/405_denoised0_' SampleNum '.mat'];

load(FileName)

ImDAPI = max(ImageStack405,[],3);
figure
imagesc(ImDAPI)
colormap gray
axis equal
title('Z stack max image')%flattened image of nuclei full xyz stack

% select an embryo FOV in the sample, perform all ds actions in this FOV
selectROI

I = ImageStack405(ceil(y1):ceil(y3),ceil(x1):ceil(x2),:);
figure
ImDAPI = max(I,[],3);
imagesc(ImDAPI)
colormap gray
axis equal
caxis([0 10000])
title('selected FOV 405')

%correct for uneven background illumination

se = strel('disk',50);
background = imopen(I,se);
figure
imagesc(max(background,[],3))
colormap gray
axis equal
caxis([0 3500])
title('selected FOV primary probe 405 bkg')

I2 = I - background;
figure
imagesc(max(I2,[],3))
colormap gray
axis equal
caxis([0 2000])
title('selected FOV primary probe 405 minus bkg')

%set min and max values to even out territory intensity

I3 = I2;
I3(I3<zeroThre405) = 0;
I3(I3>Thresh405) = Thresh405;
figure
imagesc(max(I3,[],3))
colormap gray
axis equal
caxis([0 1000])
title('selected FOV 405 thresh')

%3D binary
BW = imbinarize(I3, NucbinaryThresh);
nucBW = imfill(BW, 'holes');
figure
imshowpair(max(I3,[],3),max(nucBW,[],3),'montage')
title('before and after 3D binary')

%3D dist transform
D = -bwdist(~nucBW); 

%set background to its own catchment basin - i.e. where BW is black
D(~nucBW) = -Inf;

%3D watershed
D = imhmin(D,NucMin); %the height threshold for suppressing shallow minima
L = watershed(D);
figure
imshow(label2rgb(L(:,:,50),'jet','w'))

L = imdilate(L,true(5));
nucL = L-1;
clear L

Lmax = max(nucL,[],3);
nucLrgb = label2rgb(Lmax,'jet','w','shuffle');
figure
imshow(nucLrgb)
title('colored watershed label matrix nuclei')

figure
imagesc(ImDAPI)
colormap gray
hold on
himage = imshow(nucLrgb);
himage.AlphaData = 0.3;
title('nuclei labels superimposed transparently on DAPI')

age = max(max(max(nucL)))

nstats = regionprops3(nucL, 'all');
ncenters = nstats.Centroid;




%% segment based on primary probe 561

FileName = ['sequential/560_denoised0_' SampleNum];
load(FileName)
 

Im560 = max(ImageStack560,[],3);
figure
imagesc(Im560)
colormap gray
axis equal
title('z stack max image 561 primary ')
caxis([100 700])


I = ImageStack560(ceil(y1):ceil(y3),ceil(x1):ceil(x2),:);

%correct for uneven background illumination
Im560 = max(I,[],3);
figure
imagesc(Im560)
colormap gray
axis equal
caxis([0 1000])
title('selected FOV primary probe 560')

se = strel('disk',10);
background = imopen(I,se);
figure
imagesc(max(background,[],3))
colormap gray
axis equal
caxis([0 1000])
title('selected FOV primary probe 560 bkg')

I2 = I - background;
figure
imagesc(max(I2,[],3))
colormap gray
axis equal
caxis([0 500])
title('selected FOV primary probe 560 minus bkg')

%set min and max values to even out territory intensity

I3 = I2;
I3(I3< zeroThre560) = 0;
I3(I3>Thresh560) = Thresh560;
figure
imagesc(max(I3,[],3))
colormap gray
axis equal
caxis([0 400])
title('selected FOV primary probe 530 thresh')



se = strel('disk', erodeThresh);

Icbre = imerode(I3, se);

Icbrobr = imreconstruct(Icbre, I3);
figure
imagesc(max(Icbrobr,[],3))
colormap gray
axis equal
title('closing and opening by reconstruction');


%3D binary
BW = imbinarize(Icbrobr, TerrbinaryThresh);
BW = imfill(BW, 'holes');
BW = imclose(BW, se);
figure
imshowpair(max(Icbrobr,[],3),max(BW,[],3),'montage')
title('before and after 3Dbinary')

BW1 = BW;

se = strel('sphere',TerrDil560);
dilatedBW = imdilate(BW1,se);

%3D dist transform
D = -bwdist(~dilatedBW); 

%set background to its own catchment basin - where BW is 0
D(~dilatedBW) = -Inf;

%3D watershed
D = imhmin(D,TerrMin); %the height threshold for suppressing shallow minima
L = watershed(D);
figure
imshow(label2rgb(L(:,:,50),'jet','w'))
title('watershed slice 50')


L = imdilate(L,true(5));
L = L-1; %remove background catchment basin - set to 0

%%if pixels were bkg in nuclei label matrix and
%%if yes convert those pixels to bkg in territory label matrix (zero)

if nuc == 1 
    L(find(~nucL))=0;
end

Lmax = max(L,[],3);
Lrgb = label2rgb(Lmax,'jet','w','shuffle');
imshow(Lrgb)
title('colored watershed label matrix territory 560')

figure
imagesc(Im560)
colormap gray
hold on
himage = imshow(Lrgb);
himage.AlphaData = 0.3;
title('territory labels superimposed transparently on 560')

numterr = max(max(max(L)))

L1 = L;
clear L 


 
L1max = max(L1,[],3);
L1rgb = label2rgb(L1max,'jet','w','shuffle');
figure
imshow(L1rgb)
title('territories 560')
axis equal


%% load and crop results file to FOV
load(['result' num2str(SampleNum) '.mat']);

if ~isempty(TADsToExclude)
    for i = TADsToExclude 
        Xfit{1,i} = [];
        Yfit{1,i} = [];
        Zfit{1,i} = [];
        Xgof{1,i} = [];
        Ygof{1,i} = [];
        Zgof{1,i} = [];
        Intensity{1,i} = [];
    end
end

%convert to pixel
for i = 1:length(Intensity)
    Yfit{i} = sz*110/1000-Yfit{i}; 
    Yfit{i} = Yfit{i}/110*1000; %pxl
    Xfit{i} = Xfit{i}/110*1000; %pxl
    Zfit{i} = Zfit{i}/0.2; % z stack no.
end

figure
hold on
for i = 1:length(Intensity)
    scatter3(Xfit{i}, Yfit{i}, Zfit{i}, 'ok', 'MarkerFaceColor', Color (i,:));
    text(Xfit{i}+0.05, Yfit{i}, Zfit{i}, num2str(i), 'Color', 'Black', ...
                'FontWeight', 'Bold', 'FontSize',12); %write numbers on the image
end
hold off
axis equal
xlabel('x');
ylabel('y');
zlabel('z');
set(gca,'YDir','reverse')
title('All foci - pxl')

Xfitnew = cell(1,TotalTADNum);
Yfitnew = cell(1,TotalTADNum);
Zfitnew = cell(1,TotalTADNum);
Xgofnew = cell(1,TotalTADNum);
Ygofnew = cell(1,TotalTADNum);
Zgofnew = cell(1,TotalTADNum);
Intensitynew = cell(1,TotalTADNum);
for i=1:TotalTADNum
    q = 0;
    if ~isempty(Xfit{1,i})
        for j=1:length(Xfit{1,i})
            if Xfit{1,i}(1,j) >= x1 && Xfit{1,i}(1,j) <= x2 && Yfit{1,i}(1,j) >= y1 && Yfit{1,i}(1,j) <= y3 %inside FOV coordinates
                q = q+1;
                Xfitnew{1,i}(1,q) = Xfit{1,i}(1,j);
                Yfitnew{1,i}(1,q) = Yfit{1,i}(1,j);
                Zfitnew{1,i}(1,q) = Zfit{1,i}(1,j);
                Xgofnew{1,i}(1,q) = Xgof{1,i}(1,j);
                Ygofnew{1,i}(1,q) = Ygof{1,i}(1,j);
                Zgofnew{1,i}(1,q) = Zgof{1,i}(1,j);
                Intensitynew{1,i}(1,q) = Intensity{1,i}(1,j);
            end
        end
    end
end

clear Xfit Yfit Zfit Xgof Ygof Zgof Intensity

Xfit = Xfitnew;
Yfit = Yfitnew;
Zfit = Zfitnew;
Xgof = Xgofnew;
Ygof = Ygofnew;
Zgof = Zgofnew;
Intensity = Intensitynew;

figure
hold on
for i = 1:length(Intensity)
    scatter3(Xfit{i}, Yfit{i}, Zfit{i}, 'ok', 'MarkerFaceColor', Color (i,:));
    text(Xfit{i}+0.05, Yfit{i}, Zfit{i}, num2str(i), 'Color', 'Black', ...
                'FontWeight', 'Bold', 'FontSize',12); %write numbers on the image
end
hold off
axis equal
xlabel('x');
ylabel('y');
zlabel('z');
set(gca,'YDir','reverse')
title('All foci in FOV')

%zero the foci locations to the cropped FOV in x and y so that they align
%with the cropped image from prev steps
Xfit_new = cell(1,TotalTADNum);
Yfit_new = cell(1,TotalTADNum);
for i=1:TotalTADNum
    for j=1:length(Xfit{1,i})
        Xfit_new{1,i} = Xfit{1,i} - x1;
        Yfit_new{1,i} = Yfit{1,i} - y1;
    end
end

figure
imagesc(Im560)
colormap gray
hold on

h1 = patch(isosurface(L1,0.5));
h1.EdgeColor = 'none';
h1.FaceColor = [1,0,0];
h1.FaceAlpha = 0.25;

h2 = patch(isosurface(nucL,0.5));
h2.EdgeColor = 'none';
h2.FaceColor = [0,1,0];
h2.FaceAlpha = 0.1;


for i = 1:length(Intensity)
    scatter3(Xfit_new{i}, Yfit_new{i}, Zfit{i}, 'ok', 'MarkerFaceColor', Color (i,:));
    text(Xfit{i}+0.05, Yfit{i}, Zfit{i}, num2str(i), 'Color', 'Black', ...
                'FontWeight', 'Bold', 'FontSize',12); %write numbers on the image
end
hold off
axis equal
xlabel('x');
ylabel('y');
zlabel('z');
set(gca,'YDir','reverse')
title('All foci in FOV zeroed')

%% trace chromosomes inside 560 territories (L1)


figure
imagesc(Im560)
colormap gray
hold on
himage = imshow(nucLrgb);
himage.AlphaData = 0.1;
himage = imshow(L1rgb);
himage.AlphaData = 0.25;



for i = 1:max(max(max(L1)))
   
    [Trace, TAD_id] = traceChromosome_3D_L1(L1, i, Xfit_new, Yfit_new, Zfit, Intensity);
    if flap ~= 1
        if numel(Trace) == 1
            scatter3(Trace{1,1}(:,1), Trace{1,1}(:,2), Trace{1,1}(:,3),'or', 'MarkerFaceColor', 'r')
            plot3(Trace{1,1}(:,1), Trace{1,1}(:,2), Trace{1,1}(:,3),'-r');
            text(Trace{1,1}(:,1), Trace{1,1}(:,2),Trace{1,1}(:,3), num2str(TAD_id{1}), 'Color', 'Black', ...
                'FontWeight', 'Bold', 'FontSize',8);
        elseif numel(Trace) == 2
            scatter3(Trace{1,1}(:,1), Trace{1,1}(:,2), Trace{1,1}(:,3),'or', 'MarkerFaceColor', 'r')
            plot3(Trace{1,1}(:,1), Trace{1,1}(:,2), Trace{1,1}(:,3),'-r');
            text(Trace{1,1}(:,1), Trace{1,1}(:,2),Trace{1,1}(:,3), num2str(TAD_id{1}), 'Color', 'Black', ...
                'FontWeight', 'Bold', 'FontSize',8);
            scatter3(Trace{1,2}(:,1), Trace{1,2}(:,2), Trace{1,2}(:,3),'sb', 'MarkerFaceColor', 'b')
            plot3(Trace{1,2}(:,1), Trace{1,2}(:,2), Trace{1,2}(:,3),'-b');
            text(Trace{1,2}(:,1), Trace{1,2}(:,2),Trace{1,2}(:,3), num2str(TAD_id{2}), 'Color', 'Black', ...
                'FontWeight', 'Bold', 'FontSize',8);
        elseif numel(Trace) == 3
            scatter3(Trace{1,1}(:,1), Trace{1,1}(:,2), Trace{1,1}(:,3),'or', 'MarkerFaceColor', 'r')
            plot3(Trace{1,1}(:,1), Trace{1,1}(:,2), Trace{1,1}(:,3),'-r');
            text(Trace{1,1}(:,1), Trace{1,1}(:,2),Trace{1,1}(:,3), num2str(TAD_id{1}), 'Color', 'Black', ...
                'FontWeight', 'Bold', 'FontSize',8);
            scatter3(Trace{1,2}(:,1), Trace{1,2}(:,2), Trace{1,2}(:,3),'sb', 'MarkerFaceColor', 'b')
            plot3(Trace{1,2}(:,1), Trace{1,2}(:,2), Trace{1,2}(:,3),'-b');
            text(Trace{1,2}(:,1), Trace{1,2}(:,2),Trace{1,2}(:,3), num2str(TAD_id{2}), 'Color', 'Black', ...
                'FontWeight', 'Bold', 'FontSize',8);
            scatter3(Trace{1,3}(:,1), Trace{1,3}(:,2), Trace{1,3}(:,3),'dg', 'MarkerFaceColor', 'g')
            plot3(Trace{1,3}(:,1), Trace{1,3}(:,2), Trace{1,3}(:,3),'-g');
            text(Trace{1,3}(:,1), Trace{1,3}(:,2),Trace{1,3}(:,3), num2str(TAD_id{3}), 'Color', 'Black', ...
                'FontWeight', 'Bold', 'FontSize',8);
        elseif numel(Trace) == 4
            scatter3(Trace{1,1}(:,1), Trace{1,1}(:,2), Trace{1,1}(:,3),'or', 'MarkerFaceColor', 'r')
            plot3(Trace{1,1}(:,1), Trace{1,1}(:,2), Trace{1,1}(:,3),'-r');
            text(Trace{1,1}(:,1), Trace{1,1}(:,2),Trace{1,1}(:,3), num2str(TAD_id{1}), 'Color', 'Black', ...
                'FontWeight', 'Bold', 'FontSize',8);
            scatter3(Trace{1,2}(:,1), Trace{1,2}(:,2), Trace{1,2}(:,3),'sb', 'MarkerFaceColor', 'b')
            plot3(Trace{1,2}(:,1), Trace{1,2}(:,2), Trace{1,2}(:,3),'-b');
            text(Trace{1,2}(:,1), Trace{1,2}(:,2),Trace{1,2}(:,3), num2str(TAD_id{2}), 'Color', 'Black', ...
                'FontWeight', 'Bold', 'FontSize',8);
            scatter3(Trace{1,3}(:,1), Trace{1,3}(:,2), Trace{1,3}(:,3),'dg', 'MarkerFaceColor', 'g')
            plot3(Trace{1,3}(:,1), Trace{1,3}(:,2), Trace{1,3}(:,3),'-g');
            text(Trace{1,3}(:,1), Trace{1,3}(:,2),Trace{1,3}(:,3), num2str(TAD_id{3}), 'Color', 'Black', ...
                'FontWeight', 'Bold', 'FontSize',8);
            scatter3(Trace{1,4}(:,1), Trace{1,4}(:,2), Trace{1,4}(:,3),'^k', 'MarkerFaceColor', 'k')
            plot3(Trace{1,4}(:,1), Trace{1,4}(:,2), Trace{1,4}(:,3),'-k');
            text(Trace{1,4}(:,1), Trace{1,4}(:,2),Trace{1,4}(:,3), num2str(TAD_id{4}), 'Color', 'Black', ...
                'FontWeight', 'Bold', 'FontSize',8);
        end
    end
    if save_traces == 1 && save_age ==0
        save(['traces/' num2str(man_age) '_cell_Traces_Sample' SampleNum '_FOV' num2str(FOV) '_roi561_' num2str(i) '.mat'],...
            'Trace', 'TAD_id', 'L1', 'nucL', 'TADsToExclude', 'x1', 'x2', 'y1', 'y3');
    end
    if save_traces == 1 && save_age ==1
        save(['traces/' num2str(age) '_cell_Traces_Sample' SampleNum '_FOV' num2str(FOV) '_roi561_' num2str(i) '.mat'],...
            'Trace', 'TAD_id', 'L1', 'nucL', 'TADsToExclude', 'x1', 'x2', 'y1', 'y3');
    end
end
title('traces on segmentation stain 561 - per Terr')
axis equal
hold off


