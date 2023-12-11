% This script traces chromosomes by nearest neighbour method, using
% watershed segmentation of nuclei and primary probe signal
% ------------------------------------------------------------------------------------------------------------------------------------------
% Input: multidimensional matrices (3D image stacks) of nuclear signal,
% primary probe signal, strain-specific territories, DeltaZ, tform, tform2, DeltaZ560, DeltaZ647 and fitted foci result file for each sample
% 
% create folders before running: traces
% 
% ------------------------------------------------------------------------------------------------------------------------------------------
% Output: Traces with coordinates of each FOV
% saved file contains traces sorted by nucleus and HI and N2, as well as
% identity of nuc (segmentation label), TAD_id and strain_ID (is the Foci within HI or N2 territory) for each fitted foci. 
% ------------------------------------------------------------------------------------------------------------------------------------------
% Silvia Gutnik
% 28.04.2021 
% based on the original script obtained ANS_generated_tight_traces_primary_only

% 
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

SampleNum = '18';

FileName = ['tform_strain_' num2str(SampleNum) '.mat']; 

load(FileName) 

FileName = ['DeltaZ_strain_' num2str(SampleNum) '.mat'];
load(FileName)



load('tform.mat');

FOV = 1; %which embryo are we analyzing now (if > 1 in the sample)


nuc = 1;     %1 = take into account nuclear Vol for segmentation 0 = only in rare cases of bad nuclear stain 0
div = 20;    %threshold fraction for nuclear size exclusion (automatic removal of big and small segmented background)


save_traces = 0;   

save_age = 1;
man_age = 0;

%for primary probe segmentation
primarytosegment_561 = 1;
primarytosegment_647 = 1;

%%%%%%%%%%maxintensitythreshs%%%%%%%%%%%%%%%%%%%%%
Thresh405 = 3000;
Thresh560 = 70;
Thresh560T = 300;
Thresh647 = 70;


zeroThre560 = 10;        %strain-label HI
zeroThre560T = 80;       %primary territory

zeroThre647 = 10;        %strain-label N2
zeroThre405 = 150;      %nuclear stain

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TerrDil560T = 1;        
TerrDil560 = 2; 
TerrDil647 = 1;

%for reconstruction
dilateThresh = 5; % 
erodeThresh = 5; % 

%for binarization and watershed strength%%%change terrmin if traces are too
%many/few per territory

TerrbinaryThresh560 = 10;
TerrbinaryThresh647 = 10;
TerrbinaryThreshT = 10;
TerrMin560 = 0.45;
TerrMin647 = 0.45;
TerrMin = 0.35;
NucbinaryThresh = 1000; 
NucMin= 0.40; 

TotalTADNum = 22;
TADsToExclude =[];

Color = colormap(jet(22));

% find and segment nuclei
FileName = ['sequential/405_denoised0_' SampleNum '.mat'];

load(FileName)

ImDAPI = max(ImageStack405,[],3);
figure
imagesc(ImDAPI)
colormap gray
axis equal
caxis([0 10000])
title('Z stack max image')%flattened image of nuclei full xyz stack

% select an embryo FOV in the sample, perform all ds actions in this FOV
selectROI

I = ImageStack405(ceil(y1):ceil(y3),ceil(x1):ceil(x2),:);
figure
ImDAPI = max(I,[],3);
imagesc(ImDAPI)
colormap gray
axis equal
caxis([0 2000])
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
caxis([0 5000])
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


%
%find small nuclear segments like polar bodies and exclude from analysis -
%also find big segmented volumes if embryos are close (happens
%occasionally)

nstats = regionprops3(nucL, 'all');
nstats( any(ismissing(nstats),2), :) = [];
ncenters = nstats.Centroid;
nVolumes = nstats.Volume;


nuc_ID = unique(nucL(:));
idx = unique(nucL(:)) >0;
nuc_ID = nuc_ID(idx);


A = sort(nVolumes)

for i = 1
    
    if  max(nVolumes) > 4* A(end-i) & length(nuc_ID) > 4
        T = round(A(end-i)/4);
        S = A(end-i);
        break;
    else
        T = round(max(nVolumes)/div);
    end
end

TF = find(nVolumes < T);
indexTF = nVolumes < T;


if length(TF) > 0 
    for i = 1:length(indexTF)
        if indexTF(i) > 0
           
            nucL(find(nucL == nuc_ID(i)))= 0; 
        end
    end
end

if exist('S') == 1
    
    TF = find(nVolumes > S);
    indexTF = nVolumes > S
    if length(TF) > 0 
        for i = 1:length(indexTF)
            if indexTF(i) > 0
           
            nucL(find(nucL == nuc_ID(i) ))= 0; 
            end
        end
    end
end
    

Lmax = max(nucL,[],3);
nucLrgb = label2rgb(Lmax,'jet','w','shuffle');
figure
imagesc(ImDAPI)
colormap gray
hold on
himage = imshow(nucLrgb);
himage.AlphaData = 0.3;
title('nuclei labels superimposed transparently on DAPI - outlier')


nuc_ID = unique(nucL(:));
idx = unique(nucL(:)) >0;
nuc_ID = nuc_ID(idx);

age = length(nuc_ID)


%segment 560 territory channel and trace



FileName = ['sequential/560_denoised0_' SampleNum];
load(FileName)

Im560T = max(ImageStack560,[],3);
figure
imagesc(Im560T)
colormap gray
axis equal
title('z stack max image 561 primary')
caxis([100 700])


I = ImageStack560(ceil(y1):ceil(y3),ceil(x1):ceil(x2),:);

%correct for uneven background illumination
Im560T = max(I,[],3);
figure
imagesc(Im560T)
colormap gray
axis equal
caxis([0 1000])
title('selected FOV primary probe 560T')

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
I3(I3< zeroThre560T) = 0;
I3(I3>Thresh560T) = Thresh560T;
figure
imagesc(max(I3,[],3))
colormap gray
axis equal
caxis([0 400])
title('selected FOV primary probe 560T thresh')


se = strel('disk', erodeThresh);

Icbre = imerode(I3, se);

Icbrobr = imreconstruct(Icbre, I3);
figure
imagesc(max(Icbrobr,[],3))
colormap gray
axis equal
title('closing and opening by reconstruction');



%3D binary
BW = imbinarize(Icbrobr, TerrbinaryThreshT);
BW = imfill(BW, 'holes');
BW = imclose(BW, se);
figure
imshowpair(max(Icbrobr,[],3),max(BW,[],3),'montage')
title('before and after 3Dbinary')

BW1 = BW;

se = strel('sphere',TerrDil560T);
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

% figure
% heatmap(L(:,:,50))
% title('watershed slice 50')

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
title('colored watershed label matrix primary 560')

figure
imagesc(Im560T)
colormap gray
hold on
himage = imshow(Lrgb);
himage.AlphaData = 0.3;
title('territory labels superimposed transparently on primary 560')



L4 = L;
clear L 


terr_ID = unique(L4(:));
idx = unique(L4(:)) >0;
terr_ID = terr_ID(idx);

numterr = length(terr_ID)

L4max = max(L4,[],3);
L4rgb = label2rgb(L4max,'jet','w','shuffle');
figure
imshow(L4rgb)
title('territories 560 primary')
axis equal

figure
imagesc(Im560T)
colormap gray
hold on
himage = imshow(nucLrgb);
himage.AlphaData = 0.1;
himage = imshow(L4rgb);
himage.AlphaData = 0.25;
%
%%
% load and crop results file to FOV
load(['resultC_new_' num2str(SampleNum) '.mat']);

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
    Yfit{i} = 1608*110/1000-Yfit{i}; 
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
imagesc(Im560T)
colormap gray
hold on

h1 = patch(isosurface(L4,0.5));
h1.EdgeColor = 'none';
h1.FaceColor = [1,0,0];
h1.FaceAlpha = 0.20;


h3 = patch(isosurface(nucL,0.5));
h3.EdgeColor = 'none';
h3.FaceColor = [0,0,1];
h3.FaceAlpha = 0.2;


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
legend('Terr','nucL')

%

figure
imagesc(Im560T)
colormap gray
hold on
himage = imshow(nucLrgb);
himage.AlphaData = 0.1;
himage = imshow(L4rgb);
himage.AlphaData = 0.25;


combined_Traces = cell(length(terr_ID),4);
combined_TAD_id = cell(length(terr_ID),4);

for z = 1:length(terr_ID)
    i = terr_ID(z);
    

    [Trace, TAD_id] = traceChromosome_3D_L1(L4, i, Xfit_new, Yfit_new, Zfit, Intensity);
    if flap ~= 1
        if numel(Trace) == 1
            scatter3(Trace{1,1}(:,1), Trace{1,1}(:,2), Trace{1,1}(:,3),'or', 'MarkerFaceColor', 'r')
            plot3(Trace{1,1}(:,1), Trace{1,1}(:,2), Trace{1,1}(:,3),'-r');
            text(Trace{1,1}(:,1), Trace{1,1}(:,2),Trace{1,1}(:,3), num2str(TAD_id{1}), 'Color', 'Black', ...
                'FontWeight', 'Bold', 'FontSize',8);
            combined_Traces(z,1) = Trace(1,1);
            combined_TAD_id(z,1) = TAD_id(1);
            
            
        elseif numel(Trace) == 2
            scatter3(Trace{1,1}(:,1), Trace{1,1}(:,2), Trace{1,1}(:,3),'or', 'MarkerFaceColor', 'r')
            plot3(Trace{1,1}(:,1), Trace{1,1}(:,2), Trace{1,1}(:,3),'-r');
            text(Trace{1,1}(:,1), Trace{1,1}(:,2),Trace{1,1}(:,3), num2str(TAD_id{1}), 'Color', 'Black', ...
                'FontWeight', 'Bold', 'FontSize',8);
            scatter3(Trace{1,2}(:,1), Trace{1,2}(:,2), Trace{1,2}(:,3),'sb', 'MarkerFaceColor', 'b')
            plot3(Trace{1,2}(:,1), Trace{1,2}(:,2), Trace{1,2}(:,3),'-b');
            text(Trace{1,2}(:,1), Trace{1,2}(:,2),Trace{1,2}(:,3), num2str(TAD_id{2}), 'Color', 'Black', ...
                'FontWeight', 'Bold', 'FontSize',8);
            combined_Traces(z,1) = Trace(1,1);
            combined_TAD_id(z,1) = TAD_id(1);
            combined_Traces(z,2) = Trace(1,2);
            combined_TAD_id(z,2) = TAD_id(2);
            
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
            
            combined_Traces(z,1) = Trace(1,1);
            combined_TAD_id(z,1) = TAD_id(1);
            combined_Traces(z,2) = Trace(1,2);
            combined_TAD_id(z,2) = TAD_id(2);
            combined_Traces(z,3) = Trace(1,3);
            combined_TAD_id(z,3) = TAD_id(3);
            
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
            combined_Traces(z,1) = Trace(1,1);
            combined_TAD_id(z,1) = TAD_id(1);
            combined_Traces(z,2) = Trace(1,2);
            combined_TAD_id(z,2) = TAD_id(2);
            combined_Traces(z,3) = Trace(1,3);
            combined_TAD_id(z,3) = TAD_id(3);
            combined_Traces(z,4) = Trace(1,4);
            combined_TAD_id(z,4) = TAD_id(4);
        end
    end    
end

title('traces on segmentation stain 561 - per Terr')
axis equal
hold off




%%
%segment based on strain probe 561 HI

if primarytosegment_561 == 1
    FileName = ['sequential/560_denoised_strain_' SampleNum];
    load(FileName)
   
end



for m = 1:size(ImageStack560,3)
        ImageStack560(:,:,m) = imtransform(ImageStack560(:,:,m), tform2, 'XData', [1 1608], 'Ydata', [1 1608]);
end


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
BW = imbinarize(Icbrobr, TerrbinaryThresh560);
BW = imfill(BW, 'holes');
BW = imclose(BW, se);
figure
imshowpair(max(Icbrobr,[],3),max(BW,[],3),'montage')
title('before and after 3Dbinary')


%dilate binarymask
BW1 = BW;

se = strel('sphere',TerrDil560);
dilatedBW = imdilate(BW1,se);


%3D dist transform
D = -bwdist(~dilatedBW); 

%set background to its own catchment basin - where BW is 0
D(~dilatedBW) = -Inf;



%3D watershed

D = imhmin(D,TerrMin560); %the height threshold for suppressing shallow minima
L = watershed(D);
figure
imshow(label2rgb(L(:,:,50),'jet','w'))
title('watershed slice 50')

%

L = imdilate(L,true(5));
L = L-1; 


%adjust for sample drift in strain-imaging round

Ls = zeros(size(L));

if round(DeltaZ560) > 0
    for i = round(DeltaZ560) + 1 : size(Ls,3)
        Ls(:,:,i - round(DeltaZ560)) = L(:,:,i);
    end
else
    for i = 1:size(Ls,3) + round(DeltaZ560)
        Ls(:,:,i - round(DeltaZ560)) = L(:,:,i);
    end
end


if nuc == 1
    Ls(find(~nucL))= 0;
end


L=Ls; 


L1 = L;
clear L 

Hi_ID = unique(L1(:));
idx = unique(L1(:)) >0;
Hi_ID = Hi_ID(idx);

numterr = length(Hi_ID)


 
L1max = max(L1,[],3);
L1rgb = label2rgb(L1max,'jet','w','shuffle');
figure
imshow(L1rgb)
title('territories 560')
axis equal

figure
imagesc(Im560)
colormap gray
hold on
himage = imshow(L1rgb);
himage.AlphaData = 0.3;
title('territory labels superimposed transparently on 560')



%%segment based on strain probe 647

if primarytosegment_647 == 1
    FileName = ['sequential/647_denoised_strain_' SampleNum];
    load(FileName)
   
end


for m = 1:size(ImageStack647,3)
        ImageStack647(:,:,m) = imtransform(ImageStack647(:,:,m), tform2, 'XData', [1 1608], 'Ydata', [1 1608]);
end


for m = 1:size(ImageStack647,3)
        ImageStack647(:,:,m) = imtransform(ImageStack647(:,:,m), tform, 'XData', [1 1608], 'Ydata', [1 1608]);
end



Im647 = max(ImageStack647,[],3);
figure
imagesc(Im647)
colormap gray
axis equal
title('z stack max image 647 strain ')
caxis([100 700])


I = ImageStack647(ceil(y1):ceil(y3),ceil(x1):ceil(x2),:);

%correct for uneven background illumination
Im647 = max(I,[],3);
figure
imagesc(Im647)
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
title('selected FOV primary probe 647 bkg')

I2 = I - background;
figure
imagesc(max(I2,[],3))
colormap gray
axis equal
caxis([0 500])
title('selected FOV primary probe 647 minus bkg')

%set min and max values to even out territory intensity

I3 = I2;
I3(I3< zeroThre647) = 0;
I3(I3>Thresh647) = Thresh647;
figure
imagesc(max(I3,[],3))
colormap gray
axis equal
caxis([0 400])
title('selected FOV primary probe 647 thresh')



se = strel('disk', erodeThresh);

Icbre = imerode(I3, se);

Icbrobr = imreconstruct(Icbre, I3);
figure
imagesc(max(Icbrobr,[],3))
colormap gray
axis equal
title('closing and opening by reconstruction');


%3D binary
BW = imbinarize(Icbrobr, TerrbinaryThresh647);
BW = imfill(BW, 'holes');
BW = imclose(BW, se);
figure
imshowpair(max(Icbrobr,[],3),max(BW,[],3),'montage')
title('before and after 3Dbinary')


%
BW1 = BW;

se = strel('sphere',TerrDil647);
dilatedBW = imdilate(BW1,se);


%3D dist transform
D = -bwdist(~dilatedBW); 

%set background to its own catchment basin - where BW is 0
D(~dilatedBW) = -Inf;

%3D watershed

D = imhmin(D,TerrMin647); %the height threshold for suppressing shallow minima
L = watershed(D);
figure
imshow(label2rgb(L(:,:,50),'jet','w'))
title('watershed slice 50')






L = imdilate(L,true(5));
L = L-1;


Ls = zeros(size(L));

if round(DeltaZ647) > 0
    for i = round(DeltaZ647) + 1 : size(Ls,3)
        Ls(:,:,i - round(DeltaZ647)) = L(:,:,i);
    end
else
    for i = 1:size(Ls,3) + round(DeltaZ647)
        Ls(:,:,i - round(DeltaZ647)) = L(:,:,i);
    end
end



if nuc == 1
    Ls(find(~nucL))=0;
end


L=Ls; 


L2 = L;
clear L 

N2_ID = unique(L2(:));
idx = unique(L2(:)) >0;
N2_ID = N2_ID(idx);

numterr = length(N2_ID)
 

L2max = max(L2,[],3);
L2rgb = label2rgb(L2max,'jet','w','shuffle');
figure
imshow(L2rgb)
title('territories 647')
axis equal

figure
imagesc(Im647)
colormap gray
hold on
himage = imshow(L2rgb);
himage.AlphaData = 0.3;
title('territory labels superimposed transparently on 647')


%% remove overlapping territories from L1 (HI) and L2 (N2) segmentation


overlab = zeros(size(L1));


for i = 1:size(overlab,1)
    for j = 1:size(overlab,2)
        for k = 1:size(overlab,3)
        overlab(i,j,k) = L1(i,j,k) * L2(i,j,k);
        end
    end
end

overlab_terr = overlab; 
overlab = overlab + 1; 


overlab(find(overlab ~= 1))= 0; 


Hi_only = zeros(size(L1));


for i = 1:size(Hi_only,1)
    for j = 1:size(Hi_only,2)
        for k = 1:size(Hi_only,3)
        Hi_only(i,j,k) = L1(i,j,k) * overlab(i,j,k);
        end
    end
end

L1 = Hi_only; 


N2_only = zeros(size(L2));


for i = 1:size(N2_only,1)
    for j = 1:size(N2_only,2)
        for k = 1:size(N2_only,3)
        N2_only(i,j,k) = L2(i,j,k) * overlab(i,j,k);
        end
    end
end

L2 = N2_only; 



%%% 


nstats = regionprops3(nucL, 'all');
nstats( any(ismissing(nstats),2), :) = [];
ncenters = nstats.Centroid;


Terrstats = regionprops3(L4, 'all');
Terrstats( any(ismissing(Terrstats),2), :) = [];
Terrcenters = Terrstats.Centroid;



N2stats = regionprops3(L2, 'all');
N2stats( any(ismissing(N2stats),2), :) = [];
N2centers = N2stats.Centroid;


HIstats = regionprops3(L1, 'all');
HIstats( any(ismissing(HIstats),2), :) = [];
HIcenters = HIstats.Centroid;


%asign each primary territory (unique L4) to a nucleus 

for i =  1:length(Terrcenters)
    for j = 1:size(ncenters,1)
        dist2(i,j) = sqrt(sum((Terrcenters(i,:) - ncenters(j,:)) .^ 2, 2));
    end
end
    

%find the smallest distance and note down the ncenter of the territory
% Terrcellcenter contains all nuclei centers of each Territory

Terrcellcenter = [];

for i = 1:length(Terrcenters)   
    a = ncenters((dist2(i,:) == min(dist2(i,:)))',:);
    Terrcellcenter = [Terrcellcenter; a] ;
end



%find which nucL volume label Terrcellcenter is in and write down the
%number of nucleus in nucnum

%%
nucnums = zeros(length(Terrcenters),1);  %nucnum contains the label for nucL that the territory is in
Terrnum = zeros(size(ncenters,1),6) ;   %Territory ID each nucleus has inside; each row is a nucleus each column contains an Territory ID
for i = 1:size(ncenters,1)
    x = nuc_ID(i);
    for j = 1:length(Terrcellcenter)
        
        if Terrcellcenter(j,:)  == ncenters(i,:)
            nucnums(j) = x;
            for ii = 1: size(Terrnum,2)
                if Terrnum(i,ii) == 0
                    Terrnum(i,ii) = j;
                    break;
                end
            end
    
        end
    end
end


%%
% visulaization of all traces per nucleus


 n= 0;      
 tracecenters = zeros(1,3);
 trace_Volumes = cell(size(combined_Traces,1),size(combined_Traces,2)) ;
 

for i = 1 :size(Terrnum,1)
    nucID = nuc_ID(i);
    nucL1 = nucL;
    nucL1(find(nucL1 ~= nucID))= 0; 
    terr = Terrnum(i,:);
    idx = terr > 0;
    terr = terr(idx);
    L7 = zeros(size(L4));
    
    figure
    imagesc(ImDAPI)
    colormap gray
    hold on


    h1 = patch(isosurface(L1,0.5));
    h1.EdgeColor = 'none';
    h1.FaceColor = [1,0,0];
    h1.FaceAlpha = 0.15;

    h2 = patch(isosurface(L2,0.5));
    h2.EdgeColor = 'none';
    h2.FaceColor = [0,1,0];
    h2.FaceAlpha = 0.15; 
    
    h3 = patch(isosurface(nucL1,0.5));
    h3.EdgeColor = 'none';
    h3.FaceColor = [0,1,1];
    h3.FaceAlpha = 0.15;
    
    
    
    
    
    for iii = 1: length(terr)
        a = terr(iii);
        L4_temp = double(L4);
        terrID = terr_ID(a);
        L4_temp(find(L4_temp ~= terrID)) = 0;
        L1_tmp = L1;
        L1_tmp(find(~nucL1)) = 0;
         L2_tmp = L2;
        L2_tmp(find(~nucL1)) = 0;
        
    
        for ii = 1: size(combined_Traces,2)
            if ~isempty(combined_Traces{a,ii}) & length(combined_Traces{a,ii}) > 3
                
            n = n + 1
            

            interiorPoints = combined_Traces{a,ii};      %# Generate 200 3-D points
            DT = DelaunayTri(interiorPoints);  %# Create the tetrahedral mesh
            xyz = interiorPoints;
            tri = DT.Triangulation;
            FV = struct('faces',tri,'vertices',xyz);

            Volume = [size(I,1) size(I,2) size(I,3)];

            L5 = polygon2voxel(FV,Volume,'none');
            L5 = double(L5);
            L5(find(L5 ~= 0)) = n;
            
            trace_all{n,8} = L5;
            trace_all{n,9} = L1_tmp;
            trace_all{n,10} = L2_tmp;
            
            trace_stats = regionprops3(L5, 'all');
            trace_stats(any(ismissing(trace_stats),2), :) = [];
           
            tracecenters(n,:) = trace_stats.Centroid;            %% center point of each trace > 3 TADs
            
            trace_all(n,1) = num2cell(a);
            trace_all(n,2) = num2cell(ii);
             trace_all(n,3) = num2cell(nucID);
             trace_all{n,4} = combined_Traces{a,ii};
             trace_all{n,5} = combined_TAD_id{a,ii};
             
              
                
                    if ii == 1
                    scatter3(combined_Traces{a,ii}(:,1), combined_Traces{a,ii}(:,2), combined_Traces{a,ii}(:,3),'or', 'MarkerFaceColor', 'r')
                    plot3(combined_Traces{a,ii}(:,1), combined_Traces{a,ii}(:,2), combined_Traces{a,ii}(:,3),'-r');
                    text(combined_Traces{a,ii}(:,1), combined_Traces{a,ii}(:,2),combined_Traces{a,ii}(:,3), num2str(combined_TAD_id{a,ii}), 'Color', 'Black', ...
                        'FontWeight', 'Bold', 'FontSize',8);
                    plot3(tracecenters(n,1), tracecenters(n,2), tracecenters(n,3));
                    text(tracecenters(n,1), tracecenters(n,2),tracecenters(n,3), num2str(n), 'Color', 'Black', ...
                        'FontWeight', 'Bold', 'FontSize',16);
                    
                    elseif ii ==2

                    scatter3(combined_Traces{a,ii}(:,1), combined_Traces{a,ii}(:,2), combined_Traces{a,ii}(:,3),'sb', 'MarkerFaceColor', 'b')
                    plot3(combined_Traces{a,ii}(:,1), combined_Traces{a,ii}(:,2), combined_Traces{a,ii}(:,3),'-b');
                    text(combined_Traces{a,ii}(:,1), combined_Traces{a,ii}(:,2),combined_Traces{a,ii}(:,3), num2str(combined_TAD_id{a,ii}), 'Color', 'Black', ...
                        'FontWeight', 'Bold', 'FontSize',8);
                     plot3(tracecenters(n,1), tracecenters(n,2), tracecenters(n,3));
                    text(tracecenters(n,1), tracecenters(n,2),tracecenters(n,3), num2str(n), 'Color', 'Black', ...
                        'FontWeight', 'Bold', 'FontSize',16);
                    elseif ii == 3

                        scatter3(combined_Traces{a,ii}(:,1), combined_Traces{a,ii}(:,2), combined_Traces{a,ii}(:,3),'dg', 'MarkerFaceColor', 'g')
                    plot3(combined_Traces{a,ii}(:,1), combined_Traces{a,ii}(:,2), combined_Traces{a,ii}(:,3),'-g');
                    text(combined_Traces{a,ii}(:,1), combined_Traces{a,ii}(:,2),combined_Traces{a,ii}(:,3), num2str(combined_TAD_id{a,ii}), 'Color', 'Black', ...
                        'FontWeight', 'Bold', 'FontSize',8);
                     plot3(tracecenters(n,1), tracecenters(n,2), tracecenters(n,3));
                    text(tracecenters(n,1), tracecenters(n,2),tracecenters(n,3), num2str(n), 'Color', 'Black', ...
                        'FontWeight', 'Bold', 'FontSize',16);

                    elseif ii == 4
                        scatter3(combined_Traces{a,ii}(:,1), combined_Traces{a,ii}(:,2), combined_Traces{a,ii}(:,3),'^k', 'MarkerFaceColor', 'k')
                    plot3(combined_Traces{a,ii}(:,1), combined_Traces{a,ii}(:,2), combined_Traces{a,ii}(:,3),'-k');
                    text(combined_Traces{a,ii}(:,1), combined_Traces{a,ii}(:,2),combined_Traces{a,ii}(:,3), num2str(combined_TAD_id{a,ii}), 'Color', 'Black', ...
                        'FontWeight', 'Bold', 'FontSize',8);
                     plot3(tracecenters(n,1), tracecenters(n,2), tracecenters(n,3));
                    text(tracecenters(n,1), tracecenters(n,2),tracecenters(n,3), num2str(n), 'Color', 'Black', ...
                        'FontWeight', 'Bold', 'FontSize',16);

                    end
             end
        end
    end
    
  
    legend('HI','N2', 'nucL');
    axis equal
    set(gca,'YDir','reverse')
    title(['traces in segmentation nucleus' num2str(nucID)])
    
 end
 
 
 
%% is trace overlapping into other volume - assigning identity to all Region Foci

for i = 1:size(trace_all,1)
    
    for ii = 1: length(trace_all{i,4}(:,1))
       if L1(round(trace_all{i,4}(ii,2)), round(trace_all{i,4}(ii,1)), round(trace_all{i,4}(ii,3))) > 0 & L2(round(trace_all{i,4}(ii,2)),round(trace_all{i,4}(ii,1)), round(trace_all{i,4}(ii,3))) == 0
            trace_all{i,6}(ii,1) =  1;   % HI only TAD
       
       elseif L1(round(trace_all{i,4}(ii,2)), round(trace_all{i,4}(ii,1)), round(trace_all{i,4}(ii,3))) == 0 & L2(round(trace_all{i,4}(ii,2)),round(trace_all{i,4}(ii,1)), round(trace_all{i,4}(ii,3))) > 0
            trace_all{i,6}(ii,1) =  2;
            
       else 
            trace_all{i,6}(ii,1) =  0;
            if trace_all{i,6}(ii,1) ==  0;
                L1_tmp = trace_all{i,9};
                L2_tmp = trace_all{i,10};
                ID1 = unique(L1_tmp(:));
                idx = unique(L1_tmp(:)) >0;
                ID1 = ID1(idx);
                
                ID2 = unique(L2_tmp(:));
                idx = unique(L2_tmp(:)) >0;
                ID2 = ID2(idx);

                
                boundaries1 = cell(size(L1_tmp,3),1);
                boundaries2 = cell(size(L2_tmp,3),1);
                for x = 1:size(L1_tmp,3)

                    boundaries1{x,1} = bwboundaries(L1_tmp(:,:,x));
                end

                for x = 1:size(L2_tmp,3)

                boundaries2{x,1} = bwboundaries(L2_tmp(:,:,x));
                end


                boundarie_points1 = [];


                for iii = 1: size(boundaries1,1)
                    if ~isempty(boundaries1{iii,1})

                        boundarie_points1_y = boundaries1{iii,1}{1}(:,1);
                        boundarie_points1_x = boundaries1{iii,1}{1}(:,2);
                        boundarie_points1_z = zeros(length(boundaries1{iii,1}{1}(:,2)),1) + iii ;
                        boundarie_points1 = [boundarie_points1; [boundarie_points1_x,  boundarie_points1_y,  boundarie_points1_z]];
                    end
                end
                
                if length(ID1) > 0
                    dist5 = sqrt(sum((boundarie_points1 - trace_all{i,4}(ii,:)) .^ 2, 2));
                end
                boundarie_points2 = [];


                for iiii = 1: size(boundaries2,1)
                    if ~isempty(boundaries2{iiii,1})

                        boundarie_points2_y = boundaries2{iiii,1}{1}(:,1);
                        boundarie_points2_x = boundaries2{iiii,1}{1}(:,2);
                        boundarie_points2_z = zeros(length(boundaries2{iiii,1}{1}(:,2)),1) + iiii ;
                        boundarie_points2 = [boundarie_points2; [boundarie_points2_x,  boundarie_points2_y,  boundarie_points2_z]];
                    end
                end
                
                 if length(ID2) > 0
                    dist6 = sqrt(sum((boundarie_points2 - trace_all{i,4}(ii,:)) .^ 2, 2));
                 end
                 
                 if length(ID1) > 0 & length(ID2) > 0
                 
                     if min(dist6) < min(dist5) 
                         trace_all{i,6}(ii,1) = 20;
                     elseif min(dist6) > min(dist5)           
                          trace_all{i,6}(ii,1) = 10;
                     end
                 end
                
            end

       end
       if overlab_terr(round(trace_all{i,4}(ii,2)), round(trace_all{i,4}(ii,1)), round(trace_all{i,4}(ii,3))) > 0 
            trace_all{i,11}(ii,1) =  1;
       
       else trace_all{i,11}(ii,1) =  0;
       end
    end
end




%% second categorizer using whats the biggest contributor to the trace - HI-Foci or N2-Foci

for ii = 1:size(trace_all,1)
    
    if length(find(trace_all{ii,6} ==1)) + length(find(trace_all{ii,6} ==10))  > length(find(trace_all{ii,6} == 2)) + length(find(trace_all{ii,6} == 20))
        trace_all(ii,7) = num2cell(1);
    elseif length(find(trace_all{ii,6} ==1)) + length(find(trace_all{ii,6} ==10)) < length(find(trace_all{ii,6} ==2)) + length(find(trace_all{ii,6} == 20))
        trace_all(ii,7) = num2cell(2);
    else
        trace_all(ii,7) = num2cell(3);
    end
end

    



%%                 saving based on nucleus
 
 %%% each result file coresponds to one nucleus
 %%% for Trace_N2, Trace_HI: each row is a different territory within the
 %%% nucleus, each column is a trace (can be up to 4); if there was no
 %%% trace for HI or N2 detected within the territory it will be left
 %%% blanck
 %%% 
 
 
idx =  unique([trace_all{:,3}]);
        
    
for i = 1: length(idx)
    nuc = idx(i);
    tmp = find([trace_all{:,3}] == nuc);
    idx2 = unique([trace_all{tmp}]);     %number of unique terr in Nuc 
    u = 0;                       %defines line of HI Traces
    
    Trace_N2 = cell(length(idx2),3);
    Trace_HI = cell(length(idx2),3);
    Trace_both = cell(length(idx2),3);
    TAD_id_N2 = cell(length(idx2),3);
    TAD_id_HI = cell(length(idx2),3);
    TAD_id_both = cell(length(idx2),3);
    strain_ID_N2 = cell(length(idx2),3);
    strain_ID_HI = cell(length(idx2),3);
    strain_ID_both = cell(length(idx2),3);
    overlab_id_N2 = cell(length(idx2),3);
    overlab_id_HI = cell(length(idx2),3);
    overlab_id_both = cell(length(idx2),3);
    
    terr_id = cell(length(idx2),1);
    
    for ii = 1:length(idx2)
        terr = idx2(ii);
        tmp2 = find([trace_all{:,1}] == terr);
        
        t = 0;      %counts trace - column
        
        u = u + 1; %counts territory - row
        
        terr_id{u,1} = terr; 
        for iii = 1:length(tmp2)
            trc = tmp2(iii);
        
            if trace_all{trc,7} == 1
                
                t = t +1;
                Trace_HI(u,t)  = trace_all(trc,4);
                TAD_id_HI(u,t) = trace_all(trc,5);
                strain_ID_HI(u,t) = trace_all(trc,6);
                overlab_id_HI(u,t) = trace_all(trc,11);
               
            elseif  trace_all{trc,7} == 2
                
                t = t +1;
                Trace_N2(u,t) = trace_all(trc,4);
                TAD_id_N2(u,t) = trace_all(trc,5);
                strain_ID_N2(u,t) = trace_all(trc,6);
                overlab_id_N2(u,t) = trace_all(trc,11);
                
                elseif  trace_all{trc,7} == 3
                    t = t +1;
                    Trace_both(u,t) = trace_all(trc,4);
                    TAD_id_both(u,t) = trace_all(trc,5);
                    strain_ID_both(u,t) = trace_all(trc,6);
                    overlab_id_both(u,t) = trace_all(trc,11);
                
            end
        end
    end
    
    if save_traces == 1 && save_age == 0
        save(['traces/' num2str(man_age) '_cell_Traces_Sample' SampleNum '_FOV' num2str(FOV) '_roi561_' num2str(i) '.mat'],...
            'Trace_N2','Trace_HI','Trace_both','TAD_id_HI','TAD_id_N2','TAD_id_both','strain_ID_N2','strain_ID_HI','strain_ID_both', 'overlab_id_N2','overlab_id_HI','overlab_id_both','L4', 'nucL','terr_id', 'L1', 'L2', 'nuc', 'overlab_terr', 'TADsToExclude', 'x1', 'x2', 'y1', 'y3');
    end
    if save_traces == 1 && save_age ==1
        save(['traces/' num2str(age) '_cell_Traces_Sample' SampleNum '_FOV' num2str(FOV) '_roi561_' num2str(i) '.mat'],...
            'Trace_N2','Trace_HI','Trace_both', 'TAD_id_HI','TAD_id_N2','TAD_id_both','strain_ID_N2','strain_ID_HI','strain_ID_both', 'overlab_id_N2','overlab_id_HI','overlab_id_both', 'L4', 'nucL','terr_id','L1', 'L2', 'nuc', 'overlab_terr', 'TADsToExclude', 'x1', 'x2', 'y1', 'y3');
    end
end
    


       