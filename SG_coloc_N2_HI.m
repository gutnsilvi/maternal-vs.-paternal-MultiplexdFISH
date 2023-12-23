% This script uses stringent segmentation to calculate the overlap between
% two different chromosome territory signals 
% Version 1.0
% 
% Silvia Gutnik, V1, 05.07.2022
% % % % % % % % % % % % % % % % % % 
tic

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
FOV = 1; %which embryo are we analyzing now (if > 1 in the sample)

div = 30; 

FileName = ['tform_strain_' num2str(SampleNum) '.mat']; 
load(FileName) 
FileName = ['DeltaZ_strain_' num2str(SampleNum) '.mat'];
load(FileName)
load('tform.mat');

%threshold fraction for nuclear size exclusion (p-body-small speckles)


save_quant = 0;   % set to 1 if you want to save the result

save_age = 1;
man_age = 0;


%%%%%%%%%%maxintensitythreshs%%%%%%%%%%%%%%%%%%%%%%change 405 mainly%%%%%%%
Thresh405 = 3000;
Thresh560 = 70;

Thresh647 = 70;
zeroThre560 = 40;        %strain-label HI


zeroThre647 = 30;        %strain-label N2
zeroThre405 = 700;      %nuclear stain

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%for reconstruction
dilateThresh = 5; % 
erodeThresh = 5; % 

%for binarization and watershed strength%%%change terrmin if traces are too
%many/few per territory

TerrbinaryThresh560 = 3;
TerrbinaryThresh647 = 3;
TerrMin560 = 0.45;
TerrMin647 = 0.45;
NucbinaryThresh = 500; 
NucMin= 0.30; 


% % % find and segment nuclei

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

I_DAPI = ImageStack405(ceil(y1):ceil(y3),ceil(x1):ceil(x2),:);
figure
ImDAPI = max(I_DAPI,[],3);
imagesc(ImDAPI)
colormap gray
axis equal
caxis([0 10000])
title('selected FOV 405')

%correct for uneven background illumination

se = strel('disk',50);
background = imopen(I_DAPI,se);
figure
imagesc(max(background,[],3))
colormap gray
axis equal
caxis([0 3500])
title('selected FOV primary probe 405 bkg')

I2 = I_DAPI - background;
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



nstats = regionprops3(nucL, 'all');
nstats( any(ismissing(nstats),2), :) = [];
ncenters = nstats.Centroid;
nVolumes = nstats.Volume;
nuc_ID = unique(nucL(:));
idx = unique(nucL(:)) >0;
nuc_ID = nuc_ID(idx);

age = length(nuc_ID)

%


Lmax = max(nucL,[],3);
nucLrgb = label2rgb(Lmax,'jet','w','shuffle');
figure
imagesc(ImDAPI)
colormap gray


hold on
himage = imshow(nucLrgb);
himage.AlphaData = 0.3;
hold on

for i = 1: size(ncenters,1)

plot3(ncenters(i,1), ncenters(i,2), ncenters(i,3));
 text(ncenters(i,1), ncenters(i,2),ncenters(i,3), num2str(nuc_ID(i)), 'Color', 'Blue', ...
                'FontWeight', 'Bold', 'FontSize',16);
end
axis equal



WhetherROI = questdlg('Do you want to exclude nuclei?'); %ask question
if(strcmp(WhetherROI, 'Yes'))
    answer = inputdlg('Enter nuclei ID space-separated :',...
             'Nuclei Id', [1 50]);
    nuc_ex = str2num(answer{1});
    
end

if exist('nuc_ex')

    for i = 1: length(nuc_ex) 
        nucL(find(nucL == nuc_ex(i)))= 0; 
    end
end
    
    
nstats = regionprops3(nucL, 'all');
nstats( any(ismissing(nstats),2), :) = [];
ncenters = nstats.Centroid;
nVolumes = nstats.Volume;


nuc_ID = unique(nucL(:));
idx = unique(nucL(:)) >0;
nuc_ID = nuc_ID(idx);


Lmax = max(nucL,[],3);
nucLrgb = label2rgb(Lmax,'jet','w','shuffle');
figure
imagesc(ImDAPI)
colormap gray
hold on
himage = imshow(nucLrgb);
himage.AlphaData = 0.3;
hold on

for i = 1: size(ncenters,1)

plot3(ncenters(i,1), ncenters(i,2), ncenters(i,3));
 text(ncenters(i,1), ncenters(i,2),ncenters(i,3), num2str(nuc_ID(i)), 'Color', 'Blue', ...
                'FontWeight', 'Bold', 'FontSize',16);
end
axis equal 
title('nuclei labels superimposed transparently on DAPI - filtered')



%



%%


FileName = ['sequential/647_denoised_strain_' SampleNum  '.mat'];
load(FileName)

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

% set min and max values to even out territory intensity

I3 = I2;
I3(I3< zeroThre647) = 0;
I3(I3>Thresh647) = Thresh647;
figure
imagesc(max(I3,[],3))
colormap gray
axis equal
caxis([0 400])
title('selected FOV primary probe 647 thresh')


BW = imbinarize(I3, TerrbinaryThresh647);
figure
imshowpair(max(I,[],3),max(BW,[],3),'montage')
title('before and after 3Dbinary')


%3D dist transform
D = -bwdist(~BW); 

%set background to its own catchment basin - where BW is 0
D(~BW) = -Inf;

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

Ls(find(~nucL))=0;

L=Ls; 

L2 = L;
clear L 


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
himage = imshow(nucLrgb);
himage.AlphaData = 0.2;
himage = imshow(L2rgb);
himage.AlphaData = 0.3;
title('territory labels superimposed transparently on 647')


%%

FileName = ['sequential/560_denoised_strain_' SampleNum];
load(FileName)


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


%3D binary
BW = imbinarize(I3, TerrbinaryThresh560);
figure
imshowpair(max(I3,[],3),max(BW,[],3),'montage')
title('before and after 3Dbinary')


%3D dist transform
D = -bwdist(~BW); 

%set background to its own catchment basin - where BW is 0
D(~BW) = -Inf;



%3D watershed

D = imhmin(D,TerrMin560); %the height threshold for suppressing shallow minima
L = watershed(D);
figure
imshow(label2rgb(L(:,:,50),'jet','w'))
title('watershed slice 50')

%

L = imdilate(L,true(5));
L = L-1; 


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

Ls(find(~nucL))= 0;

L=Ls; 


L1 = L;
clear L 


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
himage = imshow(nucLrgb);
himage.AlphaData = 0.2;
himage = imshow(L1rgb);
himage.AlphaData = 0.3;
title('territory labels superimposed transparently on 560')


%
overlab = zeros(size(L1));


for i = 1:size(overlab,1)
    for j = 1:size(overlab,2)
        for k = 1:size(overlab,3)
        overlab(i,j,k) = L1(i,j,k) * L2(i,j,k);
        end
    end
end

overlab_terr = overlab; 


%
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

h3 = patch(isosurface(overlab_terr ,0.5));
h3.EdgeColor = 'none';
h3.FaceColor = [0,0,1];
h3.FaceAlpha = 0.25;

legend('HI','N2', 'overlap');
axis equal
set(gca,'YDir','reverse')
title('overlap')


% overlap percent per nucleus

for i = 1:length(nuc_ID)
    nuc = nuc_ID(i);
    nucL1 = nucL;
    nucL1(find(nucL1 ~= nuc))= 0; 
    L1_temp = L1;
    L1_temp(find(~nucL1))= 0; 
    L2_temp = L2;
    L2_temp(find(~nucL1))= 0; 
    overlab_terr_temp = overlab_terr;
    overlab_terr_temp(find(~nucL1))= 0;

    HIstats = regionprops3(L1_temp, 'Volume');
    Vol(i).Hi = sum(HIstats.Volume);
    N2stats = regionprops3(L2_temp, 'Volume');
    Vol(i).N2 = sum(N2stats.Volume);
    OverlabStats = regionprops3(overlab_terr_temp, 'Volume');
    Vol(i).overlab = sum(OverlabStats.Volume);

    figure
    imagesc(ImDAPI)
    colormap gray
    hold on
    

    h1 = patch(isosurface(L1_temp,0.5));
    h1.EdgeColor = 'none';
    h1.FaceColor = [1,0,0];
    h1.FaceAlpha = 0.15;
    
    h2 = patch(isosurface(L2_temp,0.5));
    h2.EdgeColor = 'none';
    h2.FaceColor = [0,1,0];
    h2.FaceAlpha = 0.15; 
    
    h3 = patch(isosurface(overlab_terr_temp ,0.5));
    h3.EdgeColor = 'none';
    h3.FaceColor = [0,0,1];
    h3.FaceAlpha = 0.25;

legend('HI','N2', 'overlap');
axis equal
set(gca,'YDir','reverse')
title('overlap')


end

%%
if save_quant == 1;
 save(['overlap/' num2str(age) 'cell_Sample' SampleNum '_FOV' num2str(FOV) '.mat'],...
            'Vol','x1', 'x2', 'y1', 'y3');
end
