%% this script extracts all traces from a given Chr.mat file with exactly two traces per nucleus
% and calculates the distance between homologs within the same nucleus
%%can be adapted to other questions by changing the criteria which traces
%%are taken into account (for example: all nuclei with less then 4 traces)
%distance matrixes will be generated at the end
%as a last step, you can visualize the intra-chromosomal distances of the
%filtered homologs and generate a distance matrix 


clear all 
close all 

AgeRange = '2to40cell';     %needed for labeling of your figures later

TotalTADNum = 22;           %how many regions do your traces have 
load('2to40cell/AllChromosomes561_fused_singles.mat')


%% step1 store nucID from Traces in Chr variable = in which nucleus is my trace

for ii = 1 :length(Chr)
    fileName = Chr(ii).TraceName;
    load(fileName)
    
    x = round(Trace{1,1}(1,1));
    y = round(Trace{1,1}(1,2));
    z = round(Trace{1,1}(1,3));
    
    Chr(ii).nucID = nucL(y,x,z);
end

%% extract Dataset, Sample and FOV for all traces and store it in the Chr.mat - this is needed to filter later

for ii = 1 :length(Chr)
    fileName = Chr(ii).TraceName;
    SampleIndex = strfind(fileName,'Sample');   %strfind finds the startindex of your querry within a string
    
    FOVindex = strfind(fileName,'FOV');
    
    RoiIndex = strfind(fileName,'roi');
    
    
    
    if length(SampleIndex) ~= 1
        return
    end
    
    Chr(ii).Sample = str2num(fileName(SampleIndex + 6 : FOVindex - 2));
    Chr(ii).Fov = str2num(fileName(FOVindex + 3 : RoiIndex - 2));
    
    Chr(ii).DatesetID = str2num(fileName(8:15));    %check if your datastructure fits to this - important is to extract a unique identifier for your dataset
     
end

%% this section is now extracting all traces which are full-filling the criterium that there is exactly 2 traces per nucleus
% first: all traces from 1 datasets are selected - second: all traces from
% 1 sample are selected - third: all traces from one FOV are selected -
% fourth: all traces from one nucleus are selected

datasets = unique([Chr.DatesetID]);
h = 0;

for ii = 1 :length(datasets)
    n = 0;
    for i = 1 : length(Chr)
        if Chr(i).DatesetID == datasets(ii)
            n = n + 1;
            Chr_Dataset(n) = Chr(i);
        end
 
    end
    
    Samples = unique([Chr_Dataset.Sample]);
    for a = 1 : length(Samples)
        n = 0;
        for i = 1 : length(Chr_Dataset)
            if Chr_Dataset(i).Sample == Samples(a)
                n = n + 1;
                Chr_Sample(n) = Chr_Dataset(i);
            end
        end
        
        FOV = unique([Chr_Sample.Fov]);
        for b = 1 : length(FOV)
            n = 0;
            for i = 1 : length(Chr_Sample)
                if Chr_Sample(i).Fov == FOV(b)
                    n = n + 1;
                    Chr_FOV(n) = Chr_Sample(i);
                end
            end
            
            nuc = unique([Chr_FOV.nucID]);
            for c = 1: length(nuc) 
                n = 0;
                for i = 1: length(Chr_FOV)
                    if Chr_FOV(i).nucID == nuc(c)
                        n = n + 1;
                        Chr_nuc(n) = Chr_FOV(i);
                    end
       
                end
                if length(Chr_nuc)== 2     %this selects all traces with 2 traces per nucleus and stores it in pseudo-Chr 
                    h = h+1;               %each line of Chr_h1 & Chr_h2 represents a pair of homologs for 1 nucleus
                    
                    Chr_h1(h) = Chr_nuc(1);
                    Chr_h2(h) = Chr_nuc(2);
                end
                clear Chr_nuc
            end
            clear Chr_FOV
        end
        clear Chr_Sample
    end
    clear Chr_Dataset
end
                
save('2to40cell/homologs/AllChromosomes561_homologs.mat', 'Chr_h1', 'Chr_h2')




%% this section measures the distances between homologs, but gives the trace a pseudo-directionality

for i = 1:TotalTADNum
    for j = 1:TotalTADNum
        DisList_hom = [];
        for k = 1:length(Chr_h1) 
            if Chr_h1(k).r(i) == 1 && Chr_h2(k).r(j) == 1  
                
                    DisList_hom = [DisList_hom ((Chr_h1(k).x(i)-Chr_h2(k).x(j))^2+(Chr_h1(k).y(i)-Chr_h2(k).y(j))^2+(Chr_h1(k).z(i)-Chr_h2(k).z(j))^2)^0.5];
                
            end
        end  
        DisListAll_hom{i,j} = DisList_hom;
    end
end

%% this section accounts for the fact that distances H1-R1 to H2-R2 = H1-R2 to H2-R1
for i = 1: TotalTADNum
    for j = 1: TotalTADNum
        if i ~= j & i < j 
            a = DisListAll_hom{i,j}(:);
            b = DisListAll_hom{j,i}(:);
            c = [a ; b]';
            DisListAll_hom{i,j} = c; 
            DisListAll_hom{j,i} = c;
            Mean_hom(i,j) = mean(DisListAll_hom{i,j}(:));
            Mean_hom(j,i) = mean(DisListAll_hom{j,i}(:));
        end
        if i == j 
            Mean_hom(i,j) = mean(DisListAll_hom{i,j}(:));
            Mean_hom(j,i) = mean(DisListAll_hom{j,i}(:));
        end
    end
end

% this section excludes region 11 (remove or adjust for your data)


Mean_hom(11,:) = [];
Mean_hom(:,11) = [];

TotalTADNum = 21; 



%% mean figure for homolog comparison


figure1 = figure;
colormap(hot);
% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');
% Create image
image(Mean_hom,'Parent',axes1,'CDataMapping','scaled');
title(['Mean H1 and H2 pairwise distance (um) -'  AgeRange]);
% Uncomment the following line to preserve the X-limits of the axes
xlim(axes1,[0.5 TotalTADNum + 0.5]);
% Uncomment the following line to preserve the Y-limits of the axes
ylim(axes1,[0.5 TotalTADNum + 0.5]);
box(axes1,'on');
axis(axes1,'ij');
colorbar('peer',axes1);
caxis([0 4])
PlotProp
axis square
%print('Mean_early','-dsvg')


%% this section compares distance distribution of homologous regions - plot can be adapted to your liking 
TotalTADNum = 22;

edges = linspace(0, 6, 23); %defines  
for a = 1:TotalTADNum


figure(1000)
subplot(5,5,a)
histogram(DisListAll_hom{a,a},'BinEdges',edges)
grid on;
xlim([0, 6]);
title(['Region' num2str(a)])
end
saveas(figure(1000),'distances_distribution_homologRegions' ,'fig');



for r = 1:TotalTADNum
    for a = 1:TotalTADNum
    figure(1000+r)
    subplot(5,5,a)
    histogram(DisListAll_hom{r,a},'BinEdges',edges)
    grid on;
    xlim([0, 6]);
    title(['Region' num2str(a)])
    end
  sgtitle(['Region' num2str(r)]) 
  saveas(figure(1000+r),(['distances_distribution_Region' num2str(r)]) ,'fig');

end







%% combine Chr_h1 and Chr_h2 to visualize distance matrix (intra) for all filtered chromosomes 

Chr = [Chr_h1 Chr_h2];


%%
Mean = zeros(TotalTADNum,TotalTADNum);
Std = zeros(TotalTADNum,TotalTADNum);
Median = zeros(TotalTADNum,TotalTADNum);


for i = 1:TotalTADNum
    for j = 1:TotalTADNum
        DisList = [];
        for k = 1:length(Chr) 
            if Chr(k).r(i) == 1 && Chr(k).r(j) == 1
                
                    DisList = [DisList ((Chr(k).x(i)-Chr(k).x(j))^2+(Chr(k).y(i)-Chr(k).y(j))^2+(Chr(k).z(i)-Chr(k).z(j))^2)^0.5];
                
            end
        end
        Mean(i,j) = mean(DisList);
        Std(i,j) = std(DisList);
        SEM(i,j) = std(DisList)/(length(DisList))^0.5;
        NofData(i,j) = length(DisList);
        DisListAll{i,j} = DisList;
        Median(i,j) = median(DisList);
    end
end



%% make figures

figure2 = figure;
colormap(hot);
% Create axes
axes1 = axes('Parent',figure2);
hold(axes1,'on');
% Create image
image(Mean,'Parent',axes1,'CDataMapping','scaled');
title(['Mean pairwise distance (um) -'  AgeRange]);
% Uncomment the following line to preserve the X-limits of the axes
xlim(axes1,[0.5 TotalTADNum+0.5]);
% Uncomment the following line to preserve the Y-limits of the axes
ylim(axes1,[0.5 TotalTADNum+0.5]);
box(axes1,'on');
axis(axes1,'ij');
colorbar('peer',axes1);
[r,c] = find(isnan(Mean)) ;
hold on ; 
plot(c, r, 'yo') ;
caxis([0 2])
PlotProp
axis square
%print('Mean_early','-dsvg')

figure3 = figure;
colormap(hot);
% Create axes
axes1 = axes('Parent',figure3);
hold(axes1,'on');
% Create image
image(Median,'Parent',axes1,'CDataMapping','scaled');
title(['Median pairwise distance (um) -'  AgeRange]);
% Uncomment the following line to preserve the X-limits of the axes
xlim(axes1,[0.5 TotalTADNum+0.5]);
% Uncomment the following line to preserve the Y-limits of the axes
ylim(axes1,[0.5 TotalTADNum+0.5]);
box(axes1,'on');
axis(axes1,'ij');
colorbar('peer',axes1);
%caxis([0 2.5])
PlotProp
axis square
%print('Median_early','-dsvg')

figure
imagesc(SEM)
colorbar
title(['SEM of pairwise distance - ' AgeRange]);
PlotProp
axis square
%print('SEM_early','-dsvg')

figure
imagesc(Std)
colorbar
title(['Std of pairwise distance - ' AgeRange]);
PlotProp
axis square
%print('STD_early','-dsvg')

figure
imagesc(NofData)
colorbar
title(['Number of measurements - ' AgeRange]);
PlotProp
axis square
%print('NofDATA_early','-dsvg')

histogram([Chr.TracesInTerritory])
title('wild-type segmentation results')
ylabel('number of traces')
xlabel('total traces in territory volume')
PlotProp

TracesInTerritory = [Chr.TracesInTerritory].';

TracesPerSegment = [];
TracesPerSegment(1,1) = sum(TracesInTerritory(:) == 1)/length(Chr)*100;
TracesPerSegment(1,2) = sum(TracesInTerritory(:) == 2)/length(Chr)*100;
TracesPerSegment(1,3) = sum(TracesInTerritory(:) == 3)/length(Chr)*100;
TracesPerSegment(1,4) = sum(TracesInTerritory(:) == 4)/length(Chr)*100;
figure
bar(TracesPerSegment)
title('wild-type segmentation results')
ylabel('% of total traces')
xlabel('traces within territory volume')
ylim([0 100])
PlotProp
%print('wt_segmentation_early','-dsvg')

%% save data

save([AgeRange '\homologs\MeanPairDistanceMatrix_all' num2str(SegChannel) '.mat'], 'Mean');
save([AgeRange '\homologs\MedianPairDistanceMatrix_all' num2str(SegChannel) '.mat'], 'Median');
save([AgeRange '\homologs\StdPairDistanceMatrix_all' num2str(SegChannel) '.mat'], 'Std');
save([AgeRange '\homologs\SEMPairDistanceMatrix_all' num2str(SegChannel) '.mat'], 'SEM');
save([AgeRange '\homologs\NofDataPairDistanceMatrix_all' num2str(SegChannel) '.mat'], 'NofData');
save([AgeRange '\homologs\AllChromosomes_all' num2str(SegChannel) '.mat'], 'Chr','-v7.3');
save([AgeRange '\homologs\DisListAll_all' num2str(SegChannel) '.mat'], 'DisListAll');

save([AgeRange '\homologs\DisListAll_homologs' num2str(SegChannel) '.mat'], 'DisListAll_hom');
save([AgeRange '\homologs\MeanPairDistanceMatrix_homologs' num2str(SegChannel) '.mat'], 'Mean_hom');



       