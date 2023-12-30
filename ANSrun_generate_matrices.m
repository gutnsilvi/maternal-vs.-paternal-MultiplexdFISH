%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script compiles all traces in a mat-file, filtered by age and
%%extracts mean, median, SEM, Std
%
%   Input: 
%   chromosome traces
%
%   Output: 
%   AllChromosomes file
%   mean pairwise distances
%   median pairwise distances
%   SEM pairwise distances
%   Std pairwise distances
%   Distance List of all pairwise distances of all chromosomes
%   Number of datapoints for all pairwise distances of all chromosomes
%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic;
clear all
close all

Folder{1} = 'traces1';
Folder{2} = 'traces2';
Folder{3} = 'traces3';
Folder{4} = 'traces4';
Folder{5} = 'traces5';
Folder{6} = 'traces6';


NumSamples = 16;
AgeStart = 2; %change these to specify embryo age
AgeStop = 40;
MaxROI = 200; %change to specify largest roi# 
MaxFOV = 7; %change to specify largest FOV# 

AgeRange = ([ num2str(AgeStart) 'to' num2str(AgeStop) 'cell']); 
SegChannel = 561;
sz = 1608;

path = [Folder '/' AgeRange];
mkdir(path)

TotalTADNum = 22;
n = 0;
flag = 0;

for ii = 1:length(Folder)
    for i = 0:NumSamples
        for q = AgeStart:AgeStop
            for b = 1:MaxFOV
                for j = 1:MaxROI
                    if exist([Folder '/' num2str(q) '_cell_Traces_Sample' num2str(i) '_FOV' num2str(b) '_roi' num2str(SegChannel) '_' num2str(j) '.mat'])
                        load([Folder '/' num2str(q) '_cell_Traces_Sample' num2str(i) '_FOV' num2str(b) '_roi' num2str(SegChannel) '_' num2str(j) '.mat']);
                        TraceName = num2str([Folder '/' num2str(q) '_cell_Traces_Sample' num2str(i) '_FOV' num2str(b) '_roi' num2str(SegChannel) '_' num2str(j) '.mat']);
                        flag = 1;                
                    end
                    if flag == 1
                        if length(Trace) ==4
                            for g = 1:4
                                if length(TAD_id{g})>1 && length(TAD_id{1,g}) == length(unique(TAD_id{1,g})) %skip trace if it has multiple foci per TAD/region (remove ambiguous traces)
                                    n = n+1;
                                    Chr(n).x = zeros(TotalTADNum,1);
                                    Chr(n).y = zeros(TotalTADNum,1);
                                    Chr(n).z = zeros(TotalTADNum,1);
                                    Chr(n).r = zeros(TotalTADNum,1);
                                    
                                    for k = 1: length(TAD_id{1,g})
                                        Chr(n).x(TAD_id{1,g}(k,1)) = Trace{1,g}(k,1)*110/1000; % convert to um from pixel units
                                        Chr(n).y(TAD_id{1,g}(k,1)) = sz*110/1000-Trace{1,g}(k,2)*110/1000;
                                        Chr(n).z(TAD_id{1,g}(k,1)) = Trace{1,g}(k,3)*0.2;
                                        Chr(n).r(TAD_id{1,g}(k,1)) = 1;
                                        Trace{1,g}(k,1) = Trace{1,g}(k,1)*110/1000;
                                        Trace{1,g}(k,2) = Trace{1,g}(k,2)*110/1000;
                                        Trace{1,g}(k,3) = Trace{1,g}(k,3)*0.2;
                                    end
                                    Chr(n).RoG = rog(Trace{1,g});
                                    Chr(n).Age = q;
                                    Chr(n).TraceNo = g;
                                    Chr(n).TracesInTerritory = 4;
                                    if SegChannel == 647
                                        Chr(n).L2 = L2; %561 segmentation = L1, 647 segmentation = L2
                                    else
                                        Chr(n).L1 = L1;
                                    end
                                    Chr(n).nucL = nucL;
                                    Chr(n).TraceName = TraceName;
                                end
                            end
                        end
                        if length(Trace) ==3
                            for g = 1:3
                                if length(TAD_id{g})>1 && length(TAD_id{1,g}) == length(unique(TAD_id{1,g}))
                                    n = n+1;
                                    Chr(n).x = zeros(TotalTADNum,1);
                                    Chr(n).y = zeros(TotalTADNum,1);
                                    Chr(n).z = zeros(TotalTADNum,1);
                                    Chr(n).r = zeros(TotalTADNum,1);
                                    
                                    for k = 1: length(TAD_id{1,g})
                                        Chr(n).x(TAD_id{1,g}(k,1)) = Trace{1,g}(k,1)*110/1000; % convert to um from pixel units
                                        Chr(n).y(TAD_id{1,g}(k,1)) = sz*110/1000-Trace{1,g}(k,2)*110/1000;
                                        Chr(n).z(TAD_id{1,g}(k,1)) = Trace{1,g}(k,3)*0.2;
                                        Chr(n).r(TAD_id{1,g}(k,1)) = 1;
                                        Trace{1,g}(k,1) = Trace{1,g}(k,1)*110/1000;
                                        Trace{1,g}(k,2) = Trace{1,g}(k,2)*110/1000;
                                        Trace{1,g}(k,3) = Trace{1,g}(k,3)*0.2;
                                    end
                                    Chr(n).RoG = rog(Trace{1,g});
                                    Chr(n).Age = q;
                                    Chr(n).TraceNo = g;
                                    Chr(n).TracesInTerritory = 3;
                                    if SegChannel == 647
                                        Chr(n).L2 = L2; %561 segmentation = L1, 647 segmentation = L2
                                    else
                                        Chr(n).L1 = L1;
                                    end
                                    Chr(n).nucL = nucL;
                                    Chr(n).TraceName = TraceName;
                                end
                            end
                        end
                        if length(Trace) ==2
                            for g = 1:2
                                if length(TAD_id{g})>1 && length(TAD_id{1,g}) == length(unique(TAD_id{1,g}))
                                    n = n+1;
                                    Chr(n).x = zeros(TotalTADNum,1);
                                    Chr(n).y = zeros(TotalTADNum,1);
                                    Chr(n).z = zeros(TotalTADNum,1);
                                    Chr(n).r = zeros(TotalTADNum,1);
                                    
                                    for k = 1: length(TAD_id{1,g})
                                        Chr(n).x(TAD_id{1,g}(k,1)) = Trace{1,g}(k,1)*110/1000; % convert to um from pixel units
                                        Chr(n).y(TAD_id{1,g}(k,1)) = sz*110/1000-Trace{1,g}(k,2)*110/1000;
                                        Chr(n).z(TAD_id{1,g}(k,1)) = Trace{1,g}(k,3)*0.2;
                                        Chr(n).r(TAD_id{1,g}(k,1)) = 1;
                                        Trace{1,g}(k,1) = Trace{1,g}(k,1)*110/1000;
                                        Trace{1,g}(k,2) = Trace{1,g}(k,2)*110/1000;
                                        Trace{1,g}(k,3) = Trace{1,g}(k,3)*0.2;
                                    end
                                    Chr(n).RoG = rog(Trace{1,g});
                                    Chr(n).Age = q;
                                    Chr(n).TraceNo = g;
                                    Chr(n).TracesInTerritory = 2;
                                    if SegChannel == 647
                                        Chr(n).L2 = L2; %561 segmentation = L1, 647 segmentation = L2
                                    else
                                        Chr(n).L1 = L1;
                                    end
                                    Chr(n).nucL = nucL;
                                    Chr(n).TraceName = TraceName;
                                end
                            end
                        end
                        if length(Trace) ==1 && iscell(Trace) && length(TAD_id{1})>1 && length(TAD_id{1,1}) == length(unique(TAD_id{1,1}))
                            n = n+1;
                            Chr(n).x = zeros(TotalTADNum,1);
                            Chr(n).y = zeros(TotalTADNum,1);
                            Chr(n).z = zeros(TotalTADNum,1);
                            Chr(n).r = zeros(TotalTADNum,1);
                            
                            for k = 1: length(TAD_id{1,1})
                                Chr(n).x(TAD_id{1,1}(k,1)) = Trace{1,1}(k,1)*110/1000; % convert to um from pixel units
                                Chr(n).y(TAD_id{1,1}(k,1)) = sz*110/1000-Trace{1,1}(k,2)*110/1000;
                                Chr(n).z(TAD_id{1,1}(k,1)) = Trace{1,1}(k,3)*0.2;
                                Chr(n).r(TAD_id{1,1}(k,1)) = 1;
                                Trace{1,1}(k,1) = Trace{1,1}(k,1)*110/1000;
                                Trace{1,1}(k,2) = Trace{1,1}(k,2)*110/1000;
                                Trace{1,1}(k,3) = Trace{1,1}(k,3)*0.2;
                            end
                            Chr(n).RoG = rog(Trace{1,1});
                            Chr(n).Age = q;
                            Chr(n).TraceNo = 1;
                            Chr(n).TracesInTerritory = 1;
                            if SegChannel == 647
                                Chr(n).L2 = L2; %561 segmentation = L1, 647 segmentation = L2
                            else
                                Chr(n).L1 = L1;
                            end                      
                            Chr(n).nucL = nucL;
                            Chr(n).TraceName = TraceName;
                        end
                        flag = 0;
                    end
                end
            end
        end
    end
end


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
plot(c, r, 'yo')
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







%%
save([path '\MeanPairDistanceMatrix' num2str(SegChannel) '.mat'], 'Mean');
save([path '\MedianPairDistanceMatrix' num2str(SegChannel) '.mat'], 'Median');
save([path '\StdPairDistanceMatrix' num2str(SegChannel) '.mat'], 'Std');
save([path '\SEMPairDistanceMatrix' num2str(SegChannel) '.mat'], 'SEM');
save([path '\NofDataPairDistanceMatrix' num2str(SegChannel) '.mat'], 'NofData');
save([path '\AllChromosomes' num2str(SegChannel) '.mat'], 'Chr', '-v7.3');
save([path '\DisListAll' num2str(SegChannel) '.mat'], 'DisListAll');


toc


%%

