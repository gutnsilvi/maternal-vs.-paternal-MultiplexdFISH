%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script compiles all traces in a mat-file, filtered by age and
% extracts mean, median, SEM, Std 
% the script only uses traces that meet following conditions: 
%   - 2 traces per nucleus (1 from HI, 1 from N2) 
%   Input: 
%   chromosome traces from crosses of N2 and HI
%
%   Output: (separated for both strains)
%   AllChromosomes file
%   mean pairwise distances
%   median pairwise distances
%   SEM pairwise distances
%   Std pairwise distances
%   Distance List of all pairwise distances of all chromosomes
%   Number of datapoints for all pairwise distances of all chromosomes
% 
% create folders before running: traces

%%%%%%%%%%%%%%%%%%%%%%%%%%%


tic;
clear all
close all

Folder{1} = 'traces1';
Folder{2} = 'traces2';
Folder{3} = 'traces3';
Folder{4} = 'traces4';
Folder{5} = 'traces5';



NumSamples = 38;
AgeStart = 2; %change these to specify embryo age
AgeStop = 40;
MaxROI = 50; %change to specify largest roi# 
MaxFOV = 7; %change to specify largest FOV# 

AgeRange = ([ num2str(AgeStart) 'to' num2str(AgeStop) 'cell']); 
SegChannel = 561;
sz = 1608;

%mkdir(AgeRange)

TotalTADNum = 22;
n = 0;      %counter chromosomes
flag = 0;

excluded_nuc = 0; %counts how many are excluded based on to many traces per nucleus

excluded_nuc_HI = 0; %counts how many are excluded based on to many HI traces
excluded_nuc_N2 = 0; %counts how many excluded based on to many N2 traces 


for ii = 1:length(Folder)
    for i = 0:NumSamples
        for q = AgeStart:AgeStop
            for b = 1:MaxFOV
                for j = 1:MaxROI
                    if exist([Folder{ii} '/' num2str(q) '_cell_Traces_Sample' num2str(i) '_FOV' num2str(b) '_roi' num2str(SegChannel) '_' num2str(j) '.mat'])
                        load([Folder{ii} '/' num2str(q) '_cell_Traces_Sample' num2str(i) '_FOV' num2str(b) '_roi' num2str(SegChannel) '_' num2str(j) '.mat']);
                        TraceName = num2str([Folder{ii} '/' num2str(q) '_cell_Traces_Sample' num2str(i) '_FOV' num2str(b) '_roi' num2str(SegChannel) '_' num2str(j) '.mat']);
                        flag = 1;                
                    end
                    if flag == 1
                        
                        numHI_traces = length(find(~cellfun(@isempty,Trace_HI)));
                        numN2_traces = length(find(~cellfun(@isempty,Trace_N2)));
                        numboth_traces = length(find(~cellfun(@isempty,Trace_both)));   %%to exclude nuclei where the segmentation was not correctly assigning expected NofTraces
                        
                        
                        if numHI_traces + numN2_traces + numboth_traces <= 4 && numHI_traces == 1 & numN2_traces == 1 & numboth_traces == 0 %sort only nuclei with exactly 2 traces - 1 per strain
                                                                                               
                         n = n +1;   
                                                                                                          
                         idx_terr = size(Trace_HI,1);        %number of territories detected
                            
                         for iii = 1: idx_terr
                             
                            
                            [TraceN2_1, TraceN2_2, TraceN2_3] = Trace_N2{iii, 1:3};     % traces per territory
                            [TraceHI_1, TraceHI_2, TraceHI_3] = Trace_HI{iii, 1:3};     %traces per territory
                            
                            
                            HI_trace = int8(~isempty(TraceHI_3)) + int8(~isempty(TraceHI_2)) + int8(~isempty(TraceHI_1));
                            N2_trace = int8(~isempty(TraceN2_3)) + int8(~isempty(TraceN2_2)) + int8(~isempty(TraceN2_1));        %counts traces of each category in territory we are looking at now
                            
                            
                            if HI_trace == 1 & N2_trace == 0                        %territory contains 1 HI trace
                              if ~isempty(TraceHI_1)
                                    g = 1;
                                else g = 2;
                                end
                                
                                if length(TAD_id_HI{iii,g})>1 && length(TAD_id_HI{iii,g}) == length(unique(TAD_id_HI{iii,g})) 
                                    
                                   
                                    Chr_HI(n).x = zeros(TotalTADNum,1);
                                    Chr_HI(n).y = zeros(TotalTADNum,1);
                                    Chr_HI(n).z = zeros(TotalTADNum,1);
                                    Chr_HI(n).r = zeros(TotalTADNum,1);
%                                     
                                    

                                    for k = 1: length(TAD_id_HI{iii,g})
                                        Chr_HI(n).x(TAD_id_HI{iii,g}(k,1)) = Trace_HI{iii,g}(k,1)*110/1000; % convert to um from pixel units
                                        Chr_HI(n).y(TAD_id_HI{iii,g}(k,1)) = sz*110/1000-Trace_HI{iii,g}(k,2)*110/1000;
                                        Chr_HI(n).z(TAD_id_HI{iii,g}(k,1)) = Trace_HI{iii,g}(k,3)*0.2;
                                        Chr_HI(n).r(TAD_id_HI{iii,g}(k,1)) = 1;
%                                         
                                    end
                                    
                                    Chr_HI(n).RoG = rog(Trace_HI{iii,g});
                                    Chr_HI(n).Age = q;
                                    Chr_HI(n).TraceNo = g;
                                    Chr_HI(n).TracesInTerritory = 1;
%                                    
                                    Chr_HI(n).nuc = nuc;
                                    Chr_HI(n).TraceName = [TraceName '_nuc' num2str(nuc) '_Terr' num2str(iii)] ;
                                
                                    
                                
                               
                                end
                            end
                            
                            
                            
                            
                            if HI_trace == 0 & N2_trace == 1        %territory contains 1 N2 trace
                                
                               if ~isempty(TraceN2_1)
                                    g = 1;
                                else g = 2;
                                end
                                
                                
                                if length(TAD_id_N2{iii,g})>1 && length(TAD_id_N2{iii,g}) == length(unique(TAD_id_N2{iii,g})) 
                                    
                                   
                                    Chr_N2(n).x = zeros(TotalTADNum,1);
                                    Chr_N2(n).y = zeros(TotalTADNum,1);
                                    Chr_N2(n).z = zeros(TotalTADNum,1);
                                    Chr_N2(n).r = zeros(TotalTADNum,1);
%                                    

                                    for k = 1: length(TAD_id_N2{iii,g})
                                        Chr_N2(n).x(TAD_id_N2{iii,g}(k,1)) = Trace_N2{iii,g}(k,1)*110/1000; % convert to um from pixel units
                                        Chr_N2(n).y(TAD_id_N2{iii,g}(k,1)) = sz*110/1000-Trace_N2{iii,g}(k,2)*110/1000;
                                        Chr_N2(n).z(TAD_id_N2{iii,g}(k,1)) = Trace_N2{iii,g}(k,3)*0.2;
                                        Chr_N2(n).r(TAD_id_N2{iii,g}(k,1)) = 1;
%                                         
                                    end
                                    Chr_N2(n).RoG = rog(Trace_N2{iii,g});
                                    Chr_N2(n).Age = q;
                                    Chr_N2(n).TraceNo = g;
                                    Chr_N2(n).TracesInTerritory = 1;
%                                     
                                    Chr_N2(n).nuc = nuc;
                                    Chr_N2(n).TraceName = [TraceName '_nuc' num2str(nuc) '_Terr' num2str(iii)];
                                
                                
                                end
                            end
                            
                            
                           
                            if HI_trace ==  N2_trace & N2_trace == 1        %%HI and N2 trace in the same territory     
                                
                                if ~isempty(TraceHI_1)
                                    g = 1;
                                else g = 2;
                                end
                                    
                                if length(TAD_id_HI{iii,g})>1 && length(TAD_id_HI{iii,g}) == length(unique(TAD_id_HI{iii,g})) 
                                    
                                    
                                    Chr_HI(n).x = zeros(TotalTADNum,1);
                                    Chr_HI(n).y = zeros(TotalTADNum,1);
                                    Chr_HI(n).z = zeros(TotalTADNum,1);
                                    Chr_HI(n).r = zeros(TotalTADNum,1);
%                                    

                                    for k = 1: length(TAD_id_HI{iii,g})
                                        Chr_HI(n).x(TAD_id_HI{iii,g}(k,1)) = Trace_HI{iii,g}(k,1)*110/1000; % convert to um from pixel units
                                        Chr_HI(n).y(TAD_id_HI{iii,g}(k,1)) = sz*110/1000-Trace_HI{iii,g}(k,2)*110/1000;
                                        Chr_HI(n).z(TAD_id_HI{iii,g}(k,1)) = Trace_HI{iii,g}(k,3)*0.2;
                                        Chr_HI(n).r(TAD_id_HI{iii,g}(k,1)) = 1;
%                                         
                                    end
                                    
                                    Chr_HI(n).RoG = rog(Trace_HI{iii,g});
                                    Chr_HI(n).Age = q;
                                    Chr_HI(n).TraceNo = g;
                                    Chr_HI(n).TracesInTerritory = 2;
%                                    
                                    Chr_HI(n).nuc = nuc;
                                    Chr_HI(n).TraceName = [TraceName '_nuc' num2str(nuc) '_Terr' num2str(iii)];
                                    
                                    
                                 
                                    
                                end
                                    
                                
                                if ~isempty(TraceN2_1)
                                    g = 1;
                                else g = 2;
                                end
                                
                                if length(TAD_id_N2{iii,g})>1 && length(TAD_id_N2{iii,g}) == length(unique(TAD_id_N2{iii,g})) 
                                    
                                    
                                    Chr_N2(n).x = zeros(TotalTADNum,1);
                                    Chr_N2(n).y = zeros(TotalTADNum,1);
                                    Chr_N2(n).z = zeros(TotalTADNum,1);
                                    Chr_N2(n).r = zeros(TotalTADNum,1);
%                                     

                                    for k = 1: length(TAD_id_N2{iii,g})
                                        Chr_N2(n).x(TAD_id_N2{iii,g}(k,1)) = Trace_N2{iii,g}(k,1)*110/1000; % convert to um from pixel units
                                        Chr_N2(n).y(TAD_id_N2{iii,g}(k,1)) = sz*110/1000-Trace_N2{iii,g}(k,2)*110/1000;
                                        Chr_N2(n).z(TAD_id_N2{iii,g}(k,1)) = Trace_N2{iii,g}(k,3)*0.2;
                                        Chr_N2(n).r(TAD_id_N2{iii,g}(k,1)) = 1;
%                                         
                                    end
                                    Chr_N2(n).RoG = rog(Trace_N2{iii,g});
                                    Chr_N2(n).Age = q;
                                    Chr_N2(n).TraceNo = g;
                                    Chr_N2(n).TracesInTerritory = 2;
%                                     
                                    Chr_N2(n).nuc = nuc;
                                    Chr_N2(n).TraceName = [TraceName '_nuc' num2str(nuc) '_Terr' num2str(iii)];
                                    
                                    
                                     
                                end
                            end
                         end
                        
                            
                            
                             
                            
                   
                            
                        elseif numHI_traces + numN2_traces > 4 
                            excluded_nuc = excluded_nuc + numHI_traces + numN2_traces;
                        elseif numHI_traces > 1 | numN2_traces > 1 
                            excluded_nuc_HI = excluded_nuc_HI + numHI_traces;
                        
                            excluded_nuc_N2 = excluded_nuc_N2 + numN2_traces;
                         
                        end
                        
                   
                        
                    flag = 0;
                    end
                end
            end
        end
    end
end 

%% remove empty rows and corresponding row in reciprocal dataset

toremove_N2 = find(cellfun(@isempty,{Chr_HI.Age}));
toremove_HI = find(cellfun(@isempty,{Chr_N2.Age}));

Chr_HI(toremove_HI) = [];

Chr_N2(toremove_N2) = [];

toremove_N2 = find(cellfun(@isempty,{Chr_N2.Age}));
Chr_N2(toremove_N2) = [];

toremove_HI = find(cellfun(@isempty,{Chr_HI.Age}));
Chr_HI(toremove_HI) = [];

if length(Chr_HI) ~= length(Chr_N2)
    if length(Chr_HI) < length(Chr_N2)
        Chr_N2(end) = [];
    else Chr_HI(end) = [];
    end
end

    
%%
Mean_HI = zeros(TotalTADNum,TotalTADNum);
Std_HI = zeros(TotalTADNum,TotalTADNum);
Median_HI = zeros(TotalTADNum,TotalTADNum);

Mean_N2 = zeros(TotalTADNum,TotalTADNum);
Std_N2 = zeros(TotalTADNum,TotalTADNum);
Median_N2 = zeros(TotalTADNum,TotalTADNum);

Mean_hom = zeros(TotalTADNum,TotalTADNum);
Std_hom = zeros(TotalTADNum,TotalTADNum);
Median_hom = zeros(TotalTADNum,TotalTADNum);
SEM_hom = zeros(TotalTADNum,TotalTADNum);
NofData_hom = zeros(TotalTADNum,TotalTADNum);



for i = 1:TotalTADNum
    for j = 1:TotalTADNum
        DisList_HI = [];
        for k = 1:length(Chr_HI) 
            if Chr_HI(k).r(i) == 1 && Chr_HI(k).r(j) == 1
                
                    DisList_HI = [DisList_HI ((Chr_HI(k).x(i)-Chr_HI(k).x(j))^2+(Chr_HI(k).y(i)-Chr_HI(k).y(j))^2+(Chr_HI(k).z(i)-Chr_HI(k).z(j))^2)^0.5];
                
            end
        end
        Mean_HI(i,j) = mean(DisList_HI);
        Std_HI(i,j) = std(DisList_HI);
        SEM_HI(i,j) = std(DisList_HI)/(length(DisList_HI))^0.5;
        NofData_HI(i,j) = length(DisList_HI);
        DisListAll_HI{i,j} = DisList_HI;
        Median_HI(i,j) = median(DisList_HI);
    end
end


for i = 1:TotalTADNum
    for j = 1:TotalTADNum
        DisList_N2 = [];
        for k = 1:length(Chr_N2) 
            if Chr_N2(k).r(i) == 1 && Chr_N2(k).r(j) == 1
                
                    DisList_N2 = [DisList_N2 ((Chr_N2(k).x(i)-Chr_N2(k).x(j))^2+(Chr_N2(k).y(i)-Chr_N2(k).y(j))^2+(Chr_N2(k).z(i)-Chr_N2(k).z(j))^2)^0.5];
                
            end
        end
        Mean_N2(i,j) = mean(DisList_N2);
        Std_N2(i,j) = std(DisList_N2);
        SEM_N2(i,j) = std(DisList_N2)/(length(DisList_N2))^0.5;
        NofData_N2(i,j) = length(DisList_N2);
        DisListAll_N2{i,j} = DisList_N2;
        Median_N2(i,j) = median(DisList_N2);
    end
end


for i = 1:TotalTADNum
    for j = 1:TotalTADNum
        DisList_hom = [];
        for k = 1:length(Chr_HI) 
            if Chr_HI(k).r(i) == 1 && Chr_N2(k).r(j) == 1  
                
                    DisList_hom = [DisList_hom ((Chr_HI(k).x(i)-Chr_N2(k).x(j))^2+(Chr_HI(k).y(i)-Chr_N2(k).y(j))^2+(Chr_HI(k).z(i)-Chr_N2(k).z(j))^2)^0.5];
                
            end
        end
        Mean_hom(i,j) = mean(DisList_hom);
        Std_hom(i,j) = std(DisList_hom);
        SEM_hom(i,j) = std(DisList_hom)/(length(DisList_hom))^0.5;
        NofData_hom(i,j) = length(DisList_hom);
        DisListAll_hom{i,j} = DisList_hom;
        Median_hom(i,j) = median(DisList_hom);
    end
end

Mean_hom(11, :) = [];
Mean_hom(:,11) = [];

%% make figures

figure2 = figure;
colormap(hot);
% Create axes
axes1 = axes('Parent',figure2);
hold(axes1,'on');
% Create image
image(Mean_HI,'Parent',axes1,'CDataMapping','scaled');
title(['Mean HI pairwise distance (um) -'  AgeRange]);
% Uncomment the following line to preserve the X-limits of the axes
xlim(axes1,[0.5 22.5]);
% Uncomment the following line to preserve the Y-limits of the axes
ylim(axes1,[0.5 22.5]);
box(axes1,'on');
axis(axes1,'ij');
colorbar('peer',axes1);
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
image(Mean_N2,'Parent',axes1,'CDataMapping','scaled');
title(['Mean N2 pairwise distance (um) -'  AgeRange]);
% Uncomment the following line to preserve the X-limits of the axes
xlim(axes1,[0.5 22.5]);
% Uncomment the following line to preserve the Y-limits of the axes
ylim(axes1,[0.5 22.5]);
box(axes1,'on');
axis(axes1,'ij');
colorbar('peer',axes1);
caxis([0 2])
PlotProp
axis square
%print('Mean_early','-dsvg')

%

figure4 = figure;
colormap(hot);
% Create axes
axes1 = axes('Parent',figure4);
hold(axes1,'on');
% Create image
image(Median_HI,'Parent',axes1,'CDataMapping','scaled');
title(['Median pairwise distance HI (um) -'  AgeRange]);
% Uncomment the following line to preserve the X-limits of the axes
xlim(axes1,[0.5 22.5]);
% Uncomment the following line to preserve the Y-limits of the axes
ylim(axes1,[0.5 22.5]);
box(axes1,'on');
axis(axes1,'ij');
colorbar('peer',axes1);
caxis([0 2])
PlotProp
axis square
%print('Median_early','-dsvg')

figure5 = figure;
colormap(hot);
% Create axes
axes1 = axes('Parent',figure5);
hold(axes1,'on');
% Create image
image(Median_N2,'Parent',axes1,'CDataMapping','scaled');
title(['Median pairwise distance N2 (um) -'  AgeRange]);
% Uncomment the following line to preserve the X-limits of the axes
xlim(axes1,[0.5 22.5]);
% Uncomment the following line to preserve the Y-limits of the axes
ylim(axes1,[0.5 22.5]);
box(axes1,'on');
axis(axes1,'ij');
colorbar('peer',axes1);
caxis([0 2])
PlotProp
axis square
%print('Median_early','-dsvg')


figure(6)
imagesc(SEM_HI)
colorbar
title(['SEM of pairwise distance HI - ' AgeRange]);
PlotProp
axis square
%print('SEM_early','-dsvg')

figure(7)
imagesc(SEM_N2)
colorbar
title(['SEM of pairwise distance N2 - ' AgeRange]);
PlotProp
axis square
%print('SEM_early','-dsvg')

figure(8)
imagesc(Std_HI)
colorbar
title(['Std of pairwise distance HI - ' AgeRange]);
PlotProp
axis square
%print('STD_early','-dsvg')


figure(9)
imagesc(Std_N2)
colorbar
title(['Std of pairwise distance N2 - ' AgeRange]);
PlotProp
axis square
%print('STD_early','-dsvg')

figure(10)
imagesc(NofData_HI)
colorbar
title(['Number of measurements HI - ' AgeRange]);
PlotProp
axis square
%print('NofDATA_early','-dsvg')

figure(11)
imagesc(NofData_N2)
colorbar
title(['Number of measurements N2 - ' AgeRange]);
PlotProp
axis square
%print('NofDATA_early','-dsvg')

figure(12)
histogram([Chr_HI.TracesInTerritory])
title('HI segmentation results')
ylabel('number of traces')
xlabel('total traces in territory volume')
PlotProp

TracesInTerritory_HI = [Chr_HI.TracesInTerritory].';

TracesPerSegment = [];
TracesPerSegment(1,1) = sum(TracesInTerritory_HI(:) == 1)/length(Chr_HI)*100;
TracesPerSegment(1,2) = sum(TracesInTerritory_HI(:) == 2)/length(Chr_HI)*100;
TracesPerSegment(1,3) = sum(TracesInTerritory_HI(:) == 3)/length(Chr_HI)*100;
TracesPerSegment(1,4) = sum(TracesInTerritory_HI(:) == 4)/length(Chr_HI)*100;
figure(13)
bar(TracesPerSegment)
title('HI segmentation results')
ylabel('% of total traces')
xlabel('traces within territory volume')
ylim([0 100])
PlotProp
%print('wt_segmentation_early','-dsvg')

figure(14)
histogram([Chr_N2.TracesInTerritory])
title('N2 segmentation results')
ylabel('number of traces')
xlabel('total traces in territory volume')
PlotProp

TracesInTerritory_N2 = [Chr_N2.TracesInTerritory].';

TracesPerSegment = [];
TracesPerSegment(1,1) = sum(TracesInTerritory_N2(:) == 1)/length(Chr_N2)*100;
TracesPerSegment(1,2) = sum(TracesInTerritory_N2(:) == 2)/length(Chr_N2)*100;
TracesPerSegment(1,3) = sum(TracesInTerritory_N2(:) == 3)/length(Chr_N2)*100;
TracesPerSegment(1,4) = sum(TracesInTerritory_N2(:) == 4)/length(Chr_N2)*100;
figure(15)
bar(TracesPerSegment)
title('N2 segmentation results')
ylabel('% of total traces')
xlabel('traces within territory volume')
ylim([0 100])
PlotProp
%print('wt_segmentation_early','-dsvg')




%% figure for homolog comparison


figure16 = figure;
colormap(hot);
% Create axes
axes1 = axes('Parent',figure16);
hold(axes1,'on');
% Create image
image(Mean_hom,'Parent',axes1,'CDataMapping','scaled');
title(['Mean HI and N2 pairwise distance (um) -'  AgeRange]);
% Uncomment the following line to preserve the X-limits of the axes
xlim(axes1,[0.5 21.5]);
% Uncomment the following line to preserve the Y-limits of the axes
ylim(axes1,[0.5 21.5]);
box(axes1,'on');
axis(axes1,'ij');
colorbar('peer',axes1);
caxis([0 4])
PlotProp
axis square
%print('Mean_early','-dsvg')

figure17 = figure;
colormap(hot);
% Create axes
axes1 = axes('Parent',figure17);
hold(axes1,'on');
% Create image
image(Median_hom,'Parent',axes1,'CDataMapping','scaled');
title(['Median pairwise distance HI and N2 (um) -'  AgeRange]);
% Uncomment the following line to preserve the X-limits of the axes
xlim(axes1,[0.5 22.5]);
% Uncomment the following line to preserve the Y-limits of the axes
ylim(axes1,[0.5 22.5]);
box(axes1,'on');
axis(axes1,'ij');
colorbar('peer',axes1);
caxis([0 4])
PlotProp
axis square
%print('Median_early','-dsvg')


figure(18)
imagesc(SEM_hom)
colorbar
title(['SEM of pairwise distance HI and N2 - ' AgeRange]);
PlotProp
axis square
%print('SEM_early','-dsvg')

figure(9)
imagesc(Std_hom)
colorbar
title(['Std of pairwise distance HI and N2 - ' AgeRange]);
PlotProp
axis square
%print('STD_early','-dsvg')

figure(10)
imagesc(NofData_hom)
colorbar
title(['Number of measurements HI and N2 - ' AgeRange]);
PlotProp
axis square
%print('NofDATA_early','-dsvg')

%%

edges = linspace(0, 6, 23);
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





%%
save([AgeRange '\homologs\MeanPairDistanceMatrix_HI' num2str(SegChannel) '.mat'], 'Mean_HI');
save([AgeRange '\homologs\MedianPairDistanceMatrix_HI' num2str(SegChannel) '.mat'], 'Median_HI');
save([AgeRange '\homologs\StdPairDistanceMatrix_HI' num2str(SegChannel) '.mat'], 'Std_HI');
save([AgeRange '\homologs\SEMPairDistanceMatrix_HI' num2str(SegChannel) '.mat'], 'SEM_HI');
save([AgeRange '\homologs\NofDataPairDistanceMatrix_HI' num2str(SegChannel) '.mat'], 'NofData_HI');
save([AgeRange '\homologs\AllChromosomes_HI' num2str(SegChannel) '.mat'], 'Chr_HI', '-v7.3');
save([AgeRange '\homologs\DisListAll_HI' num2str(SegChannel) '.mat'], 'DisListAll_HI');


save([AgeRange '\homologs\MeanPairDistanceMatrix_N2' num2str(SegChannel) '.mat'], 'Mean_N2');
save([AgeRange '\homologs\MedianPairDistanceMatrix_N2' num2str(SegChannel) '.mat'], 'Median_N2');
save([AgeRange '\homologs\StdPairDistanceMatrix_N2' num2str(SegChannel) '.mat'], 'Std_N2');
save([AgeRange '\homologs\SEMPairDistanceMatrix_N2' num2str(SegChannel) '.mat'], 'SEM_N2');
save([AgeRange '\homologs\NofDataPairDistanceMatrix_N2' num2str(SegChannel) '.mat'], 'NofData_N2');
save([AgeRange '\homologs\AllChromosomes_N2' num2str(SegChannel) '.mat'], 'Chr_N2', '-v7.3');
save([AgeRange '\homologs\DisListAll_N2' num2str(SegChannel) '.mat'], 'DisListAll_N2');



toc