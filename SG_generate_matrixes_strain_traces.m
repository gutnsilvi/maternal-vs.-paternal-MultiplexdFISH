%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script compiles all traces in a mat-file, filtered by age and
%%extracts mean, median, SEM, Std
%
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


NumSamples = 38;
AgeStart = 17; %change these to specify embryo age
AgeStop = 24;
MaxROI = 40; %change to specify largest roi# 
MaxFOV = 7; %change to specify largest FOV# 

AgeRange = ([ num2str(AgeStart) 'to' num2str(AgeStop) 'cell']); 
SegChannel = 561;
sz = 1608;

mkdir(AgeRange)

TotalTADNum = 22;
n = 0;      %counter for n2 chromosomes
h = 0;      %counter for Hi chromosomes 
p = 0;      %count homolog comparisons
flag = 0;

for ii = 1:length(Folder)
    for s = 0:NumSamples
        for q = AgeStart:AgeStop
            for b = 1:MaxFOV
                for j = 1:MaxROI
                    if exist([Folder{ii} '/' num2str(q) '_cell_Traces_Sample' num2str(s) '_FOV' num2str(b) '_roi' num2str(SegChannel) '_' num2str(j) '.mat'])
                        load([Folder{ii} '/' num2str(q) '_cell_Traces_Sample' num2str(s) '_FOV' num2str(b) '_roi' num2str(SegChannel) '_' num2str(j) '.mat']);
                        TraceName = num2str([Folder{ii} '/' num2str(q) '_cell_Traces_Sample' num2str(i) '_FOV' num2str(b) '_roi' num2str(SegChannel) '_' num2str(j) '.mat']);
                        flag = 1;                
                    end
                    if flag == 1
                        
                        numHI_traces = length(find(~cellfun(@isempty,Trace_HI)));
                        numN2_traces = length(find(~cellfun(@isempty,Trace_N2)));   %%to exclude nuclei where the segmentation was not correctly assigning expected NofTraces
                        numboth_traces = length(find(~cellfun(@isempty,Trace_both)));
                        
                        if numHI_traces + numN2_traces + numboth_traces <= 4 & numHI_traces <= 2 & numN2_traces <= 2
                     
                            %exludes all results where to many traces where found per nucleus, can only be maximum 2 for each strain
                                                                                                          
                         idx_terr = size(Trace_HI,1);        %number of territories detected
                            
                         for iii = 1: idx_terr
                            
                             if length(Trace_N2(iii,:)) == 4
                            
                                [TraceN2_1, TraceN2_2, TraceN2_3, TraceN2_4] = Trace_N2{iii, 1:4};  % traces per territory
                                 N2_trace = int8(~isempty(TraceN2_2)) + int8(~isempty(TraceN2_1))+ int8(~isempty(TraceN2_3)) + int8(~isempty(TraceN2_4));
                             else  [TraceN2_1, TraceN2_2, TraceN2_3] = Trace_N2{iii, 1:3};
                                   N2_trace = int8(~isempty(TraceN2_2)) + int8(~isempty(TraceN2_1))+ int8(~isempty(TraceN2_3));
                             end
                           
                             if length(Trace_HI(iii,:)) == 4
                             
                            [TraceHI_1, TraceHI_2, TraceHI_3, TraceHI_4] = Trace_HI{iii, 1:4}; 
                            HI_trace = int8(~isempty(TraceHI_2)) + int8(~isempty(TraceHI_1)) + int8(~isempty(TraceHI_3)) + int8(~isempty(TraceHI_4));
                            
                             else [TraceHI_1, TraceHI_2, TraceHI_3] = Trace_HI{iii, 1:3};   %traces per territory
                                    HI_trace = int8(~isempty(TraceHI_2)) + int8(~isempty(TraceHI_1)) + int8(~isempty(TraceHI_3));
                                    %counts traces of each category in territory we are looking at now
                             end
                            if HI_trace == 1 & N2_trace == 0                        %territory contains 1 HI trace
                               
                               t = [];              

                                for i = 1:size(Trace_HI,2)

                                    if ~isempty(Trace_HI{iii,i})
                                        t = [t, i];
                                    end
                                end
                                
                                g = t; 

                               
                                if length(TAD_id_HI{iii,g})>1 && length(TAD_id_HI{iii,g}) == length(unique(TAD_id_HI{iii,g}))
                                    h = h+1;
                                    Chr_HI(h).x = zeros(TotalTADNum,1);
                                    Chr_HI(h).y = zeros(TotalTADNum,1);
                                    Chr_HI(h).z = zeros(TotalTADNum,1);
                                    Chr_HI(h).r = zeros(TotalTADNum,1);

                                    for k = 1: length(TAD_id_HI{iii,g})
                                        Chr_HI(h).x(TAD_id_HI{iii,g}(k,1)) = Trace_HI{iii,g}(k,1)*110/1000; % convert to um from pixel units
                                        Chr_HI(h).y(TAD_id_HI{iii,g}(k,1)) = sz*110/1000-Trace_HI{iii,g}(k,2)*110/1000;
                                        Chr_HI(h).z(TAD_id_HI{iii,g}(k,1)) = Trace_HI{iii,g}(k,3)*0.2;
                                        Chr_HI(h).r(TAD_id_HI{iii,g}(k,1)) = 1;
%                                         
                                    end
                                    
                                    Chr_HI(h).RoG = rog(Trace_HI{iii,g});
                                    Chr_HI(h).Age = q;
                                    Chr_HI(h).TraceNo = g;
                                    Chr_HI(h).TracesInTerritory = 1;
%                                    
                                    Chr_HI(h).nuc = nuc;
                                    Chr_HI(h).TraceName = [TraceName '_nuc' num2str(nuc) '_Terr' num2str(iii)] ;
                                
                               
                                end
                            end
                            
                            if HI_trace == 0 & N2_trace == 1        %territory contains 1 N2 trace
                               
                                t = [];              

                                for i = 1:size(Trace_N2,2)

                                    if ~isempty(Trace_N2{iii,i})
                                        t = [t, i];
                                    end
                                end
                                
                                g = t;
                                
                                if length(TAD_id_N2{iii,g})>1 && length(TAD_id_N2{iii,g}) == length(unique(TAD_id_N2{iii,g}))
                                    n = n+1;
                                    Chr_N2(n).x = zeros(TotalTADNum,1);
                                    Chr_N2(n).y = zeros(TotalTADNum,1);
                                    Chr_N2(n).z = zeros(TotalTADNum,1);
                                    Chr_N2(n).r = zeros(TotalTADNum,1);

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
                              
                               t = [];              

                                for i = 1:size(Trace_HI,2)

                                    if ~isempty(Trace_HI{iii,i})
                                        t = [t, i];
                                    end
                                end
                                
                                g = t;
                               
                                
                                if length(TAD_id_HI{iii,g})>1 && length(TAD_id_HI{iii,g}) == length(unique(TAD_id_HI{iii,g}))
                                    h = h+1;
                                    Chr_HI(h).x = zeros(TotalTADNum,1);
                                    Chr_HI(h).y = zeros(TotalTADNum,1);
                                    Chr_HI(h).z = zeros(TotalTADNum,1);
                                    Chr_HI(h).r = zeros(TotalTADNum,1);

                                    for k = 1: length(TAD_id_HI{iii,g})
                                        Chr_HI(h).x(TAD_id_HI{iii,g}(k,1)) = Trace_HI{iii,g}(k,1)*110/1000; % convert to um from pixel units
                                        Chr_HI(h).y(TAD_id_HI{iii,g}(k,1)) = sz*110/1000-Trace_HI{iii,g}(k,2)*110/1000;
                                        Chr_HI(h).z(TAD_id_HI{iii,g}(k,1)) = Trace_HI{iii,g}(k,3)*0.2;
                                        Chr_HI(h).r(TAD_id_HI{iii,g}(k,1)) = 1;
%                                         
                                    end
                                    
                                    Chr_HI(h).RoG = rog(Trace_HI{iii,g});
                                    Chr_HI(h).Age = q;
                                    Chr_HI(h).TraceNo = g;
                                    Chr_HI(h).TracesInTerritory = 2;
%                                    
                                    Chr_HI(h).nuc = nuc;
                                    Chr_HI(h).TraceName = [TraceName '_nuc' num2str(nuc) '_Terr' num2str(iii)];
                                end
                                    
                                t = [];              

                                for i = 1:size(Trace_N2,2)

                                    if ~isempty(Trace_N2{iii,i})
                                        t = [t, i];
                                    end
                                end
                                
                                g = t;
                                
                                if length(TAD_id_N2{iii,g})>1 && length(TAD_id_N2{iii,g}) == length(unique(TAD_id_N2{iii,g}))
                                    n = n+1;
                                    Chr_N2(n).x = zeros(TotalTADNum,1);
                                    Chr_N2(n).y = zeros(TotalTADNum,1);
                                    Chr_N2(n).z = zeros(TotalTADNum,1);
                                    Chr_N2(n).r = zeros(TotalTADNum,1);

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
                             
                            
                             if HI_trace ==  N2_trace & N2_trace == 2        %%HI and N2 trace in the same territory     
                               
                               t = [];              

                                for i = 1:size(Trace_HI,2)

                                    if ~isempty(Trace_HI{iii,i})
                                        t = [t, i];
                                    end
                                end

                               
                                for iv = 1:length(t)
                                    g = t(iv);
                                
                                
                                    if length(TAD_id_HI{iii,g})>1 && length(TAD_id_HI{iii,g}) == length(unique(TAD_id_HI{iii,g}))
                                        h = h+1;
                                        Chr_HI(h).x = zeros(TotalTADNum,1);
                                        Chr_HI(h).y = zeros(TotalTADNum,1);
                                        Chr_HI(h).z = zeros(TotalTADNum,1);
                                        Chr_HI(h).r = zeros(TotalTADNum,1);

                                        for k = 1: length(TAD_id_HI{iii,g})
                                            Chr_HI(h).x(TAD_id_HI{iii,g}(k,1)) = Trace_HI{iii,g}(k,1)*110/1000; % convert to um from pixel units
                                            Chr_HI(h).y(TAD_id_HI{iii,g}(k,1)) = sz*110/1000-Trace_HI{iii,g}(k,2)*110/1000;
                                            Chr_HI(h).z(TAD_id_HI{iii,g}(k,1)) = Trace_HI{iii,g}(k,3)*0.2;
                                            Chr_HI(h).r(TAD_id_HI{iii,g}(k,1)) = 1;
%                                             
                                        end
                                    
                                        Chr_HI(h).RoG = rog(Trace_HI{iii,g});
                                        Chr_HI(h).Age = q;
                                        Chr_HI(h).TraceNo = g;
                                        Chr_HI(h).TracesInTerritory = 4;
%                                        
                                        Chr_HI(h).nuc = nuc;
                                        Chr_HI(h).TraceName = [TraceName '_nuc' num2str(nuc) '_Terr' num2str(iii)];
                                    end
                                end
                                   
                                t = [];              

                                for i = 1:size(Trace_N2,2)

                                    if ~isempty(Trace_N2{iii,i})
                                        t = [t, i];
                                    end
                                end
                                
                                for iv = 1:length(t)
                                    g = t(iv);
                                
                                    
                                    if length(TAD_id_N2{iii,g})>1 && length(TAD_id_N2{iii,g}) == length(unique(TAD_id_N2{iii,g}))
                                        n = n+1;
                                        Chr_N2(n).x = zeros(TotalTADNum,1);
                                        Chr_N2(n).y = zeros(TotalTADNum,1);
                                        Chr_N2(n).z = zeros(TotalTADNum,1);
                                        Chr_N2(n).r = zeros(TotalTADNum,1);

                                        for k = 1: length(TAD_id_N2{iii,g})
                                            Chr_N2(n).x(TAD_id_N2{iii,g}(k,1)) = Trace_N2{iii,g}(k,1)*110/1000; % convert to um from pixel units
                                            Chr_N2(n).y(TAD_id_N2{iii,g}(k,1)) = sz*110/1000-Trace_N2{iii,g}(k,2)*110/1000;
                                            Chr_N2(n).z(TAD_id_N2{iii,g}(k,1)) = Trace_N2{iii,g}(k,3)*0.2;
                                            Chr_N2(n).r(TAD_id_N2{iii,g}(k,1)) = 1;
%                                             
                                        end
                                        Chr_N2(n).RoG = rog(Trace_N2{iii,g});
                                        Chr_N2(n).Age = q;
                                        Chr_N2(n).TraceNo = g+2;
                                        Chr_N2(n).TracesInTerritory = 4;
%                                        
                                        Chr_N2(n).nuc = nuc;
                                        Chr_N2(n).TraceName = [TraceName '_nuc' num2str(nuc) '_Terr' num2str(iii)];
                                    end
                                end
                             end   
                             
                             if HI_trace == 2  & N2_trace == 1
                                 
                                t = [];              

                                for i = 1:size(Trace_HI,2)

                                    if ~isempty(Trace_HI{iii,i})
                                        t = [t, i];
                                    end
                                end

                               
                                for iv = 1:length(t)
                                    g = t(iv);
                                
                                    if length(TAD_id_HI{iii,g})>1 && length(TAD_id_HI{iii,g}) == length(unique(TAD_id_HI{iii,g}))
                                        h = h+1;
                                        Chr_HI(h).x = zeros(TotalTADNum,1);
                                        Chr_HI(h).y = zeros(TotalTADNum,1);
                                        Chr_HI(h).z = zeros(TotalTADNum,1);
                                        Chr_HI(h).r = zeros(TotalTADNum,1);

                                        for k = 1: length(TAD_id_HI{iii,g})
                                            Chr_HI(h).x(TAD_id_HI{iii,g}(k,1)) = Trace_HI{iii,g}(k,1)*110/1000; % convert to um from pixel units
                                            Chr_HI(h).y(TAD_id_HI{iii,g}(k,1)) = sz*110/1000-Trace_HI{iii,g}(k,2)*110/1000;
                                            Chr_HI(h).z(TAD_id_HI{iii,g}(k,1)) = Trace_HI{iii,g}(k,3)*0.2;
                                            Chr_HI(h).r(TAD_id_HI{iii,g}(k,1)) = 1;
%                                             
                                        end
                                    
                                        Chr_HI(h).RoG = rog(Trace_HI{iii,g});
                                        Chr_HI(h).Age = q;
                                        Chr_HI(h).TraceNo = g;
                                        Chr_HI(h).TracesInTerritory = 3;
%                                         
                                        Chr_HI(h).nuc = nuc;
                                        Chr_HI(h).TraceName = [TraceName '_nuc' num2str(nuc) '_Terr' num2str(iii)];
                                    end
                                end
                                 
                                   
                                t = [];              

                                for i = 1:size(Trace_N2,2)

                                    if ~isempty(Trace_N2{iii,i})
                                        t = [t, i];
                                    end
                                end
                                
                                g = t; 
                                    
                                    if length(TAD_id_N2{iii,g})>1 && length(TAD_id_N2{iii,g}) == length(unique(TAD_id_N2{iii,g}))
                                        n = n+1;
                                        Chr_N2(n).x = zeros(TotalTADNum,1);
                                        Chr_N2(n).y = zeros(TotalTADNum,1);
                                        Chr_N2(n).z = zeros(TotalTADNum,1);
                                        Chr_N2(n).r = zeros(TotalTADNum,1);

                                        for k = 1: length(TAD_id_N2{iii,g})
                                            Chr_N2(n).x(TAD_id_N2{iii,g}(k,1)) = Trace_N2{iii,g}(k,1)*110/1000; % convert to um from pixel units
                                            Chr_N2(n).y(TAD_id_N2{iii,g}(k,1)) = sz*110/1000-Trace_N2{iii,g}(k,2)*110/1000;
                                            Chr_N2(n).z(TAD_id_N2{iii,g}(k,1)) = Trace_N2{iii,g}(k,3)*0.2;
                                            Chr_N2(n).r(TAD_id_N2{iii,g}(k,1)) = 1;
%                                             
                                        end
                                        Chr_N2(n).RoG = rog(Trace_N2{iii,g});
                                        Chr_N2(n).Age = q;
                                        Chr_N2(n).TraceNo = g+2;
                                        Chr_N2(n).TracesInTerritory = 3;
%                                         
                                        Chr_N2(n).nuc = nuc;
                                        Chr_N2(n).TraceName = [TraceName '_nuc' num2str(nuc) '_Terr' num2str(iii)];
                                    end
                             end
                                
                             if HI_trace == 1  & N2_trace == 2
                                 
                                 t = [];              

                                for i = 1:size(Trace_HI,2)

                                    if ~isempty(Trace_HI{iii,i})
                                        t = [t, i];
                                    end
                                end
                                
                                g = t; 
                                
                                    if length(TAD_id_HI{iii,g})>1 && length(TAD_id_HI{iii,g}) == length(unique(TAD_id_HI{iii,g}))
                                        h = h+1;
                                        Chr_HI(h).x = zeros(TotalTADNum,1);
                                        Chr_HI(h).y = zeros(TotalTADNum,1);
                                        Chr_HI(h).z = zeros(TotalTADNum,1);
                                        Chr_HI(h).r = zeros(TotalTADNum,1);

                                        for k = 1: length(TAD_id_HI{iii,g})
                                            Chr_HI(h).x(TAD_id_HI{iii,g}(k,1)) = Trace_HI{iii,g}(k,1)*110/1000; % convert to um from pixel units
                                            Chr_HI(h).y(TAD_id_HI{iii,g}(k,1)) = sz*110/1000-Trace_HI{iii,g}(k,2)*110/1000;
                                            Chr_HI(h).z(TAD_id_HI{iii,g}(k,1)) = Trace_HI{iii,g}(k,3)*0.2;
                                            Chr_HI(h).r(TAD_id_HI{iii,g}(k,1)) = 1;
%                                             
                                        end
                                    
                                        Chr_HI(h).RoG = rog(Trace_HI{iii,g});
                                        Chr_HI(h).Age = q;
                                        Chr_HI(h).TraceNo = g;
                                        Chr_HI(h).TracesInTerritory = 3;
%                                         
                                        Chr_HI(h).nuc = nuc;
                                        Chr_HI(h).TraceName = [TraceName '_nuc' num2str(nuc) '_Terr' num2str(iii)];
                                    end
                                 
                                    
                                     t = [];              

                                    for i = 1:size(Trace_N2,2)

                                        if ~isempty(Trace_N2{iii,i})
                                            t = [t, i];
                                        end
                                    end


                                    for iv = 1:length(t)
                                        g = t(iv);

                                    
                                        if length(TAD_id_N2{iii,g})>1 && length(TAD_id_N2{iii,g}) == length(unique(TAD_id_N2{iii,g}))
                                            n = n+1;
                                            Chr_N2(n).x = zeros(TotalTADNum,1);
                                            Chr_N2(n).y = zeros(TotalTADNum,1);
                                            Chr_N2(n).z = zeros(TotalTADNum,1);
                                            Chr_N2(n).r = zeros(TotalTADNum,1);

                                            for k = 1: length(TAD_id_N2{iii,g})
                                                Chr_N2(n).x(TAD_id_N2{iii,g}(k,1)) = Trace_N2{iii,g}(k,1)*110/1000; % convert to um from pixel units
                                                Chr_N2(n).y(TAD_id_N2{iii,g}(k,1)) = sz*110/1000-Trace_N2{iii,g}(k,2)*110/1000;
                                                Chr_N2(n).z(TAD_id_N2{iii,g}(k,1)) = Trace_N2{iii,g}(k,3)*0.2;
                                                Chr_N2(n).r(TAD_id_N2{iii,g}(k,1)) = 1;
%                                                
                                            end
                                            Chr_N2(n).RoG = rog(Trace_N2{iii,g});
                                            Chr_N2(n).Age = q;
                                            Chr_N2(n).TraceNo = g+1;
                                            Chr_N2(n).TracesInTerritory = 3;
%                                             
                                            Chr_N2(n).nuc = nuc;
                                            Chr_N2(n).TraceName = [TraceName '_nuc' num2str(nuc) '_Terr' num2str(iii)];
                                        end
                                    end
                                end
                             
                            end
                   
                        end
                          
                     flag = 0;
                    end
            end
        end
    end
 end
end 


toc

%%
Mean_HI = zeros(TotalTADNum,TotalTADNum);
Std_HI = zeros(TotalTADNum,TotalTADNum);
Median_HI = zeros(TotalTADNum,TotalTADNum);

Mean_N2 = zeros(TotalTADNum,TotalTADNum);
Std_N2 = zeros(TotalTADNum,TotalTADNum);
Median_N2 = zeros(TotalTADNum,TotalTADNum);


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


%%
save([AgeRange '\MeanPairDistanceMatrix_HI' num2str(SegChannel) '.mat'], 'Mean_HI');
save([AgeRange '\MedianPairDistanceMatrix_HI' num2str(SegChannel) '.mat'], 'Median_HI');
save([AgeRange '\StdPairDistanceMatrix_HI' num2str(SegChannel) '.mat'], 'Std_HI');
save([AgeRange '\SEMPairDistanceMatrix_HI' num2str(SegChannel) '.mat'], 'SEM_HI');
save([AgeRange '\NofDataPairDistanceMatrix_HI' num2str(SegChannel) '.mat'], 'NofData_HI');
save([AgeRange '\AllChromosomes_HI' num2str(SegChannel) '.mat'], 'Chr_HI', '-v7.3');
save([AgeRange '\DisListAll_HI' num2str(SegChannel) '.mat'], 'DisListAll_HI');


save([AgeRange '\MeanPairDistanceMatrix_N2' num2str(SegChannel) '.mat'], 'Mean_N2');
save([AgeRange '\MedianPairDistanceMatrix_N2' num2str(SegChannel) '.mat'], 'Median_N2');
save([AgeRange '\StdPairDistanceMatrix_N2' num2str(SegChannel) '.mat'], 'Std_N2');
save([AgeRange '\SEMPairDistanceMatrix_N2' num2str(SegChannel) '.mat'], 'SEM_N2');
save([AgeRange '\NofDataPairDistanceMatrix_N2' num2str(SegChannel) '.mat'], 'NofData_N2');
save([AgeRange '\AllChromosomes_N2' num2str(SegChannel) '.mat'], 'Chr_N2', '-v7.3');
save([AgeRange '\DisListAll_N2' num2str(SegChannel) '.mat'], 'DisListAll_N2');


