%% this script sorts out distances between homologs which are larger than a threshold


TotalTADNum = 22;

load('2to40cell/homologs/AllChromosomes561_homologs.mat', 'Chr_h1', 'Chr_h2')

DisList_chr = cell(length(Chr_h1),1);
for k = 1 :length(Chr_h1) 
   DisList_chr{k,2} =  sum(Chr_h1(k).r);
   DisList_chr{k,3} =  sum(Chr_h2(k).r);
   Dis_chr = zeros(TotalTADNum,TotalTADNum);
    for i = 1:TotalTADNum
        for j = 1:TotalTADNum

            if Chr_h1(k).r(i) == 1 && Chr_h2(k).r(j) == 1  
                Dis_chr(i,j) = ((Chr_h1(k).x(i)-Chr_h2(k).x(j))^2+(Chr_h1(k).y(i)-Chr_h2(k).y(j))^2+(Chr_h1(k).z(i)-Chr_h2(k).z(j))^2)^0.5;
                DisList_chr{k,1} = Dis_chr;
            end
        end
    end
   DisList_chr{k,4} =    max(max(DisList_chr{k,1}));
   DisList_chr{k,5} = Chr_h1(k).TraceName;
   DisList_chr{k,6} = Chr_h2(k).TraceName;
end

%% creates a list for all homolog comparissons which are longer then 6µm (this means the nuclei were not well separated)

sortOut = [];
cut_off = 6;

for i = 1:length(DisList_chr)

    if DisList_chr{i,4} > cut_off
        sortOut = [sortOut  i];
    end
end



Chr_h1(sortOut) = [];
Chr_h2(sortOut) = [];


%save('2to40cell/homologs/AllChromosomes561_homologs_filtered.mat', 'Chr_HI', 'Chr_N2')
                

%%
