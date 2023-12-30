% This script quantifies the overlap of N2 and HI territory signals 
% watershed segmentation of nuclei and primary probe signal
% ------------------------------------------------------------------------------------------------------------------------------------------
% Input: overlaps 
% 
%
% ------------------------------------------------------------------------------------------------------------------------------------------
% Output: compiled txt file of all overlap volumes


tic;
clear all
close all


Vol_all = [];

Folder{1} = 'xx';
Folder{2} = 'xx';
Folder{3} = 'xx';
Folder{4} = 'xx';

for x = 1:length(Folder)
    a=dir([Folder{x} '/*.mat']);
    n = size(a,1);

    for ii = 1:n
        load([Folder{x} '/' a(ii).name])
        AgeIndex = strfind(a(ii).name,'cell');
        age = str2num(a(ii).name(1:AgeIndex-1));
        for i = 1:size(Vol,2)
            Vol(i).age = age;
            Vol(i).ratio = Vol(i).overlab/(Vol(i).N2 + Vol(i).Hi);
            Vol(i).name = [Folder{x} '/' a(ii).name];
    
        end
        Vol_all = [Vol_all, Vol];
    
    end
end

b = struct2table(Vol_all);
writetable(b,'overlap.txt')
        