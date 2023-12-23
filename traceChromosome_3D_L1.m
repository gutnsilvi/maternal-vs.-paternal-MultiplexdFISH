function [Trace, TAD_id] = traceChromosome_3D_L1(L1,w, Xfit_new, Yfit_new, Zfit, Intensity)
% Trace is a cell array containing 2-4 trace coordinates; TAD_id is a cell
% array containing the TAD ID number of the foci in the traces
global flap 
global x1
global x2
global y1
global y3

NumImgRounds = length(Intensity); % number of imaging rounds.

% cycle through all images to find foci that are within the imaged volume
for i = 1:NumImgRounds
    Index = [];
    for j = 1:length(Xfit_new{i})
        if Zfit{i}(j)>1 && Zfit{i}(j)<151 && Yfit_new{i}(j)>1 && Yfit_new{i}(j)< ceil(y3-y1) ...
                && Xfit_new{i}(j)>1 && Xfit_new{i}(j)< ceil(x2-x1)
            if L1(round(Yfit_new{i}(j)), round(Xfit_new{i}(j)), round(Zfit{i}(j))) == w
                Index = [Index j];
            end
        end
    end
    
    foci{i}.x = Xfit_new{i}(Index);
    foci{i}.y = Yfit_new{i}(Index);
    foci{i}.z = Zfit{i}(Index);
    foci{i}.intensity = Intensity{i}(Index);
    NumFoci(i) = length(Index); % number of foci in this imaging round

% add condition if more than 4 foci detected (i.e. if there is some background signal), use only top 4 foci in intensity
% for tracing (sort descend foci using intensity indices)
    if NumFoci(1, i) > 4
        [sortfoci{i}.intensity,ind] = sort(foci{i}.intensity, 'descend');
        sortfoci{i}.x = foci{i}.x(ind);
        sortfoci{i}.y = foci{i}.y(ind);
        sortfoci{i}.z = foci{i}.z(ind);

        foci{1,i}.x = sortfoci{1,i}.x(1,1:4);
        foci{1,i}.y = sortfoci{1,i}.y(1,1:4);
        foci{1,i}.z = sortfoci{1,i}.z(1,1:4);
        foci{1,i}.intensity = sortfoci{1,i}.intensity(1,1:4);
        NumFoci(1,i) = length(foci{1,i}.intensity);
    end
end


% determine how many traces are there in the regions defined by L and
% I: see if the majority of images have 1, 2, 3 or 4 foci in this region
N1 = length(find(NumFoci == 1));
N2 = length(find(NumFoci == 2));
N3 = length(find(NumFoci == 3));
N4 = length(find(NumFoci == 4));

if N4==N3 && N4==N2 && N4==N1 && N4==0
    NumTraces = 0;
elseif N4>=N3 && N4>=N2 && N4>=N1
    NumTraces = 4;
elseif N3>=N4 && N3>=N2 && N3>=N1
    NumTraces = 3;
elseif N2>=N3 && N2>=N4 && N2>=N1
    NumTraces = 2;
elseif N1>=N3 && N1>=N4 && N1>=N2
    NumTraces = 1;
else
    NumTraces = 0;
end

% if number of traces is 0, stop the tracing and exit; otherwise find the
% first image that show exactly the same number of foci as the identified number of traces.
% start the tracing from that imaging round.


if NumTraces == 0
    Trace = nan;
    TAD_id = nan;
    display(['No traces in regions ' num2str(w)]);
    flap = 1;
else
    flap = 0;
    % find the starting round for tracing analysis
    Idx = find(NumFoci == NumTraces);
    StartingRound = Idx(1);
    % start the forward tracing
    for i = 1:NumTraces
        Trace{i} = [foci{StartingRound}.x(i), foci{StartingRound}.y(i), foci{StartingRound}.z(i)]; %xyz coordinates of first foci in each trace
        TAD_id{i} = StartingRound; %tad id of each coordinate in corresponding trace cell
        Merged(i) = 0; % this is for keeping track of which trace has merged with which trace
    end
    MergedOrNot = 0; % 0 means no merging event
    for i = StartingRound:(NumImgRounds-1) % i is the current round
        NextRound = i+1;
        if isempty(foci{NextRound}.x)
        end
        % when no merging event has happened, try to mantain the number of
        % traces, do not allow split, but allow pairs of traces to merge,
        % do not allow 3 or 4 traces to merge into one
        if MergedOrNot == 0
            for j = 1:NumTraces
                if ~isempty(foci{NextRound}.x)
                    for k = 1:length(foci{NextRound}.x)
                        disList(k) = dis(Trace{j}(end,:),[foci{NextRound}.x(k) foci{NextRound}.y(k) foci{NextRound}.z(k)]);
                    end
                    [M, Idx] = min(disList); %not min(disList(k))
                    clear disList
                    Trace{j} = [Trace{j}; [foci{NextRound}.x(Idx) foci{NextRound}.y(Idx) foci{NextRound}.z(Idx)]];
                    TAD_id{j} = [TAD_id{j}; NextRound]; %changed from i to j
                end
            end
            % determine if traces have merged
            flag = 0;
            for j = 1:NumTraces-1
                for k = j+1:NumTraces
                    if dis(Trace{j}(end,:),Trace{k}(end,:)) == 0
                        MergedOrNot = 1;
                        if Merged(j) == 0 && Merged(k) == 0 && max(Merged) == 0
                            Merged(j) = 1;
                            Merged(k) = 1;
                        elseif Merged(j) == 0 && Merged(k) == 0 && max(Merged) == 1
                            Merged(j) = 2;
                            Merged(k) = 2;
                        else
                            flag = 1;
                        end
                    end
                end
            end
            if flag == 1 %more than 2 traces are merged into one, throw away
                Merged = zeros(NumTraces,1);
                MergedOrNot = 0;
                for j = 1:NumTraces
                    Trace{j} = Trace{j}(1:end-1,:);
                    TAD_id{j} = TAD_id{j}(1:end-1,:);
                end
            end
        else
            % when merging event has happened, try to split the traces to go back to the orginial
            % number of traces
            if length(foci{NextRound}.x)>0
                if length(foci{NextRound}.x) > NumTraces
                    [B, Idx] = sort(foci{NextRound}.intensity,'descend');
                    foci{NextRound}.x = foci{NextRound}.x(Idx(1:NumTraces));
                    foci{NextRound}.y = foci{NextRound}.y(Idx(1:NumTraces));
                    foci{NextRound}.z = foci{NextRound}.z(Idx(1:NumTraces));
                end
                for j = 1:length(foci{NextRound}.x)
                    for k = 1:NumTraces
                        disList = [];
                        for m = 1:size(Trace{k},1)
                            disList = [disList dis(Trace{k}(m,:),[foci{NextRound}.x(j) foci{NextRound}.y(j) foci{NextRound}.z(j)])];
                        end
                        disToTrace(k) = mean(disList);
                        clear disList
                    end
                    [M, Idx] = min(disToTrace);
                    DestinyTrace(j) = Idx;
                    clear disToTrace
                end
                for j = 1:length(foci{NextRound}.x)
                    Trace{DestinyTrace(j)} = [Trace{DestinyTrace(j)}; [foci{NextRound}.x(j) foci{NextRound}.y(j) foci{NextRound}.z(j)]];
                    TAD_id{DestinyTrace(j)} = [TAD_id{DestinyTrace(j)}; NextRound];
                end
            end
            % determine if traces have merged
            MergedOrNot = 0;
            Merged = zeros(NumTraces,1);
            for j = 1:NumTraces-1
                for k = j+1:NumTraces
                    if dis(Trace{j}(end,:),Trace{k}(end,:)) == 0
                        MergedOrNot = 1;
                        if Merged(j) == 0 && Merged(k) == 0 && max(Merged) == 0
                            Merged(j) = 1;
                            Merged(k) = 1;
                        elseif Merged(j) == 0 && Merged(k) == 0 && max(Merged) == 1
                            Merged(j) = 2;
                            Merged(k) = 2;
                        end
                    end
                end
            end
        end
    end
end
