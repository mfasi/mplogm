function [handler] = plotBarStackGroups_bck_fwd(stackData, groupLabels)
%% Plot a set of stacked bars, but group them according to labels provided.
%%
%% Params: 
%%      stackData is a 3D matrix (i.e., stackData(i, j, k) => (Group, Stack, StackElement)) 
%%      groupLabels is a CELL type (i.e., { 'a', 1 , 20, 'because' };)
%%
%% Copyright 2011 Evan Bollig (bollig at scs DOT fsu ANOTHERDOT edu
%%
%% 
NumGroupsPerAxis = size(stackData, 1);
NumStacksPerGroup = size(stackData, 2);


% Count off the number of bins
groupBins = 1:NumGroupsPerAxis;
MaxGroupWidth = 0.45; % Fraction of 1. If 1, then we have all bars in groups touching
groupOffset = MaxGroupWidth/NumStacksPerGroup;
hold on; 
for i=1:NumStacksPerGroup

    Y = squeeze(stackData(:,i,:));
    
    % Center the bars:
    
    internalPosCount = i - ((NumStacksPerGroup+1) / 2);
    
    % Offset the group draw positions:
    groupDrawPos = (internalPosCount)* groupOffset + groupBins;
    
    h(i,:) = bar(Y, 'stacked');
    colormap(jet);
    C = colormap();
    colours = C(1:10:end, :);
    colours = colours([1,6,2,5,3,4],:);
    if (i == 1)
        for k = 1:3
            h(i,k).FaceColor = colours(k, :);
            legend(h(1,:),'s','s','m')
        end
    elseif (i == 2)
        for k = 1:3
            h(i,k).FaceColor = colours(k+3, :);
            legend(h(2,:),'s','s','m')
        end
        legend([h(1,1), h(1,2),...
            h(2,1), h(2,2)],...
            '$s$', 'm','$s$', '$m$',...
            'Location', 'NE');
    end
    set(h(i,:),'BarWidth',groupOffset);
    set(h(i,:),'XData',groupDrawPos);
end
handler = h(NumStacksPerGroup,:);
hold off;
set(gca,'XTickMode','manual');
set(gca,'XTick',1:NumGroupsPerAxis);
set(gca,'XTickLabelMode','manual');
set(gca,'XTickLabel',groupLabels);
end