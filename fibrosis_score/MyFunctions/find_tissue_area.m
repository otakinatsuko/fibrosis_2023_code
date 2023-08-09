function [output] = find_tissue_area(input_mask,s)
    %% dilate the mask and find large contours
    se = strel('square',s*2);
    temp = imdilate(input_mask,se);
    % figure;imagesc(temp)
    [x,y] = size(temp);
    CC = bwconncomp(temp);
    numPixels = cellfun(@numel,CC.PixelIdxList);
    idx = find((numPixels>5*10^3)==1);
    temp2 = zeros(x,y);
    for l=1:length(idx)
        temp2(CC.PixelIdxList{idx(l)})=1;
    end
    % figure;imagesc(temp2);

    temp3 = imfill(temp2,'holes');
    % figure;imagesc(temp3);
    
    %% connect the edge regions
    s2 = 100;
    se2 = strel('line',s2,0);
    temp4 = temp3(1,:);
    idx1 = imerode(movmean(temp4,s2),se2)>0;

    temp4 = temp3(end,:);
    idx2 = imerode(movmean(temp4,s2),se2)>0;

    temp4 = temp3(:,1);
    idx3 = imerode(movmean(temp4,s2),se2)>0;

    temp4 = temp3(:,1);
    idx4 = imerode(movmean(temp4,s2),se2)>0;

    temp3(1,idx1)=1;
    temp3(end,idx2)=1;
    temp3(idx3,1)=1;
    temp3(idx4,1)=1;
    
    %% output
    temp3 = imfill(temp3,'holes');
    se = strel('square',s*2);
    output = imerode(temp3,se);
end