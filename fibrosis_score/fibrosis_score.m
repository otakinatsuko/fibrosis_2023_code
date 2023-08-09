clear all
close all

%% read files
addpath('MyFunctions')
% path = '/Users/natsuko/Desktop/input_example/';
path = './input/';
files = dir([path,'/*.tif']);
s = 10;
roi=100;
th = 0.5;
th_final = 5;

for K = 1:length(files)
    %% read each image
    temp_file = files(K).name;
    img = double(imread(strcat(path,'/', temp_file)));
    img = img(:,:,1:3);
    [x, y, z] = size(img);

    %% fill edges with mean white value
    temp_nonwhite = (((mean(img,3))./std(img,[],3))<30);
    mymean = mean(img(~temp_nonwhite));
    for k=1:3
        img(1:s,:,k) = mymean;
        img(end-s+1:end,:,k) = mymean;
        img(:,1:s,k) = mymean;
        img(:,end-s+1:end,k) = mymean;
    end

    %% calculate nonwhite region
    th1 = 180;
    temp =mean(img,3);
    mask = (temp<th1);
    th2 = sum(sum(temp .* mask))/sum(sum(mask))*1.3;
    mask_nonwhite = (temp< th2);
    
    %% calculate area of tissue section
    mask_tissue = find_tissue_area(mask_nonwhite,s);
    figure(1);
    subplot(2,1,1);imagesc(uint8(img));daspect([1 1 1]);
    subplot(2,1,2);imagesc(uint8(img.*mask_tissue));daspect([1 1 1])
    drawnow

    %% calculate fibrotic region
    temp = movmean(movmean(mask_nonwhite,roi, 1),roi, 2);
    temp2 = temp(logical((temp>0.1).*(temp<0.4)));
    temp2= round(100*temp2)/100;
    th3 = mode(temp2)*1.6;
    mask_fibrotic_movmean = movmean(movmean(mask_nonwhite,roi, 1),roi, 2);
    mask_fibrotic = bwareaopen(mask_nonwhite.*(mask_fibrotic_movmean>th3),10000);
    img_fibrotic = img.*mask_fibrotic;

    figure(2)
    subplot(3,1,1);imagesc(uint8(img));daspect([1 1 1]);drawnow
    subplot(3,1,2);imagesc(uint8(img_fibrotic));daspect([1 1 1])
    subplot(3,1,3);imagesc(mask_fibrotic);daspect([1 1 1])
    drawnow
    % waitforbuttonpress

    %% color discrimination of fibrotic region
    clear input output
    input = img_fibrotic;
    [img_fibrotic_blue,img_fibrotic_red] = color_discrimination(input, th_final);
    figure(3);
    subplot(1,3,1);imagesc(uint8(input));daspect([1 1 1])
    subplot(1,3,2);imagesc(uint8(img_fibrotic_red));daspect([1 1 1])
    subplot(1,3,3);imagesc(uint8(img_fibrotic_blue));daspect([1 1 1])
    ratio(K) = sum(sum(img_fibrotic_blue(:,:,1)>0))/sum(sum(mask_tissue>0));

    %% save images
    imwrite(uint8(img_fibrotic_red), strcat('./output/',temp_file,'_pca_red_','.tiff'))
    imwrite(uint8(img_fibrotic_blue), strcat('./output/',temp_file,'_pca_blue_','.tiff')) 
    imwrite(uint8(mask_nonwhite*255), strcat('./output/',temp_file,'_mask_nonwhite_','.tiff'))
    imwrite(uint8(mask_tissue*255), strcat('./output/',temp_file,'_mask_tissue_','.tiff'))
    imwrite(uint8(mask_fibrotic_movmean/max(mask_fibrotic_movmean(:))*255), strcat('./output/',temp_file,'_mask_fibrotic_movmean_','.tiff'))
    imwrite(uint8(img_fibrotic), strcat('./output/',temp_file,'_img_fibrotic_','.tiff'))
    
    K
end

%%
writematrix(ratio*100,"./output/scores.txt","Delimiter","\t")
