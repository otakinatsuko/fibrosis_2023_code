% clear all
% close all

path = './input/';
files = dir([path,'*.tif']);
mag = 1/30;
savepath ='./output_v04/';

for kk = 1:length(files)
    %% read image
    img = imread(strcat(path,files(kk).name));
    img = double(img(50:3250,50:3250,:));
    [x y z] = size(img);

    %% circular area detection and segmentation
    img_gray = double(rgb2gray(uint8(img)));
    img_gray_lowres = imresize(single(img_gray), mag);
    [gx,gy] = gradient(img_gray_lowres);
    grad = sqrt(gx.^2+gy.^2);    
    grad_mean = mean(grad,'all');
    grad_sd = std(grad,0,'all');  
    [centers, radii, metric] = imfindcircles(grad>(grad_mean + 2*grad_sd),[40 70],Sensitivity=0.95);    
    [xg yg] = size(img_gray);
    [xl yl] = size(img_gray_lowres);   
    center = (centers(1,:)-xl/2)*30 + xg/2;
    radius = radii(1)*30*0.7;
    [xt,yt]=meshgrid(1:x,1:y);    
    mask_circle = sqrt((xt-center(1,1)).^2+(yt-center(1,2)).^2)<=radius;
    img_central = img.*mask_circle;
    % figure;imagesc(uint8(img));

    %% illumination normalization
    ill = imgaussfilt(movmean(movmean(img_gray,100,1),100,2),20);
    ill = ill / max(max(ill.*mask_circle));
    img_spatial_normalized = img./ill.*mask_circle;
    % figure;imagesc(ill)
    % figure;imagesc(uint8(img_spatial_normalized));

    %% normalize each color with background
    mask_nonwhite = (((mean(img_spatial_normalized,3))./std(img_spatial_normalized,[],3))<100);
    img_bg = img_spatial_normalized.*(1-mask_nonwhite);
    % figure;imagesc(uint8(img_bg))
    img_bg(isnan(img_bg))=0;
    for k=1:3
        temp_img = img_bg(:,:,k);
        bg(k) = mean(temp_img(temp_img>0));
        img_bg_normalized(:,:,k) = img_spatial_normalized(:,:,k)/bg(k);
    end
    img_bg_normalized = img_bg_normalized*max(bg);
    % figure;imagesc(uint8(img_bg_normalized))
    
    %% create cell mask
    img_absorption = (max(bg) - img_bg_normalized).*mask_circle;
    [x,y,z] = size(img_absorption);
    input = zeros(x*y,z);
    for k=1:y
      input(x*(k-1)+1:x*k,:) = img_absorption(1:x,k,:);
    end
    [coeff, score] = pca(input);
      output = zeros(x,y,z);
    for k=1:y
      output(1:x,k,:) = score(x*(k-1)+1:x*k,:);
    end
    mask_cells = output(:,:,1)>0;
    % figure(201);imagesc(output(:,:,1));daspect([1 1 1])


    %% create output images
    for k=1:3
        img_cells_abs(:,:,k) = output(:,:,1)*coeff(k,1);
    end
    img_noncells_abs = -img_cells_abs.*(1-mask_cells);
    img_cells_abs = img_cells_abs.*mask_cells;
    img_cells_trans = img_bg_normalized.*mask_cells;
    img_noncells_trans = img_bg_normalized.*(1-mask_cells);
    figure(110);
    subplot(1,3,1);imagesc(uint8(img_cells_trans));daspect([1 1 1])
    subplot(1,3,2);imagesc(uint8(img_noncells_trans));daspect([1 1 1])
    subplot(1,3,3);imagesc(uint8(img_cells_abs*5));daspect([1 1 1])
    drawnow
    % waitforbuttonpress

    %% metrics
    input = img_cells_trans;
    intR = sum(sum(input(:,:,1)))/sum(sum(input(:,:,1)>0));
    intG = sum(sum(input(:,:,2)))/sum(sum(input(:,:,2)>0));
    intB = sum(sum(input(:,:,3)))/sum(sum(input(:,:,3)>0));
    input = img_cells_abs/max(bg);
    absR = sum(sum(input(:,:,1)))/sum(sum(input(:,:,1)>0));
    absG = sum(sum(input(:,:,2)))/sum(sum(input(:,:,2)>0));
    absB = sum(sum(input(:,:,3)))/sum(sum(input(:,:,3)>0));
    area_ratio = sum(sum(squeeze(input(:,:,1)>0)))/sum(mask_circle(:));
    output_metrics(:,kk) = [intR,intG,intB,absR,absG,absB,area_ratio];

    imwrite(uint8(img_central), strcat(savepath, files(kk).name,'_central_','.tiff'))
    imwrite(uint8(img_bg_normalized), strcat(savepath, files(kk).name,'_norm_','.tiff'))
    imwrite(uint8(img_cells_trans), strcat(savepath, files(kk).name,'_cells_trans_','.tiff'))
    imwrite(uint8(img_noncells_trans), strcat(savepath, files(kk).name,'_nocells_trans_','.tiff'))
    imwrite(uint8(img_absorption/max(img_absorption(:))*255*5), strcat(savepath, files(kk).name,'_cells_absorption_','.tiff'))
    imwrite(uint8(img_cells_abs/max(img_cells_abs(:))*255*5), strcat(savepath, files(kk).name,'_cells_abs_','.tiff'))
    imwrite(uint8(img_noncells_abs/max(img_cells_abs(:))*255*5), strcat(savepath, files(kk).name,'_noncells_abs_','.tiff'))
    imwrite(uint8(mask_cells*255), strcat(savepath, files(kk).name,'_cellmask','.tiff'))
    figure(200);plot(coeff(:,1),'ko-');xlim([0.5,3.5]);ylim([0,1]);saveas(gcf,strcat(savepath, files(kk).name,'_pca'),'epsc')
end
% fasdfa
for k=1:7
    strength = [squeeze(output_metrics(k,:))];
    alloy = {'5','5','5','5','3','3','3','3','4','4','4','4','1','1','1','1','2','2','2','2'};
    [~,~,stats] = anova1(strength,alloy);
end

%% export table
T = array2table(output_metrics);
T.Properties.VariableNames = ["5_1","5_2","5_3","5_4","3_1","3_2","3_3","3_4","4_1","4_2","4_3","4_4","1_1","1_2","1_3","1_4","2_1","2_2","2-3","2_4"];
T.Properties.RowNames = ["intR","intG","intB","absR","absG","absB","area_ratio"];

writetable(T, strcat(savepath, "sirius_red.txt"),'Delimiter','\t','WriteRowNames',true)