clear all
close all

path = '/Users/natsuko/Desktop/sirius_red_20230530/';
files = dir([path,'*.tif']);
P=50;
mag = 1/30; % image magnification
threshold = 1.00;

for kk = 1:length(files)
    img = imread(strcat(path,files(kk).name));
    img = double(img(50:3250,50:3250,:)); % crop the image to the same size
    [x y z] = size(img);
    
    %% circular area detection and segmentation
    img_gray = rgb2gray(uint8(img));
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

    %% rough detection of non-white region
    mask_nonwhite = (((mean(img_central,3))./std(img_central,[],3))<30);
    mask_nonwhite = bwareaopen(mask_nonwhite,P);
    img_nonwhite = img_central.*mask_nonwhite;

    %% PCA and color matrix retrieval
    [x,y,z] = size(img_nonwhite);
    input = zeros(x*y,z);
    for k=1:y
      input(x*(k-1)+1:x*k,:) = img_nonwhite(1:x,k,:);
    end
    [coeff, score] = pca(input);
      output = zeros(x,y,z);
    for k=1:y
      output(1:x,k,:) = score(x*(k-1)+1:x*k,:);
    end
    mat = [coeff(1,2), 0, 0; 0, coeff(2,2), 0; 0,0,coeff(3,2)];
    invmat =inv(mat);
    % figure;imagesc(output(:,:,2))
    % figure;plot(coeff(:,2))

    %% illumination variation estimation and normalization
    temp3 = invmat(1,1)*double(img(:,:,1))+invmat(2,2)*double(img(:,:,2))+invmat(3,3)*double(img(:,:,3));
    % figure;imagesc(-temp3)
    ill = imgaussfilt(movmean(movmean(temp3,100,1),100,2),20);
    ill = ill .* mask_circle;
    if mean2(ill)<0
        ill = -ill;
        temp3 = -temp3;
    end
    temp = temp3./ill.*mask_circle;
    img_norm = img_central./ill*mean(ill(ill>0));
    mask_cells = (temp<threshold);
    img_cells = img_norm.*mask_cells;
    img_cells(isnan(img_cells))=0;
    
    %% mask to deep red region
    mask_deep_red = (temp>0.8);
    
    % figure;imagesc(uint8(img_central));daspect([1 1 1]);drawnow
    % % bg = movmean(movmean(img./ill*mean2(ill),100,1),100,2);
    % % figure;imagesc(uint8(img./ill*mean2(ill)-bg+mean2(bg)))
    % figure;imagesc(uint8(img_cells));daspect([1 1 1]);drawnow
    % figure;imagesc(uint8(img_cells.*mask_deep_red));daspect([1 1 1]);drawnow
    figure(7);
    subplot(1,3,1);imagesc(uint8(img_central));daspect([1 1 1]);drawnow
    subplot(1,3,2);imagesc(uint8(img_cells));daspect([1 1 1]);drawnow
    subplot(1,3,3);imagesc(uint8(img_cells.*mask_deep_red));daspect([1 1 1]);drawnow

    %% metrics
    input = img_cells;
    intR = sum(sum(input(:,:,1)))/sum(sum(input(:,:,1)>0));
    intG = sum(sum(input(:,:,2)))/sum(sum(input(:,:,2)>0));
    intB = sum(sum(input(:,:,3)))/sum(sum(input(:,:,3)>0));
    area_ratio = sum(sum(sum(input>0)))/sum(mask_circle(:));
    output_metrics(:,kk) = [intR,intG,intB,area_ratio];

    input = img_cells.*mask_deep_red;
    intR = sum(sum(input(:,:,1)))/sum(sum(input(:,:,1)>0));
    intG = sum(sum(input(:,:,2)))/sum(sum(input(:,:,2)>0));
    intB = sum(sum(input(:,:,3)))/sum(sum(input(:,:,3)>0));
    area_ratio = sum(sum(sum(input>0)))/sum(mask_circle(:));
    output_metrics_wo_deep_red(:,kk) = [intR,intG,intB,area_ratio];
   
end

for k=1:4
    strength = [squeeze(output_metrics(k,:))];
    alloy = {'5','5','5','5','3','3','3','3','4','4','4','4','1','1','1','1','2','2','2','2'};
    [~,~,stats] = anova1(strength,alloy);
end
for k=1:4
    strength = [squeeze(output_metrics_wo_deep_red(k,:))];
    alloy = {'5','5','5','5','3','3','3','3','4','4','4','4','1','1','1','1','2','2','2','2'};
    [~,~,stats] = anova1(strength,alloy);
end