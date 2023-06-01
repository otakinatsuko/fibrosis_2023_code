clear all
close all

%% Read file path
path = pwd;
files = dir([path,'/*_ch01.tif']);

%% Output matrix
T_whole = table({'IL33_pos_cell_num';'mean_intensity'});
T_whole.Properties.VariableNames{1} = 'variable name';

%% 
for N = 1:length(files)
%for N = 15:15
    
    %% Read each file
    temp_file = files(N).name;
    img = imread(strcat(path,'/', temp_file));
    
    %figure;imagesc(img);
    %size(img)

    %% Set output colname
    temp_col_name = strsplit(temp_file,'_');
    
    %% Set threshold of fluorescnece intensity
    % set threshtold of the intensity
    th = 30;
    I = img(:,:,2) > th;
    
    %figure;imagesc(I);

    %% Dilution, removal of small objects and errosion
    
    se90 = strel('line',3,90);
    se0 = strel('line',3,0);
    
    I_2 = imdilate(I,[se90 se0]);  
    %figure;imagesc(I_2);
    %figure;imagesc(~I_2);
    
    I_3 = imerode(I_2,[se90 se0],'full');
    %figure;imagesc(I_3);
    
    % set threshold of the size of the small objects
    th2 = 20;
    I_4 = bwareaopen(I_3, th2);
    %figure;imagesc(I_4);
    
    %% Count the number of the cells
    
    s = regionprops(I_4(2:513,2:513),img(:,:,2),'centroid','MeanIntensity');
    %s = regionprops(I_4(2:513,2:513),img(:,:,2),'centroid');
    centroids = cat(1,s.Centroid);
    %mean(cat(1,s.MeanIntensity))

    ob_n = size(centroids);
    ob_n = ob_n(1);

    %% write images
    imwrite(I_4(2:513, 2:513), strcat('output/',temp_col_name{5},'_ch01_processed','.tiff'))
    
    %% Write tables
    T_temp = [ob_n;mean(cat(1,s.MeanIntensity))];
    T_temp = table(T_temp);
    T_temp.Properties.VariableNames{1} = temp_col_name{5};
    T_whole = horzcat(T_whole, T_temp);

end

%% export table
writetable(T_whole,strcat('output/', 'out_ch01.xlsx'))