function [out_blue,out_red] = color_discrimination(input, th_final)
    [x,y,z]=size(input);
    for k=1:y
        input2(x*(k-1)+1:x*k,:) = input(1:x,k,:);
    end
    [coeff, score] = pca(input2);
    output = zeros(x,y,z);
    for k=1:y
        output(1:x,k,:) = score(x*(k-1)+1:x*k,:);
    end

    %% color discrimination based on PCA. Since the second and third comonents of the PCA can be swapped depending on the input image, we need to determine which component to use for each input image.
    %generate color-discrimanted images with second component of PCA
    clear temp1 temp2 temp3 temp4
    th=0;
    mask = (output(:,:,2))>th;
    for k=1:z
        temp1(:,:,k) = input(:,:,k).*mask;
    end
    for k=1:z
        temp2(:,:,k) = input(:,:,k).*(1-mask);
    end

    %generate color-discrimanted images with third component of PCA
    mask = (output(:,:,3))>th;
    for k=1:z
        temp3(:,:,k) = input(:,:,k).*mask;
    end
    for k=1:z
        temp4(:,:,k) = input(:,:,k).*(1-mask);
    end
%     figure;
%     subplot(3,2,3);imagesc(uint8(temp1));daspect([1 1 1])
%     subplot(3,2,4);imagesc(uint8(temp2));daspect([1 1 1])
%     subplot(3,2,5);imagesc(uint8(temp3));daspect([1 1 1])
%     subplot(3,2,6);imagesc(uint8(temp4));daspect([1 1 1])

    %determine which component (second or third) to use. Smaller std in each masked image means color discrimination is going better. also, adjust the polarity of the PCA coefficient to determine which channel is red or blue.
    cr2 = abs(((sum(sum(std(temp2,[],3)))-sum(sum(std(temp1,[],3))))/((sum(sum(std(temp2,[],3)))+sum(sum(std(temp1,[],3)))))));
    cr3 = abs(((sum(sum(std(temp3,[],3)))-sum(sum(std(temp4,[],3))))/((sum(sum(std(temp3,[],3)))+sum(sum(std(temp4,[],3)))))));
    if cr2<cr3 %std of component 2 < std of component 3
        if coeff(1,2) > coeff(3,2) %component 2: red > blue channel
            mask = (output(:,:,2))>th_final;
        else %component 2: red < blue channel
            mask = -(output(:,:,2))>th_final;
        end
    else %std of component 3 < std of component 2
        if coeff(1,3) > coeff(3,3) %component 3: red > blue channel
            mask = (output(:,:,3))>th_final;
        else %component 3: red < blue channel
            mask = -(output(:,:,3))>th_final;
        end
    end
    for k=1:z
        out_red(:,:,k) = input(:,:,k).*mask;
    end
    for k=1:z
        out_blue(:,:,k) = input(:,:,k).*(1-mask);
    end
end