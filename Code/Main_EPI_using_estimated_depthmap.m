%% info

% Variable 1 (LF)
% size(LF) =   9     9   512   512     3 
%
%
% Variable 2 (depth_output)
% ground truth depth map 
% slope for -4:0.1:4
%
%
% Variable 3 (img_dir) is for direction of EPI 
% 1. horizontal EPI
% 2. veritcal EPI
% 3. depth map EPI 
N1=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28];
% N1=6;
for m=1
    
    N=N1(m);
    
%% synthetic images
%load(['C:\MATLAB\Add-Ons\lytro_data\Depth-Estimation-Light-Field-master\LF\DepthEstimation(Ch4)\dataset\images for testing\synthetic LF image matlab files\' num2str(N) '.mat'] , 'LF');


%% sythetic images, but using estimated depth not groundtruth

%load(['C:\MATLAB\Add-Ons\lytro_data\Depth-Estimation-Light-Field-master\LF\DepthEstimation(Ch4)\dataset\images for testing\selected images\all_synthetic_image_comparing_Gtruth_medfilt_13\' num2str(N) '.mat'], 'FinalImg2');

[LF,depth_Our,Slope,FinalImg1_reconstructed,Shift2]=Our_depth_code(N); %synthetic light field image number
% depth_output=depth_Our;
% depth_output=zeros(size(depth_Our,1),size(depth_Our,2));

a = parsePfm( 'gt_disp_lowres.pfm');
depth_output=im2single(round(flip(a,1),1));

%Slope=-4:0.2:4;

% for k=1:size(Slope,2)
%     for i=1:size(depth_output,1)
%         for j=1:size(depth_output,2)
%             if(depth_Our(i,j)==k)
%                 depth_output(i,j)=-Slope(1,k);
%                 %the '_ve sign is put so the 
%                 %code 'EPI_analysis_test8' doesnt have to be altered.
%             end
%         end
%     end
% end

% p1 = 0.89534;
% p2 = -2.313e-16;


    
%corr_slope1(1,:) = round(p1*Slope(1,:) + p2,1); 
    
% corr_slope=[-3.6,-3.5,-3.4,-3.3,-3.2,-3.2,-3.1,-3,-2.9,-2.8,-2.7,-2.6,-2.5,-2.4,-2.3,-2.2,-2.2,-2.1,-2,-1.9,-1.8,-1.7,-1.6,-1.5,-1.4,-1.3,-1.2,-1.1,-1,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.4,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.4,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,2.1,2.2,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3,3.1,3.2,3.2,3.3,3.4,3.5,3.6];
% 
% slope_1=-4:0.1:4;
% 
% corr_slope=corr_slope.*10;
% slope_1=round(slope_1*10);
% 
% depth_Our_1=round(depth_Our.*10);
% 
% for k=1:size(slope_1,2)
%     for i=1:size(depth_Our,1)
%         for j=1:size(depth_Our,2)
%             if ((depth_Our_1(i,j))==corr_slope(1,k))
%                 depth_Our_2(i,j)= slope_1(1,k)/10;
%                 Slope1(1,k) = slope_1(1,k)/10;
%             end
%         end
%     end
% end

depth_output=depth_Our.*-1;
figure;imagesc(depth_output),colormap(gray),pause(1)

tstart1 = tic();

[a1,a2,img_h_view,img_h_depth] = EPI_analysis_test12_freq_shift(LF,depth_output,1,Slope,Shift2,0); 

[a3,a4,img_v_view,img_v_depth] = EPI_analysis_test12_freq_shift(LF,depth_output,2,Slope,Shift2,0);

for i=1:9
    a2{i}=a2{i}.*-1; 
end

for i=1:9
    a4{i}=a4{i}.*-1; 
end

%% Other views synthetic images
%horizontal

    [a1_1,a1_1_depth,img_h_view_1,img_h_depth_1] = EPI_analysis_test12_freq_shift(a3{1},a4{1},1,Slope,Shift2,1); 
    
    [a1_2,a1_2_depth,img_h_view_2,img_h_depth_2] = EPI_analysis_test12_freq_shift(a3{2},a4{2},1,Slope,Shift2,2); 
    
    [a1_3,a1_3_depth,img_h_view_3,img_h_depth_3] = EPI_analysis_test12_freq_shift(a3{3},a4{3},1,Slope,Shift2,3); 
    
    [a1_4,a1_4_depth,img_h_view_4,img_h_depth_4] = EPI_analysis_test12_freq_shift(a3{4},a4{4},1,Slope,Shift2,4); 
    
    [a1_6,a1_6_depth,img_h_view_6,img_h_depth_6] = EPI_analysis_test12_freq_shift(a3{6},a4{6},1,Slope,Shift2,-4); 
    
    [a1_7,a1_7_depth,img_h_view_7,img_h_depth_7] = EPI_analysis_test12_freq_shift(a3{7},a4{7},1,Slope,Shift2,-3); 
    
    [a1_8,a1_8_depth,img_h_view_8,img_h_depth_8] = EPI_analysis_test12_freq_shift(a3{8},a4{8},1,Slope,Shift2,-2); 
    
    [a1_9,a1_9_depth,img_h_view_9,img_h_depth_9] = EPI_analysis_test12_freq_shift(a3{9},a4{9},1,Slope,Shift2,-1); 
    
    
    
%vertical    
    [b1_1,b1_1_depth,img_v_view_1,img_v_depth_1] = EPI_analysis_test12_freq_shift(a1{1},a2{1},2,Slope,Shift2,1); 
    
    [b1_2,b1_2_depth,img_v_view_2,img_v_depth_2] = EPI_analysis_test12_freq_shift(a1{2},a2{2},2,Slope,Shift2,2); 
    
    [b1_3,b1_3_depth,img_v_view_3,img_v_depth_3] = EPI_analysis_test12_freq_shift(a1{3},a2{3},2,Slope,Shift2,3); 
    
    [b1_4,b1_4_depth,img_v_view_4,img_v_depth_4] = EPI_analysis_test12_freq_shift(a1{4},a2{4},2,Slope,Shift2,4); 
    
    [b1_6,b1_6_depth,img_v_view_6,img_v_depth_6] = EPI_analysis_test12_freq_shift(a1{6},a2{6},2,Slope,Shift2,-4); 
    
    [b1_7,b1_7_depth,img_v_view_7,img_v_depth_7] = EPI_analysis_test12_freq_shift(a1{7},a2{7},2,Slope,Shift2,-3); 
    
    [b1_8,b1_8_depth,img_v_view_8,img_v_depth_8] = EPI_analysis_test12_freq_shift(a1{8},a2{8},2,Slope,Shift2,-2); 
    
    [b1_9,b1_9_depth,img_v_view_9,img_v_depth_9] = EPI_analysis_test12_freq_shift(a1{9},a2{9},2,Slope,Shift2,-1); 

clear img_v_view_1 img_v_view_2 img_v_view_3 img_v_view_4 img_v_view_5 img_v_view_6 img_v_view_7 img_v_view_8 img_v_view_9;

clear img_v_depth_1 img_v_depth_2 img_v_depth_3 img_v_depth_4 img_v_depth_5 img_v_depth_6 img_v_depth_7 img_v_depth_8 img_v_depth_9;

clear img_h_view_1 img_h_view_2 img_h_view_3 img_h_view_4 img_h_view_5 img_h_view_6 img_h_view_7 img_h_view_8 img_h_view_9;

clear img_h_depth_1 img_h_depth_2 img_h_depth_3 img_h_depth_4 img_h_depth_5 img_h_depth_6 img_h_depth_7 img_h_depth_8 img_h_depth_9;

telapsed1 = toc(tstart1);

LF1 = im2double(zeros(9,9,size(a3{1},1),size(a3{1},2),3));
for i=1:9
    for j=1:9
        if(i==1) LF1(i,j,1:size(a3{1},1),1:size(a3{1},2),1:3)= a1_1{j}; end
        if(i==2) LF1(i,j,1:size(a3{1},1),1:size(a3{1},2),1:3)= a1_2{j}; end
        if(i==3) LF1(i,j,1:size(a3{1},1),1:size(a3{1},2),1:3)= a1_3{j}; end
        if(i==4) LF1(i,j,1:size(a3{1},1),1:size(a3{1},2),1:3)= a1_4{j}; end
        if(i==5) LF1(i,j,1:size(a3{1},1),1:size(a3{1},2),1:3)= a1{j};   end
        if(i==6) LF1(i,j,1:size(a3{1},1),1:size(a3{1},2),1:3)= a1_6{j}; end
        if(i==7) LF1(i,j,1:size(a3{1},1),1:size(a3{1},2),1:3)= a1_7{j}; end
        if(i==8) LF1(i,j,1:size(a3{1},1),1:size(a3{1},2),1:3)= a1_8{j}; end
        if(i==9) LF1(i,j,1:size(a3{1},1),1:size(a3{1},2),1:3)= a1_9{j}; end
    end
end

LF1_depth = im2double(zeros(9,9,size(a3{1},1),size(a3{1},2)));
for i=1:9
    for j=1:9
        if(i==1) LF1_depth(i,j,1:size(a3{1},1),1:size(a3{1},2))= a1_1_depth{j}; end
        if(i==2) LF1_depth(i,j,1:size(a3{1},1),1:size(a3{1},2))= a1_2_depth{j}; end
        if(i==3) LF1_depth(i,j,1:size(a3{1},1),1:size(a3{1},2))= a1_3_depth{j}; end
        if(i==4) LF1_depth(i,j,1:size(a3{1},1),1:size(a3{1},2))= a1_4_depth{j}; end
        if(i==5) LF1_depth(i,j,1:size(a3{1},1),1:size(a3{1},2))= a2{j};   end
        if(i==6) LF1_depth(i,j,1:size(a3{1},1),1:size(a3{1},2))= a1_6_depth{j}; end
        if(i==7) LF1_depth(i,j,1:size(a3{1},1),1:size(a3{1},2))= a1_7_depth{j}; end
        if(i==8) LF1_depth(i,j,1:size(a3{1},1),1:size(a3{1},2))= a1_8_depth{j}; end
        if(i==9) LF1_depth(i,j,1:size(a3{1},1),1:size(a3{1},2))= a1_9_depth{j}; end
    end
end


LF2 = im2double(zeros(9,9,size(a1{1},1),size(a1{1},2),3));
for i=1:9
    for j=1:9
        if(i==1) LF2(j,i,1:size(a1{1},1),1:size(a1{1},2),1:3)= b1_1{j}; end
        if(i==2) LF2(j,i,1:size(a1{1},1),1:size(a1{1},2),1:3)= b1_2{j}; end
        if(i==3) LF2(j,i,1:size(a1{1},1),1:size(a1{1},2),1:3)= b1_3{j}; end
        if(i==4) LF2(j,i,1:size(a1{1},1),1:size(a1{1},2),1:3)= b1_4{j}; end
        if(i==5) LF2(j,i,1:size(a1{1},1),1:size(a1{1},2),1:3)= a3{j};   end
        if(i==6) LF2(j,i,1:size(a1{1},1),1:size(a1{1},2),1:3)= b1_6{j}; end
        if(i==7) LF2(j,i,1:size(a1{1},1),1:size(a1{1},2),1:3)= b1_7{j}; end
        if(i==8) LF2(j,i,1:size(a1{1},1),1:size(a1{1},2),1:3)= b1_8{j}; end
        if(i==9) LF2(j,i,1:size(a1{1},1),1:size(a1{1},2),1:3)= b1_9{j}; end
    end
end

LF2_depth = im2double(zeros(9,9,size(a3{1},1),size(a3{1},2)));
for i=1:9
    for j=1:9
        if(i==1) LF2_depth(j,i,1:size(a3{1},1),1:size(a3{1},2))= b1_1_depth{j}; end
        if(i==2) LF2_depth(j,i,1:size(a3{1},1),1:size(a3{1},2))= b1_2_depth{j}; end
        if(i==3) LF2_depth(j,i,1:size(a3{1},1),1:size(a3{1},2))= b1_3_depth{j}; end
        if(i==4) LF2_depth(j,i,1:size(a3{1},1),1:size(a3{1},2))= b1_4_depth{j}; end
        if(i==5) LF2_depth(j,i,1:size(a3{1},1),1:size(a3{1},2))= a4{j};   end
        if(i==6) LF2_depth(j,i,1:size(a3{1},1),1:size(a3{1},2))= b1_6_depth{j}; end
        if(i==7) LF2_depth(j,i,1:size(a3{1},1),1:size(a3{1},2))= b1_7_depth{j}; end
        if(i==8) LF2_depth(j,i,1:size(a3{1},1),1:size(a3{1},2))= b1_8_depth{j}; end
        if(i==9) LF2_depth(j,i,1:size(a3{1},1),1:size(a3{1},2))= b1_9_depth{j}; end
    end
end



for i=1:9
    for j=1:9
        [s_LF1(i,j),s8{i,j}] = ssim(im2double(LFDisp(LF(i,j,:,:,:))),im2double(LFDisp(LF1(i,j,:,:,:))));
    end
end

for i=1:9
    for j=1:9
        [s_LF2(i,j),s8{i,j}] = ssim(im2double(LFDisp(LF(i,j,:,:,:))),im2double(LFDisp(LF2(i,j,:,:,:))));
    end
end

LF_avg = im2uint8(zeros(9,9,size(a1{1},1),size(a1{1},2),3));
for i=1:9
    for j=1:9
        LF_avg(i,j,:,:,:) = im2uint8((LF1(i,j,:,:,:) + LF2(i,j,:,:,:))/2);
        [s_LF_avg(i,j),s8{i,j}] = ssim(rgb2gray(im2uint8(LFDisp(LF(i,j,16:512-15,16:512-15,:)))),rgb2gray(im2uint8(LFDisp(LF_avg(i,j,16:512-15,16:512-15,:)))));
    end
end

s_LF_avg_note = 'measuring accuracy by removing 15 pixels off the edge on 4 sides ssim(im2double(LFDisp(LF(i,j,15:512-15,15:512-15,:))),im2double(LFDisp(LF_avg(i,j,15:512-15,15:512-15,:))))';

% for i=1:9
%     for j=1:9
%     LF_temp1= LFDisp(LF1(i,j,:,:,:));
%     LF_temp2= LFDisp(LF2(i,j,:,:,:));
%     LF1_temp(i,j,:,:,:)=imsharpen(imgaussfilt((LF_temp1),1));
%     LF2_temp(i,j,:,:,:)=imsharpen(imgaussfilt((LF_temp2),1));
%     end
% end
% 
% Slope=-4:0.2:4;
% 
% 
%     f1=@(Slope) LFFiltShiftSum(LF1,Slope);
%     f2=@(Slope) LFFiltShiftSum(LF2,Slope);
% for i = 1:size(Slope,2)
%     ShiftImg1 = feval(f1,Slope(i));
%     ShiftImg2 = feval(f2,Slope(i));
%     if (size(ShiftImg1,3) == 4)
%     ShiftImg1 = squeeze(ShiftImg1(:,:,1:end-1));
%     ShiftImg2 = squeeze(ShiftImg2(:,:,1:end-1));
%     end
%     Shift1{i}=im2double(ShiftImg1);
%     Shift2{i}=im2double(ShiftImg2);
%     figure,imshow(Shift1{i}),pause(1)
% end


save(['C:\MATLAB\Add-Ons\lytro_data\Depth-Estimation-Light-Field-master\LF\DepthEstimation(Ch4)\dataset\images for testing\EPI_reconstruction\synthetic_LF_result\' num2str(N)],...
        'depth_Our', 'depth_output', 'FinalImg1_reconstructed' , 'LF1', 'LF2', 'LF1_depth' ,...
        'LF2_depth', 'Shift2', 'Slope', 's_LF1', 's_LF2', 'telapsed1', 's_LF_avg', 'LF_avg');

close all;


end