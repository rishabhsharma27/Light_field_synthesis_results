%EPI_analysis_test6:(Previous version) 
%1. this code generates the horizontal and vertical EPI with their depth maps
%2. this code also is working with priority i.e the slope from -4 to 0 are
%   filled first and then the code runs for slope 1 to 4 to maintain the occlusion
%3. The code is still not able to resolve where the slope between two adjoining pixels 
%   changes from -ve slope to +ve slope and vis versa

%EPI_analysis_test7(This version): 
%1. every slope is done separately from -4 to 4 making sure the occlusion
%is handeld 
%
%problem statement:
%


function [a1,a2,img_h_view1,img_h_depth] = EPI_analysis_test12_freq_shift(LF,depth_output,img_dir,Slope,Shift2,var_fill)
%img_dir is for direction of EPI 
% 1. horizontal EPI
% 2. veritcal EPI
% 3. depth map EPI
tic


% Slope_1=(40:-1:0);
% num_pixel_per_slope_1= [45,36,36,36,36,36,36,36,36,27,27,27,27,27,27,27,27,27,18,18,18,18,18,18,18,18,18,9,9,9,9,9,8,7,6,5,5,4,3,3,1];
% 
% Slope_2=0:1:40;
% num_pixel_per_slope= [1,3,3,4,5,5,6,7,8,9,9,9,9,9,18,18,18,18,18,18,18,18,18,27,27,27,27,27,27,27,27,27,36,36,36,36,36,36,36,36,45];

Slope_1=(40:-1:0);
num_pixel_per_slope_1= [33,31,31,31,29,29,27,27,27,25,25,23,23,23,21,21,19,19,19,17,17,15,15,15,15,13,11,11,11,9,9,9,7,7,5,5,3,3,3,3,1];

Slope_2=0:1:40;
num_pixel_per_slope= [1,3,3,3,3,5,5,7,7,9,9,9,11,11,11,13,15,15,15,15,17,17,19,19,19,21,21,23,23,23,25,25,27,27,27,29,29,31,31,31,33];
%% old pattern files 
% % % Slope_1=(35:-1:0);
% % % num_pixel_per_slope_1=[29,29,27,26,25,25,24,23,22,22,21,20,19,18,17,17,17,16,15,14,13,12,11,10,10,9,9,8,7,6,5,5,4,3,3,1];
% % % 
% % % 
% % % Slope_2=0:1:35;
% % % num_pixel_per_slope=[1,3,3,4,5,5,6,7,8,9,9,10,10,11,12,13,14,15,16,17,17,17,18,19,20,21,22,22,23,24,25,25,26,27,29,29];
%%

depth=round(depth_output.*-10,1);
   
depth(1:10,:) =500;
depth(503:512,:) =500;
depth(:,1:10) =500;
depth(:,503:512) =500;

mask= depth==500;

depth = round(regionfill(depth,mask));

Slope1=round(Slope.*10);
C = unique(depth);
for j=1:size(C,1)
    [slope_cor,s_cor_ind] = min(abs(Slope1 - C(j,1)));
    slope_cor_val(j) = Slope1(s_cor_ind);
end

    a=0;
    a = round(slope_cor_val - C',1);
    for m=1:size(a,2)
        if (a(1,m)~=0)
            cor_var = 0;
            cor_var = C(m);
            for i=1:size(depth,1)
                for j=1:size(depth,2)
                    if(depth(i,j)==cor_var)
                        depth(i,j) = slope_cor_val(m);
                    else
                        depth(i,j) = depth(i,j);
                    end
                end
            end
        end
    end

clear slope_cor s_cor_ind  slope_cor_val cor_var a C ;

% depth=depth_output;

% figure,imagesc(depth),colormap(gray),pause(1)
num_sub=9;

if (size(LF,4)==1)
  size_temp = [1,2];
else 
  size_temp = [3,4];
end

if(img_dir==1)
    c_view=LFDisp(LF);
    c_view=im2double(c_view);
    img_h_view=(NaN(num_sub*size(LF,size_temp(1)),size(LF,size_temp(2)),3));
    img_h_depth=double(NaN(num_sub*size(LF,size_temp(1)),size(LF,size_temp(2))));
    double img_depth;
elseif (img_dir==2)
    img_h_view=(NaN(num_sub*size(LF,size_temp(2)),size(LF,size_temp(1)),3));
    depth = imrotate(depth,-90); 
    c_view=LFDisp(LF);
    c_view=im2double(c_view);
    c_view = imrotate(c_view,-90);
    img_h_depth=double(NaN(num_sub*size(LF,size_temp(2)),size(LF,size_temp(1))));
    double img_depth;
 end
count=1;

%     img_h_view1=img_h_view;
%     img_h_depth1=img_h_depth;
    
% figure;imagesc(img_h_view),colormap(gray)

%%

count=1;
for j=1:9
    for i=j:9:size(img_h_view,1)
        a1{j}(count,:,:)= img_h_view(i,:,:);
        a2_d{j}(count,:)= img_h_depth(i,:);
        count=count+1;
    end
    count=1;
end


C = unique(depth);
count = 1; 
count1 = 1; 
C_neg = [];
C_pos = [];
for i = 1:size(C,1)
    if (C(i,1)<0)
        C_neg(count,1) = C(i,1);
        count = count + 1;
    elseif(C(i,1)>0)
        C_pos(count1,1) = C(i,1);
        count1 = count1 + 1;
    end
    
end

count=1;
for j=1:9
    for i=j:9:size(img_h_view,1)
        a1{j}(count,:,:)= img_h_view(i,:,:);
        a2_d{j}(count,:)= img_h_depth(i,:);
        count=count+1;
    end
    count=1;
end

w = hann(512);
w=w.*w';


a_depth=depth;

[a,b,~]=size(a_depth);
Nr = ifftshift([-fix(a/2):ceil(a/2)-1]);
Nc = ifftshift([-fix(b/2):ceil(b/2)-1]);
[xF,yF] = meshgrid(Nr,Nc);
tic
pattern_shift = -40:1:0;

if (~isempty(C_neg))
    for m = 1:size(C_neg,1)
        for n = 1:size(pattern_shift,2)
            if (C_neg(m,1)==pattern_shift(1,n))
                mask_depth=depth==pattern_shift(1,n);
                in_d=(ifftshift(fftshift(fft2(depth.*mask_depth)).*w));

                a_rgb=c_view;
                
                in_rgb(:,:,1)=(ifftshift(fftshift(fft2(a_rgb(:,:,1)))));
                in_rgb(:,:,2)=(ifftshift(fftshift(fft2(a_rgb(:,:,2)))));
                in_rgb(:,:,3)=(ifftshift(fftshift(fft2(a_rgb(:,:,3)))));
                
                for k=1:4
                    x0=pattern_shift(1,n)*(5-k)/(10); %// Define shifts
                    y0=0;

                    H=exp(1i*2*pi*(x0*xF/a+y0*yF/b));

                    IF_d = round(real(ifft2(in_d.*H)));

                    IF_rgb(:,:,1) = round(real(ifft2(in_rgb(:,:,1).*H)),4);
                    IF_rgb(:,:,2) = round(real(ifft2(in_rgb(:,:,2).*H)),4);
                    IF_rgb(:,:,3) = round(real(ifft2(in_rgb(:,:,3).*H)),4);

                    x0=round(x0);

                    if (x0>0)
                       IF_rgb(:,end-abs(x0)+1:end,:)=0;
                       IF_d(:,end-abs(x0)+1:end)=0;
                    elseif(x0<0)
                        IF_rgb(:,1:abs(x0),:)=0;
                        IF_d(:,1:abs(x0),:)=0;
                    end
                    mask1 = (IF_d<(C_neg(m,1)-10) | IF_d<(C_neg(m,1)+10));
                    mask3 = IF_d~=0;
                    IF_d = C_neg(m,1).*(mask1.*mask3);
                    mask2 = IF_d==C_neg(m,1);
                    IF_rgb = IF_rgb.*mask2;
                    for i=1:size(a1{5},1)
                        for j=1:size(a1{5},2)
                            if(isnan(a2_d{k}(i,j)) && (IF_d(i,j)~=0))
                                a2_d{k}(i,j) = IF_d(i,j);
                                a1{k}(i,j,:) = IF_rgb(i,j,:);
                            else
                                a2_d{k}(i,j) = a2_d{k}(i,j);
                                a1{k}(i,j,:) = a1{k}(i,j,:);
                            end
                        end
                    end

    %                 figure,imagesc(IF_rgb),pause(1)
    %                 figure,imagesc(IF_d),pause(1)
                end


                for k=6:9
                    x0=-pattern_shift(1,n)*(k-5)/(10); %// Define shifts
                    y0=0;

                    H=exp(1i*2*pi*(x0*xF/a+y0*yF/b));

                    IF_d = round(real(ifft2(in_d.*H)));

                    IF_rgb(:,:,1) = round(real(ifft2(in_rgb(:,:,1).*H)),4);
                    IF_rgb(:,:,2) = round(real(ifft2(in_rgb(:,:,2).*H)),4);
                    IF_rgb(:,:,3) = round(real(ifft2(in_rgb(:,:,3).*H)),4);

                    x0=round(x0);

                    if (x0>0)
                       IF_rgb(:,end-abs(x0)+1:end,:)=0;
                       IF_d(:,end-abs(x0)+1:end)=0;
                    elseif(x0<0)
                        IF_rgb(:,1:abs(x0),:)=0;
                        IF_d(:,1:abs(x0),:)=0;
                    end
                    mask1 = (IF_d<(C_neg(m,1)-10) | IF_d<(C_neg(m,1)+10));
                    mask3 = IF_d~=0;
                    IF_d = C_neg(m,1).*(mask1.*mask3);
    %                 IF_d = C_neg(m,1).*mask1;
                    mask2 = IF_d==C_neg(m,1);
                    IF_rgb = IF_rgb.*mask2;
                    for i=1:size(a1{5},1)
                        for j=1:size(a1{5},2)
                            if(isnan(a2_d{k}(i,j)) && (IF_d(i,j)~=0))
                                a2_d{k}(i,j) = IF_d(i,j);
                                a1{k}(i,j,:) = IF_rgb(i,j,:);
                            else
                                a2_d{k}(i,j) = a2_d{k}(i,j);
                                a1{k}(i,j,:) = a1{k}(i,j,:);
                            end
                        end
                    end

    %                 figure,imagesc(IF_rgb),pause(1)
    %                 figure,imagesc(IF_d),pause(1)
                end
            end
        end
    end
end

for k=1:9
    for i=1:size(a1{5},1)
        for j=1:size(a1{5},2)
            if(isnan(a2_d{k}(i,j)) && (depth(i,j)==0))
                a2_d{k}(i,j) = 0;
                a1{k}(i,j,:) = c_view(i,j,:);
            else
                a2_d{k}(i,j) = a2_d{k}(i,j);
                a1{k}(i,j,:) = a1{k}(i,j,:);
            end
        end
    end
end

pattern_shift = 0:1:40;

if (~isempty(C_pos))
    for m = 1:size(C_pos,1)
        for n = 1:size(pattern_shift,2)
            if (C_pos(m,1)==pattern_shift(1,n))
                mask_depth=depth==pattern_shift(1,n);
                in_d=(ifftshift(fftshift(fft2(depth.*mask_depth)).*w));

                a_rgb=c_view;
                
                in_rgb(:,:,1)=(ifftshift(fftshift(fft2(a_rgb(:,:,1)))));
                in_rgb(:,:,2)=(ifftshift(fftshift(fft2(a_rgb(:,:,2)))));
                in_rgb(:,:,3)=(ifftshift(fftshift(fft2(a_rgb(:,:,3)))));
                
                for k=1:4
                    x0=pattern_shift(1,n)*(5-k)/(10); %// Define shifts
                    y0=0;

                    H=exp(1i*2*pi*(x0*xF/a+y0*yF/b));

                    IF_d = round(real(ifft2(in_d.*H)));

                    IF_rgb(:,:,1) = round(real(ifft2(in_rgb(:,:,1).*H)),4);
                    IF_rgb(:,:,2) = round(real(ifft2(in_rgb(:,:,2).*H)),4);
                    IF_rgb(:,:,3) = round(real(ifft2(in_rgb(:,:,3).*H)),4);

                    x0=round(x0);

                    if (x0>0)
                       IF_rgb(:,end-abs(x0)+1:end,:)=0;
                       IF_d(:,end-abs(x0)+1:end)=0;
                    elseif(x0<0)
                        IF_rgb(:,1:abs(x0),:)=0;
                        IF_d(:,1:abs(x0),:)=0;
                    end
                    mask1 = (IF_d<(C_pos(m,1)-10) | IF_d<(C_pos(m,1)+10));
                    mask3 = IF_d~=0;
                    IF_d = C_pos(m,1).*(mask1.*mask3);
                    mask2 = IF_d==C_pos(m,1);
                    IF_rgb = IF_rgb.*mask2;
                    for i=1:size(a1{5},1)
                        for j=1:size(a1{5},2)
                            if(isnan(a2_d{k}(i,j)) && (IF_d(i,j)~=0))
                                a2_d{k}(i,j) = IF_d(i,j);
                                a1{k}(i,j,:) = IF_rgb(i,j,:);
                            else
                                a2_d{k}(i,j) = a2_d{k}(i,j);
                                a1{k}(i,j,:) = a1{k}(i,j,:);
                            end
                        end
                    end

    %                 figure,imagesc(IF_rgb),pause(1)
    %                 figure,imagesc(IF_d),pause(1)
                end


                for k=6:9
                    x0=-pattern_shift(1,n)*(k-5)/(10); %// Define shifts
                    y0=0;

                    H=exp(1i*2*pi*(x0*xF/a+y0*yF/b));

                    IF_d = round(real(ifft2(in_d.*H)));

                    IF_rgb(:,:,1) = round(real(ifft2(in_rgb(:,:,1).*H)),4);
                    IF_rgb(:,:,2) = round(real(ifft2(in_rgb(:,:,2).*H)),4);
                    IF_rgb(:,:,3) = round(real(ifft2(in_rgb(:,:,3).*H)),4);

                    x0=round(x0);

                    if (x0>0)
                       IF_rgb(:,end-abs(x0)+1:end,:)=0;
                       IF_d(:,end-abs(x0)+1:end)=0;
                    elseif(x0<0)
                        IF_rgb(:,1:abs(x0),:)=0;
                        IF_d(:,1:abs(x0),:)=0;
                    end
                    mask1 = (IF_d<(C_pos(m,1)-10) | IF_d<(C_pos(m,1)+10));
                    mask3 = IF_d~=0;
                    IF_d = C_pos(m,1).*(mask1.*mask3);
                    mask2 = IF_d==C_pos(m,1);
                    IF_rgb = IF_rgb.*mask2;
                    for i=1:size(a1{5},1)
                        for j=1:size(a1{5},2)
                            if(isnan(a2_d{k}(i,j)) && (IF_d(i,j)~=0))
                                a2_d{k}(i,j) = IF_d(i,j);
                                a1{k}(i,j,:) = IF_rgb(i,j,:);
                            else
                                a2_d{k}(i,j) = a2_d{k}(i,j);
                                a1{k}(i,j,:) = a1{k}(i,j,:);
                            end
                        end
                    end

    %                 figure,imagesc(IF_rgb),pause(1)
    %                 figure,imagesc(IF_d),pause(1)
                end
            end
        end
    end
end
a1{5} = c_view;
a2_d{5} = depth;
toc



%% depth map and RGB holes filling

[a2] = Depth_fill_test_2022(a2_d);

for k=1:9
     for i=1:size(a2{1},1)
        for j=1:size(a2{1},2)
            if(isnan(a2_d{k}(i,j)) && a2{k}(i,j)==500)
                a2_d_1{k}(i,j) = 500;
            else
                a2_d_1{k}(i,j) = a2_d{k}(i,j);
            end
        end
     end
 end

[a2_1] = Depth_fill(a2_d);

for i=1:9
    mask = a2{i}==500;
    temp_d = a2_1{i}.*mask;
    a2{i} = a2{i}.*~mask;
    a2{i} = a2{i} + temp_d;
    a2_d_1{i} = a2_d_1{i}.*~mask;
    a2_d_1{i} = a2_d_1{i} + temp_d;
end




count=1;
for j=1:size(a1{1},1)
    for i=1:9
        img_h_depth(count,:)= a2_d{i}(j,:);
        count=count+1;
    end
end

[a2_temp] = Depth_fill_test_1(img_h_depth);

 for k=1:9
     for i=1:size(a2{1},1)
        for j=1:size(a2{1},2)
            if(isnan(a2{k}(i,j)))
                a2{k}(i,j) = a2_temp{k}(i,j);
            else
                a2{k}(i,j) = a2{k}(i,j);
            end
        end
     end
 end
 


for i=1:9
    a2_d_temp{i} = a2_d_1{i}.*-1;
end

[a2_opp] = Depth_fill_test_2022(a2_d_temp);

  for i=1:9
    a2_opp{i} = a2_opp{i}.*-1;
  end
  for k=1:9
     for i=1:size(a2{1},1)
        for j=1:size(a2{1},2)
            if(isnan(a2_opp{k}(i,j)))
                a2_opp{k}(i,j) = a2_temp{k}(i,j);
            else
                a2_opp{k}(i,j) = a2_opp{k}(i,j);
            end
        end
     end
  end
 


% [a2]=Depth_fill_test_2(a2_d);




 for k=1:9
    for i=1:size(a2{1},1)
        for j=1:size(a1{1},2)
            if(isnan(a1{k}(i,j,1)))
                a1{k}(i,j,:)=0; 
            else
                a1{k}(i,j,:)=a1{k}(i,j,:);
            end
        end
    end
 end

        

 
  


if (img_dir==1)
     for j=1:9
        a2{j}= (a2{j}(:,:))/10; 
        a2_opp{j}= (a2_opp{j}(:,:))/10;
%         a2{j}=(medfilt2(a2{j}, [ 3 3 ]));
%         a2{j}=round(a2{j});
%         mask=isnan(a2{j}); 
%         a2{j}=(regionfill(a2{j},mask));
     end
elseif (img_dir==2)
    for j=1:9
        a2{j}= (a2{j}(:,:))/10; 
        a1{j}= imrotate(a1{j},180);    
        a2{j}= imrotate(a2{j},180);
        a2_opp{j}= imrotate(a2_opp{j},180);
        a2_opp{j}= (a2_opp{j}(:,:))/10;
%         a2{j}=(medfilt2(a2{j}, [ 3 3 ]));
%         a2{j}=round(a2{j});
%         mask=isnan(a2{j}); 
%         a2{j}=(regionfill(a2{j},mask));
     end
end

if (img_dir==2)
    for j=1:9
        a4{j}= (a1{10 - j}); 
        a4_depth{j}= (a2{10 - j});
        a5_depth{j}= (a2_opp{10 - j});
     end
end

if (img_dir==2)
    for j=1:9
        a1{j}= a4{j};
        a2{j}= a4_depth{j};
        a2_opp{j}= a5_depth{j};
    end
end

clear C;
for i=1:9
    C{i} = unique(a2{i});
    for j=1:size(C{i},1)
        slope_cor=0;
        s_cor_ind=0;
        [slope_cor,s_cor_ind] = min(abs(Slope - C{i}(j,1)));
        slope_cor_val{i}(j) = Slope(s_cor_ind);
    end
end

for k=1:9
    a=0;
    a = round(slope_cor_val{k} - C{k}',1);
    for m=1:size(a,2)
        if (a(1,m)~=0)
            cor_var = 0;
            cor_var = round(C{k}(m).*10,1);
            for i=1:size(a2{1},1)
                for j=1:size(a2{1},2)
                    if(round(a2{k}(i,j).*10,1)==cor_var)
                        a2{k}(i,j) = slope_cor_val{k}(m);
                    else
                        a2{k}(i,j) = a2{k}(i,j);
                    end
                end
            end
        end
    end
end

clear C;
clear slope_cor_val;
for i=1:9
    C{i} = unique(a2_opp{i});
    for j=1:size(C{i},1)
        slope_cor=0;
        s_cor_ind=0;
        [slope_cor,s_cor_ind] = min(abs(Slope - C{i}(j,1)));
        slope_cor_val{i}(j) = Slope(s_cor_ind);
    end
end


for k=1:9
    a=0;
    a = round(slope_cor_val{k} - C{k}',1);
    for m=1:size(a,2)
        if (a(1,m)~=0)
            cor_var = 0;
            cor_var = round(C{k}(m).*10,1);
            for i=1:size(a2_opp{1},1)
                for j=1:size(a2_opp{1},2)
                    if(round(a2_opp{k}(i,j).*10,1)==cor_var)
                        a2_opp{k}(i,j) = slope_cor_val{k}(m);
                    else
                        a2_opp{k}(i,j) = a2_opp{k}(i,j);
                    end
                end
            end
        end
    end
end






a3=a1;
a4=a2;
if (img_dir==1)
    [a1] = RGB_fill_horizontal_test_2022(a3,a2,Slope,Shift2,a2_opp,var_fill) ;
%     [a1] = RGB_fill_horizontal_test_2022_less_fs(a3,a2,Slope,Shift2,a2_opp,var_fill) ;
elseif(img_dir==2)
    for j=1:size(Shift2,2)
        Shift2{j}= imrotate(Shift2{j},90);     
    end
    [a1] = RGB_fill_horizontal_test_2022(a3,a2,Slope,Shift2,a2_opp,var_fill) ;
%     [a1] = RGB_fill_horizontal_test_2022_less_fs(a3,a2,Slope,Shift2,a2_opp,var_fill) ;
    for j=1:9
        a1{j}= imrotate(a1{j},-90);    
        a2{j}= imrotate(a2{j},-90); 
        a2_opp{j}= imrotate(a2_opp{j},-90); 
    end
    
end



% for j=1:9
%     figure(j),imagesc(a1{j}),pause(1)
% end
img_h_view1=0;
% 
% for j=1:9
%     figure(j),imagesc(a2{j}),colormap(gray),pause(1)
% end
toc
clear d;

end


