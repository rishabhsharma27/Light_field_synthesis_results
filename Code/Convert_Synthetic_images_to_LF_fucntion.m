function [LF]=Convert_Synthetic_images_to_LF_fucntion(i)

% if error check if input variable i 

%cd(['C:\MATLAB\Add-Ons\lytro_data\Depth-Estimation-Light-Field-master\LF\DepthEstimation(Ch4)\dataset\images for testing\synthetic images\' num2str(i)]);


cd(['G:\MATLAB\synthetic_LF_image1']);

%%



count=1;
count1=1;
% for N=1:11
%     a=imresize(imread(['input_Cam00' num2str(N) '.png']),0.5);
%     LF(1,N+1,:,:,:)=a;
% LF(1,N,:,:,:)=imread(['A_' num2str(N) '.png']);
% end

files = dir(fullfile('*input*'));


no_sub=sqrt(size(files,1));

for N=1:no_sub^2
%         LF(count,count1,:,:,:)=imread(['out_ (' num2str(N) ').png']);
    if ((N-1)<10)
        img=imresize(im2double(imread(['input_Cam00' num2str(N-1) '.png'])),1);
%         img_gy=rgb2gray(img);
%         img_1= flip(img_gy,2);
%         [Gx1,Gy1] = imgradientxy(img_gy,'central');
%         [Gx2,Gy2] = imgradientxy(img_1,'central');
%         
%         img_3=(Gx1+Gy1);
%         img_4=(Gx2+Gy2);
% 
%         img_4= flip(img_4,2);
%         
%         img=im2uint8(img + 2*img_3 + 2*img_4);
        LF(count,count1,:,:,:)=img;
%         LF(count,count1,:,:,:)=imread(['input_Cam00' num2str(N-1) '.png']);
    elseif((N-1)<100)
        img=imresize(im2double(imread(['input_Cam0' num2str(N-1) '.png'])),1);
%         img_gy=rgb2gray(img);
%         img_1= flip(img_gy,2);
%         [Gx1,Gy1] = imgradientxy(img_gy,'central');
%         [Gx2,Gy2] = imgradientxy(img_1,'central');
%         
%         img_3=(Gx1+Gy1);
%         img_4=(Gx2+Gy2);
% 
%         img_4= flip(img_4,2);
%         
%         img=im2uint8(img + 2*img_3 + 2*img_4);
        LF(count,count1,:,:,:)=img;
%         LF(count,count1,:,:,:)=imread(['input_Cam0' num2str(N-1) '.png']);
    else
        try
            img=imresize(im2double(imread(['input_Cam' num2str(N-1) '.png'])),1);
        end
%         try
%             img=imresize(im2double(imread(['input_Cam0' num2str(N-1) '.png'])),1);  %for real images
%         end
%         img_gy=rgb2gray(img);
%         img_1= flip(img_gy,2);
%         [Gx1,Gy1] = imgradientxy(img_gy,'central');
%         [Gx2,Gy2] = imgradientxy(img_1,'central');
%         
%         img_3=(Gx1+Gy1);
%         img_4=(Gx2+Gy2);
% 
%         img_4= flip(img_4,2);
%         
%         img=im2uint8(img + 2*img_3 + 2*img_4);
        LF(count,count1,:,:,:)=img;
%         LF(count,count1,:,:,:)=imread(['input_Cam' num2str(N-1) '.png']);
    end
        if(mod((N),no_sub)==0)
            count=count+1;
        end
        count1=count1+1;
        if((count1)==(no_sub+1) && N~= (no_sub^2-1))
            count1=1;
        end

end




%%
% count=2;
% count1=1;
% for N=0:8
% %     a=imresize(imread(['input_Cam00' num2str(N) '.png']),0.5);
% %     LF(1,N+1,:,:,:)=a;
%         img=imresize(im2double(imread(['input_Cam00' num2str(N) '.png'])),1);
%         img_gy=rgb2gray(img);
%         img_1= flip(img_gy,2);
%         [Gx1,Gy1] = imgradientxy(img_gy,'central');
%         [Gx2,Gy2] = imgradientxy(img_1,'central');
%         
%         img_3=(Gx1+Gy1);
%         img_4=(Gx2+Gy2);
% 
%         img_4= flip(img_4,2);
%         
%         img=im2uint8(img + 4*img_3 + 4*img_4);
% LF(1,N+1,:,:,:)=img;
% end
% 
% for N=9:80
%     if (N==9)
%         img=imresize(im2double(imread(['input_Cam00' num2str(N) '.png'])),1);
%         img_gy=rgb2gray(img);
%         img_1= flip(img_gy,2);
%         [Gx1,Gy1] = imgradientxy(img_gy,'central');
%         [Gx2,Gy2] = imgradientxy(img_1,'central');
%         
%         img_3=(Gx1+Gy1);
%         img_4=(Gx2+Gy2);
% 
%         img_4= flip(img_4,2);
%         
%         img=im2uint8(img + 4*img_3 + 4*img_4);
%         LF(count,count1,:,:,:)=img;
% %             a=imresize(imread(['input_Cam00' num2str(N) '.png']),0.5);
% %             LF(count,count1,:,:,:)=a;
%         count1=count1+1;
%     else
%         img=imresize(im2double(imread(['input_Cam0' num2str(N) '.png'])),1);
%         img_gy=rgb2gray(img);
%         img_1= flip(img_gy,2);
%         [Gx1,Gy1] = imgradientxy(img_gy,'central');
%         [Gx2,Gy2] = imgradientxy(img_1,'central');
%         
%         img_3=(Gx1+Gy1);
%         img_4=(Gx2+Gy2);
% 
%         img_4= flip(img_4,2);
%         
%         img=im2uint8(img + 4*img_3 + 4*img_4);
%         
%         
%         LF(count,count1,:,:,:)=img;
% %             a=imresize(imread(['input_Cam0' num2str(N) '.png']),0.5);
% %             LF(count,count1,:,:,:)=a;
%         if(mod((N+1),9)==0)
%             count=count+1;
%         end
%         count1=count1+1;
%         if((count1)==10 && N~=80 )
%             count1=1;
%         end
%     end
% 
% end
% save(['C:\MATLAB\Add-Ons\lytro_data\Depth-Estimation-Light-Field-master\LF\DepthEstimation(Ch4)\dataset\images for testing\synthetic LF image matlab files\' num2str(i)], ...
%       'LF' )
% 
% a=fn_ViewLightField(LF);

end