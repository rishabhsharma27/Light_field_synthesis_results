function [LF,depth_Our,Slope,FinalImg1_reconstructed,Shift2]=Our_depth_code(N)
    
%28th sep 2021 for result evaluation

%for no-noise image, generate 'xyz' after gradient additon,
%where as for noisy image generate 'xyz' before gradient additon
    
% dividing and patch depth estimation in single step. 
    

% [LF,xyz]=Convert_Synthetic_images_to_LF_fucntion_add_noise(N,0.1); %added noise

[LF]=Convert_Synthetic_images_to_LF_fucntion(N); %original

% depth_output = parsePfm( 'gt_disp_lowres.pfm');
% depth_output=im2single(round(flip(depth_output,1),1));

xyz = im2single(LFDisp(LF));
figure,imagesc(xyz),colormap(gray);

% [range_depth]= initial_depth_freq_fstack(LF,xyz);

% for i=1:size(LF,1)
%     for j=1:size(LF,2)
% %         LF_sep{i,j}=im2single(squeeze(LF(i,j,:,:,:)));
%         img=im2double(squeeze(LF(i,j,:,:,:)));
% 
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
%         img=im2double(img + 2*img_3 + 2*img_4);
%         LF(i,j,:,:,:)=img;
% 
% 
%     end
% end

[range_depth]= initial_depth_freq_fstack(LF,xyz);

step = 0.01;

Slope = range_depth(1):step:range_depth(2);
% Slope = -3.2:step:3.4;

% step=round(((max(max(depth_output)) + ( -1.* min(min(depth_output))))/4),1);
% 
% if (abs(max(max(depth_output))) > abs(min(min(depth_output))))
%     temp_s1=0:step:max(max(depth_output));
%     B = flip(temp_s1);
%     Slope=cat(2,-B,temp_s1);
%     [k1,k2] = find(Slope==0);
%     Slope(k2(1)) = [];
% else
%     temp_s1 = 0:step:abs(min(min(depth_output)));
%     B = flip(temp_s1);
%     Slope=cat(2,-B,temp_s1);
%     [k1,k2] = find(Slope==0);
%     Slope(k2(1)) = [];
% end


%Slope=(range_depth(1):step:range_depth(2));
% temp_s1 = min(min(depth_output)):step:0;
% temp_s2=0:step:max(max(depth_output));
% Slope(1,1:size(temp_s1,2))=min(min(depth_output)):step:0;
% Slope(1,size(temp_s1,2)+1:size(temp_s1,2)+size(temp_s2,2))=0:step:max(max(depth_output));
% Slope = round(Slope,1);


% Slope=0;
% 
% if (abs(max(max(depth_output))) > abs(min(min(depth_output))))
% 
%     step=(2.*abs(max(max(depth_output))))/10;
%     Slope = -abs(max(max(depth_output))):step:abs(max(max(depth_output)));
%     Slope = round(Slope,1);
% else 
%     step=(2.*abs(min(min(depth_output))))/10;
%     Slope = -abs(min(min(depth_output))):step:abs(min(min(depth_output)));
%     Slope = round(Slope,1);
% end



%Slope=(range_depth(1):step:range_depth(2));



% Slope = round(Slope,1);
% depth_level=100;
% % 
% step=(abs(range_depth(1)) + abs(range_depth(2)))/depth_level;
% Slope=(range_depth(1):step:range_depth(2));

% [LF]=Convert_Synthetic_images_to_LF_fucntion_add_noise; %added noise




% LF=LF(:,:,1:256,257:512,:);
% LF=LF(2:8,2:8,:,:,:);

%% time 1 (total time)


tstart1 = tic();

w = hann(512);
% w = blackman(512);
w=w.*w';

% net = denoisingNetwork('dncnn');
tic
for i=1:size(LF,1)
    for j=1:size(LF,2)
%         LF_sep{i,j}=im2single(squeeze(LF(i,j,:,:,:)));
        a1=im2single(squeeze(LF(i,j,:,:,1)));
        a2=im2single(squeeze(LF(i,j,:,:,2)));
        a3=im2single(squeeze(LF(i,j,:,:,3)));
        LF_sep{i,j}(:,:,1)=ifftshift(fftshift(fft2(a1),3).*w); 
        LF_sep{i,j}(:,:,2)=ifftshift(fftshift(fft2(a2),3).*w); 
        LF_sep{i,j}(:,:,3)=ifftshift(fftshift(fft2(a3),3).*w); 

    end
end
toc


% tstart1 = tic();
Num_iter=5;


abc=xyz;
xyz=imresize(xyz,1);
% xyz=xyz(257:512,1:256,:);
xyz=imfill(xyz,26, 'holes');


xyz=(padarray(xyz,[ 7 7 ], 'both'));   %for window size of 4


xyz=imsharpen(xyz,'Radius',2,'Amount',1);

a=xyz;
b = flip(xyz ,2); 
c = flip(xyz ,1);
a1=rgb2gray(a);b1=rgb2gray(b);c1=rgb2gray(c);
[Gx1,Gy1] = imgradientxy(a1,'central');
[Gx2,Gy2] = imgradientxy(b1,'central');
[Gx3,Gy3] = imgradientxy(c1,'central');
a2=Gx1+Gy1;
b2=Gx2+Gy2;
c2=Gx3+Gy3;

b3 = flip(b2 ,2);
c3 = flip(c2 ,1);

xyz1=(a2*2);
xyz2=(b3*2);

% xyz_add = xyz;
xyz=( xyz + xyz1/4 + xyz2/4 );

clear a2 b2 c2 Gx1 Gy1 Gx2 Gy2 Gx3 Gy3 a b c a1 b1 c1 ;
% xyz=(xyz + xyz1 + xyz2 + xyz/2)/1.5;


figure(1);imshow(xyz);



%% time 2
tstart2 = tic();

NumView = 9;  % 9 

% NumView = 7;  % 9 
% LF=LF(2:8,2:8,:,:,:);

count=size(Slope,2)+1;


LF=im2single(LF);


f=@(Slope) LFFiltShiftSum_rish_mod_freq_domain_individ_subap( LF,LF_sep, Slope ); % change LFFiltShiftSum_updated/LFFiltShiftSum_updated_rish_mod

%f=@(Slope) LFFiltShiftSum_updated( LF, Slope );

Shift2{1,size(Slope,2)} = [];

for i = 1:size(Slope,2)
    Shift2{1,i}=zeros(512,512,3);
end
% tic
for i = 1:size(Slope,2)
%    tic
   [ShiftImg] = feval(f,Slope(i));
   % ShiftImg = squeeze(ShiftImg(:,:,1:end-1));
   ShiftImg = squeeze(ShiftImg(:,:,:));
    Shift2{i}=squeeze(ShiftImg);

%    toc
end
% toc
t_fstack_gen = toc(tstart2)

%% time 3

% tstart3 = tic();

h = size(xyz,1);
w = size(xyz,2);

%octagonal window size 
w_size1= 4;
w_width1= 1;
w_vary= 1;  %when the w_size1 changes to maintain the size w_size1 the adjust is made


a=floor(h/w_size1);     %height of image divide by window size
b=floor(w/w_size1);     %width of image divide by window size
c= a*b;       %number of windows i.e a*b



%% Octagonal window1 size 
%2 octagonal windows are used the size of one defines the size of the other window ) 

[g,g1, w_size2 ,w_width2]=octagon_gen(w_size1,w_width1);

w_start2=w_size1-((w_size2)/2)+1;

w_end2=w_start2+w_size2-1;
g2=g1;
g1=padarray(g1,[w_vary w_vary], 'both');

H = fspecial('disk',1);
 


[Shift,grad_Shift] = grad_add(Shift2,count); 

Shift3=Shift;




%change scale to 2
% [Shift_1,grad_Shift_1] = grad_add(Shift3,count);

% Shift=Shift2;

%clear grad_Shift Shift2;


%% time 3

tstart3 = tic();
count2=1;
for m=1:a
    for n=1:b
        imgfocus = (xyz((m-1)*w_size1+1:(m-1)*w_size1+w_size1,(n-1)*w_size1+1:(n-1)*w_size1+w_size1,: ));
        imgfocus1=(xyz((m-1)*w_size1+w_start2-w_vary:(m-1)*w_size1+w_end2+w_vary,(n-1)*w_size1+w_start2-w_vary:(n-1)*w_size1+w_end2 +w_vary,:));
        f=fftshift(abs(fft2(imgfocus)));
        f1=fftshift(abs(fft2(imgfocus1)));

        for i=1:size(Slope,2)
            img1{i} = (Shift{i}((m-1)*w_size1+1:(m-1)*w_size1+w_size1,(n-1)*w_size1+1:(n-1)*w_size1+w_size1,: ));
            img2{i} = (Shift{i}((m-1)*w_size1+w_start2-w_vary:(m-1)*w_size1+w_end2+w_vary,(n-1)*w_size1+w_start2-w_vary:(n-1)*w_size1+w_end2 +w_vary,:));
            f2{i}=fftshift(abs(fft2(img1{i})));
            f3{i}=fftshift(abs(fft2(img2{i})));
            Sharpness1(i) = immse_mex( f , f2{i});
            Sharpness2(i) = immse_mex( f1  , f3{i}); 
        end
%         figure(1),plot(Sharpness1),pause(0.3);
       [x1, y1(count2)]= min(Sharpness1,[],2); %for IMMSE function
       [x2, y2(count2)]= min(Sharpness2,[],2);
       count2=count2+1;
%        Sharpness1 = 0;
%        Sharpness2 = 0; 
    end
end
t_depth_gen = toc(tstart3)


count10=1;
for m=1:a
    for n=1:b

      FinalImg51((m-1)*w_size1+1:(m-1)*w_size1+w_size1,(n-1)*w_size1+1:(n-1)*w_size1+w_size1)=  y1(count10).*g;
      FinalImg61((m-1)*w_size1+w_start2-w_vary:(m-1)*w_size1+w_end2+w_vary,(n-1)*w_size1+w_start2-w_vary:(n-1)*w_size1+w_end2 +w_vary)= y2(count10).*g1;        
           

      count10=count10+1; 
     
     
    end
end
sizeFinalImg5=size(FinalImg61,1);
sizeFinalImg6=size(FinalImg61,2);
FinalImg51(sizeFinalImg5,sizeFinalImg6)=0;

FinalImg2=(FinalImg61 + FinalImg51);

figure;imagesc(FinalImg2),colormap(gray)

FinalImg3=FinalImg2;

FinalImg3=round(medfilt2(FinalImg3, [  9 9 ]));


FinalImg3=FinalImg3(8:size(FinalImg3,1)-7,8:size(FinalImg3,2)-7);
FinalImg2=FinalImg2(8:size(FinalImg2,1)-7,8:size(FinalImg2,2)-7);


for i=1:Num_iter
    FinalImg2=round(medfilt2(FinalImg2, [7 7]));
end

for i=1:(size(FinalImg3,1))
    for j=1:(size(FinalImg3,2))
        if (FinalImg2(i,j)==0)
            FinalImg2(i,j)=FinalImg3(i,j);
        end
    end
end

%FinalImg2=FinalImg2(corner1(1,1):corner1(1,2),corner1(1,3):corner1(1,4));

% telapsed3 = toc(tstart3)

for k=1:size(Slope,2)
    for i=1:size(FinalImg2,1)
        for j=1:size(FinalImg2,2)
            if ((FinalImg2(i,j))==k)
                depth_Our(i,j)= Slope(1,k);
            end
        end
    end
end

telapsed1 = toc(tstart1)

FinalImg1_reconstructed=im2single(reconstruction_parallel(FinalImg2,Shift2));



close all;

end


