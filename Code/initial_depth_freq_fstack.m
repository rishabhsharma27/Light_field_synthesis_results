function [range_depth]= initial_depth_freq_fstack(LF,xyz)
tic
% [LF]=Convert_Synthetic_images_to_LF_fucntion; %original 

% [LF]=Convert_Synthetic_images_to_LF_fucntion_add_noise; %added noise

if (size(LF,1)==9)
    LF=LF(3:7,3:7,:,:,:);
else
    LF=LF(3:13,3:13,:,:,:);
end



% LF=LF(:,:,1:256,257:512,:);
% LF=LF(2:8,2:8,:,:,:);

%% time 1 (total time)


% tstart1 = tic();

w = hann(512);
% w = blackman(512);
w=w.*w';


% tic
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
% toc


% depth_level=20;



Slope=(-4:0.2:4);


% tstart1 = tic();
Num_iter=5;
% xyz = im2single(LFDisp(LF));

abc=xyz;
figure(1);imshow(xyz);

% N1=N;
% N=10;


%% time 2
% tstart2 = tic();
NumView = 5; 
% test
count = 0;
[height, width, nB] = size(LF_sep{1,1});

        LF1 = im2single(zeros(height*NumView,width*NumView,nB)); 
        for i = 1:NumView
            for j = 1:NumView
                num2 = mod(count,10);
                num1 = floor(count/10);
                LF1(i:NumView:end,j:NumView:end,:) = LF_sep{i,j};
                count = count+1;
        %         imshow(img)
            end
        end
        LF1=im2single(LF1);

count=size(Slope,2)+1;


LF=im2single(LF);
% LF(:,:,:,:,4) = ones(size(LF(:,:,:,:,1)), 'single');
% because mex code doesn't support variable type change

f=@(Slope) LFFiltShiftSum_rish_mod_freq_domain_individ_subap( LF,LF_sep, Slope ); % change LFFiltShiftSum_updated/LFFiltShiftSum_updated_rish_mod


Shift2{1,size(Slope,2)} = [];

for i = 1:size(Slope,2)
    Shift2{1,i}=zeros(512,512,3);
end
% tic
for i = 1:size(Slope,2)
%    tic
   [ShiftImg] = feval(f,Slope(i));

    Shift2{i}=ShiftImg;

%    toc
end
% toc
% telapsed2 = toc(tstart2)

%% time 3


% tstart3 = tic();
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

% xyz(:,:,1)=(medfilt2(xyz(:,:,1), [5 5]));
% xyz(:,:,2)=(medfilt2(xyz(:,:,2), [5 5]));
% xyz(:,:,3)=(medfilt2(xyz(:,:,3), [5 5]));



clear a2 b2 c2 Gx1 Gy1 Gx2 Gy2 Gx3 Gy3 a b c a1 b1 c1 ;
% xyz=(xyz + xyz1 + xyz2 + xyz/2)/1.5;

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
 


[Shift,grad_Shift] = grad_add(Shift2,count); %change scale to 2
% [Shift_1,grad_Shift_1] = grad_add(Shift3,count);

% Shift=Shift2;

count2=1;
for m=1:a
    for n=1:b 
        imgfocus = (xyz((m-1)*w_size1+1:(m-1)*w_size1+w_size1,(n-1)*w_size1+1:(n-1)*w_size1+w_size1,: ));
        f=fftshift(abs(fft2(imgfocus)));
        
        imgfocus1=(xyz((m-1)*w_size1+w_start2-w_vary:(m-1)*w_size1+w_end2+w_vary,(n-1)*w_size1+w_start2-w_vary:(n-1)*w_size1+w_end2 +w_vary,:));
        f1=fftshift(abs(fft2(imgfocus1)));

        for i=1:size(Slope,2)
            img1 = (Shift{i}((m-1)*w_size1+1:(m-1)*w_size1+w_size1,(n-1)*w_size1+1:(n-1)*w_size1+w_size1,: ));
            f2{i}=fftshift(abs(fft2(img1)));
            Sharpness1(i) = immse_mex( f , f2{i});
            
            img2 = (Shift{i}((m-1)*w_size1+w_start2-w_vary:(m-1)*w_size1+w_end2+w_vary,(n-1)*w_size1+w_start2-w_vary:(n-1)*w_size1+w_end2 +w_vary,:));
            f3{i}=fftshift(abs(fft2(img2)));
            Sharpness2(i) = immse_mex( f1  , f3{i}); 
            
        end
%         figure(1),plot(Sharpness1),pause(0.3);
       [x1, y1(count2)]= min(Sharpness1,[],2); %for IMMSE function
       
       [x2, y2(count2)]= min(Sharpness2,[],2);
       count2=count2+1;
       Sharpness1 = 0;
%        Sharpness2 = 0; 
    end
end



count10=1;
for m=1:a
    for n=1:b

      FinalImg51(m,n)=  y1(count10);
      FinalImg61(m,n)= y2(count10);        
           

      count10=count10+1; 
     
     
    end
end


FinalImg2=(FinalImg51);

FinalImg2=FinalImg2(2:end-1,2:end-1);



for i=1:Num_iter
    FinalImg2=round(medfilt2(FinalImg2, [3 3]));
end

FinalImg2=FinalImg2(2:end-1,2:end-1);

% for k=1:size(Slope,2)
%     for i=1:size(FinalImg2,1)
%         for j=1:size(FinalImg2,2)
%             if ((FinalImg2(i,j))==k)
%                 depth_Our(i,j)= Slope(1,k);
%             end
%         end
%     end
% end

C=unique(FinalImg2);

for i=1:size(C,1)
    nnz1(i)=nnz(FinalImg2 == C(i));
end

%% removing error pixels
for i=1:size(C,1)
    if(nnz1(i)<30)
        L = bwlabel(FinalImg2==C(i));
        c1=unique(L);
        c1=c1(2:end);
        nnz2=0;
        for j=1:size(c1,1)
            nnz2(j)=nnz(L== j );
        end
        if (max(nnz2)>10)
            range(1,1)=(C(i));
        else
            continue;
        end
    else
        range(1,1)=(C(i));
        break;
    end
end

for i=size(C,1):-1:1
    if(nnz1(i)<30)
        L = bwlabel(FinalImg2==C(i));
        c1=unique(L);
        c1=c1(2:end);
        nnz2=0;
        for j=1:size(c1,1)
            nnz2(j)=nnz(L== j );
        end
        if (max(nnz2)>10)
            range(1,2)=(C(i));
        else
            continue;
        end
    else
        range(1,2)=(C(i));
        break;
    end
end
%%
toc

% xyz = im2single(LFDisp(LF));
% a=0;
% for i=1:41
% a=Shift2{i}+a;
% end
% xyz1=a/41;
% xyz = im2single(LFDisp(LF));
% [s1,s2]=ssim(xyz1,xyz);


range_depth(1,1)=Slope(1,range(1,1))-0.4;
range_depth(1,2)=Slope(1,range(1,2))+0.4;

if (range_depth(1,1)<-4)
    range_depth(1,1)=-4;
end

if (range_depth(1,2)>4)
    range_depth(1,2)=4;
end



toc
end