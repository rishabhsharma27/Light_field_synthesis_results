function [Shift_1,Shift_2] = grad_add(Shift3,count)

 c=cell(1,count-1);

 Shift_1 = coder.nullcopy(c);
 Shift_2 = coder.nullcopy(c);

    for i=1:count-1
        Shift2_temp_1=imresize(Shift3{i},1); %if resize required else 
%             Shift2_temp_1=(Shift3{i});
            
%         Shift2_temp_1=Shift2_temp_1(257:512,1:256,:); %test remember to comment

%         Shift2_temp_1=im2uint8(Shift3{i});
%         F = griddedInterpolant(double(Shift2_temp_1));
%         [sx,sy,sz] = size(Shift2_temp_1);
%         xq = (1:0.4995:sx)';
%         yq = (1:0.4995:sy)';
%         zq = (1:sz)';
%         F.Method = 'linear';
%         Shift2_temp_1 = im2double(F({xq,yq,zq}));
%         Shift2_temp_1=im2single((Shift2_temp_1)/255);

        
        Shift2_temp_2= imsharpen((padarray(Shift2_temp_1,[ 7 7 ], 'both')),'Radius',2,'Amount',1);
%           Shift2_temp_2= imsharpen(Shift2_temp_1,'Radius',2,'Amount',1);
 
%         Shift2_temp_2= (padarray(Shift2_temp_1,[ 7 7 ], 'both'));

        Shift2_temp_1= rgb2gray(Shift2_temp_2) ;

        Shift5_1= flip(Shift2_temp_1 ,2);

    [Gx1,Gy1] = imgradientxy(Shift2_temp_1,'central');
    [Gx2,Gy2] = imgradientxy(Shift5_1,'central');


        Shift7_1=Gx1+Gy1;
        Shift8_1=Gx2+Gy2;

        Shift5_1= flip(Shift8_1 ,2);

        Shift_1{i}= (Shift2_temp_2 + Shift7_1/4 + Shift5_1/4 );

        Shift_2{i}= 2*Shift7_1+ 2*Shift5_1;
%     (Shift2{i} + 2*Shift7{i}+ 2*Shift5{i} +xyz/2 )/1.5;

%         Shift_1{i}=imresize(Shift_1{i},2); %resize to one

    end

end
