function [FinalImg1_reconstructed] = reconstruction_parallel(blurred,Shift2)

FinalImg1_reconstructed=zeros(size(blurred,1),size(blurred,2),3);
if(size(blurred,1)==size(Shift2{1},1))
    shift_temp=0;
else
   shift_temp=7;
end

for i=1:size(blurred,1)   %only for lorikeet image i=1:size(blurred,1) -1
        for j=1:size(blurred,2) 
           if (isnan(blurred(i,j)))
               FinalImg1_reconstructed(i,j,:)=0;
           else
               a=round(blurred(i,j));
               if (blurred(i,j)==0 )
                   FinalImg1_reconstructed(i,j,:)=0;
               else
    %          FinalImg1_reconstructed(i,j,:)= Shift2{a}(i+7,j+7,:); %for old code
                   FinalImg1_reconstructed(i,j,:)= Shift2{a}(i+shift_temp,j+shift_temp,:);%for new code

    %               FinalImg1_reconstructed(i,j,:)= Shift2{a}(i+6,j+6,:); %resize
               end
           end
        end
end



end