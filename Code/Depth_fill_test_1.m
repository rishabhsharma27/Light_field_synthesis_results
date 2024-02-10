function [a3] = Depth_fill_test_1(img_h_depth)

error_count=0;

img_h_a2_d= img_h_depth;


    for k=1
        for j=1:size(img_h_a2_d,2)-2

            s2=round(img_h_a2_d((j-1)*9+1:(j-1)*9+9,:));
            
            mask = isnan(s2);
            
            s2(isnan(s2))=0;

%             s2=imfill(s2);
            
            s2=round(inpaintCoherent(s2,mask,'Radius',1));

            img_h_a2_d((j-1)*9+1:(j-1)*9+9,:)=(s2);
        end
    end

%     img_h_a2_d=(medfilt2(img_h_a2_d, [  3 3 ]));

    
count=1;
for j=1:9
    for i=j:9:size(img_h_a2_d,1)
        a3{j}(count,:)= img_h_a2_d(i,:);
        a6{j}(count,:)= img_h_depth(i,:,:);
%         a4{j}(count,:,:)= img_h_view2(i,:,:);  
%         a5{j}(count,:,:)= img_h_view3(i,:,:);
        count=count+1;
    end
    count=1;
end

% for k=1:9
%     a3{k}=round(imgaussfilt3(a3{k},0.5),3);
% end

for k=1:4
     for i=1:size(a3{1},1)
        for j=1:size(a3{1},2)
            if(isnan(a6{k}(i,j,1)) )
                a6{k}(i,j,:)=a3{1}(i,j,:); 
            else
                a6{k}(i,j,:)=a6{k}(i,j,:);
            end
        end
    end
end

for k=6:9
     for i=1:size(a3{1},1)
        for j=1:size(a3{1},2)
            if(isnan(a6{k}(i,j,1)) )
                a6{k}(i,j,:)=a3{9}(i,j,:); 
            else
                a6{k}(i,j,:)=a6{k}(i,j,:);
            end
        end
    end
end



for i=1:9
    a6{i}=(medfilt2(a6{i}, [ 3 3 ]));
    a6{i}=round(a6{i});
    mask=isnan(a6{i}); 
    a6{i}=round(regionfill(a6{i},mask));
end

a3=a6;
% for j=1:9
%     figure(j),imagesc(a3{j}),colormap(gray),pause(1)
% end

end

% B =zeros(size(s,1),size(s,2),3);
% for l=1:9
%    B(l,:,:) = repmat(s(l,:,:),1); 
% end


