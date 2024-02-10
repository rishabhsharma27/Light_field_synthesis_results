function [a1] = RGB_fill_horizontal_test_2022(a1,a2,Slope,Shift2,a2_opp,var_fill)  
load('C:\MATLAB\Add-Ons\lytro_data\Depth-Estimation-Light-Field-master\LF\DepthEstimation(Ch4)\dataset\blender_LF_render_file\blur_reduction_var.mat', 'PSF')
    for I=1:9
        for i=1:512
            for j=1:512
                if(isnan(a1{I}(i,j,:)))
                    a1{I}(i,j,:)=0;
                end
            end
        end
    end    
    Slope1 = round(Slope.*10);

    %y = p1*x + p2;
      p1 = 4;
      p2 = -0;

    x = 0:abs(Slope(1,1)-Slope(1,2)):4.1;
%     x = 0:abs(Slope(1,1)-Slope(1,2)):Slope(1,end);
    
    fx = double(round(p1*x + p2));
    if (var_fill ==0) 
        fx2(1,1:(size(Slope1,2)+1)/2) = 0;
        fx2 = double(fx2);
        for I=1:4
        fx1=double(round(fx/I));
        
        
        temp = round(a2{I}.*10);
        temp_opp = round(a2_opp{I}.*10);
        for i=1:(size(Slope1,2)+1)/2-1
%             abc=0;
            mask1 = (temp==Slope1(i));
            if (mask1==0)
                continue;
            else
                mask1 = mask1.*(a1{I}==0);
                mask_opp = temp_opp.*mask1;
                C = unique(mask_opp);
                abc = Shift2{i}(1:512-fx2(1,size(fx2,2)-i),1:512-fx1(1,size(fx1,2)-i),:);
                abc =(padarray(abc,[fx2(1,size(fx2,2)-i) fx1(1,size(fx1,2)-i)], 'pre'));
%                 X = find(C==0);
%                 C(X)=[];
                if isempty(C)
                    continue;
                else
                    
                    for j=1:size(C,1)
                        if (C(j)==Slope1(i))
%                             continue;
                            mask1 = mask1.*(a1{I}==0);
                            temp_rgb = abc.*mask1;
                            a1{I} = a1{I} + temp_rgb;
                        else
                            mask2 = (temp_opp==C(j)).*mask1;
                            if (mask2==0)
                                continue;
                            else
                                X = find(C(j)==Slope1);
                                X1 = (find(abs(C(j))==Slope1)-1)/2;
                                abc1 = Shift2{X}(1:512-fx2(1,size(fx2,2)-i),1:512-fx1(1,size(fx1,2)-i),:);
                                abc1 =(padarray(abc1,[fx2(1,size(fx2,2)-i) fx1(1,size(fx1,2)-i)], 'pre'));
                                C1(:,:,1) = conv2(abc1(:,:,1),PSF{X1*2});
                                C1(:,:,2) = conv2(abc1(:,:,2),PSF{X1*2});
                                C1(:,:,3) = conv2(abc1(:,:,3),PSF{X1*2});
                                size_C1 = ((size(C1,1) - size(abc1,1))+1)/2;
                                abc1 = C1(size_C1:size(C1,1)-size_C1,size_C1:size(C1,2)-size_C1,:);
                                clear C1;
                                temp_rgb = abc.*mask2 - (abc1./1.5).*mask2;
                                temp_rgb = imhistmatchn(temp_rgb,abc.*mask2);
                                temp_rgb = temp_rgb.*mask2;
                                a1{I} = a1{I} + temp_rgb;
                            end
                       end
                    end
                end
            end
        end
        %figure,imagesc(a1{I}),colormap(gray),axis off

        temp = round(a2{I}.*10);
        temp_opp = round(a2_opp{I}.*10);
        for i=(size(Slope1,2)+1)/2+1:size(Slope1,2)
%         for i=size(Slope1,2) %test change back to above
%             abc=0;
            
            mask1 = (temp==Slope1(i));
            if (mask1==0)
                continue;
            else
                mask1 = mask1.*(a1{I}==0);
                mask_opp = temp_opp.*mask1;
                C = unique(mask_opp);
                abc = Shift2{i}(1+fx2(1,i-size(fx2,2)):512,1+fx1(1,i-size(fx1,2)):512,:);
                abc =(padarray(abc,[fx2(1,i-size(fx2,2)) fx1(1,i-size(fx1,2))], 'post'));
%                 X = find(C==0);
%                 C(X)=[];
                if isempty(C)
                    continue;
                else
                        for j=1:size(C,1)
                            if (C(j)==Slope1(i))
%                                 continue;
                                mask1 = mask1.*(a1{I}==0);
                                temp_rgb = abc.*mask1;
                                a1{I} = a1{I} + temp_rgb;
                            else
                                mask2 = (temp_opp==C(j)).*mask1;
                                if (mask2==0)
                                    continue;
                                else
                                    X = find(C(j)==Slope1);
                                    X1 = (find(abs(C(j))==Slope1)-1)/2;
                                    abc1 = Shift2{X}(1+fx2(1,i-size(fx2,2)):512,1+fx1(1,i-size(fx1,2)):512,:);
                                    abc1 =(padarray(abc1,[fx2(1,i-size(fx2,2)) fx1(1,i-size(fx1,2))], 'post'));
                                    C1(:,:,1) = conv2(abc1(:,:,1),PSF{X1*2});
                                    C1(:,:,2) = conv2(abc1(:,:,2),PSF{X1*2});
                                    C1(:,:,3) = conv2(abc1(:,:,3),PSF{X1*2});
                                    size_C1 = ((size(C1,1) - size(abc1,1))+1)/2;
                                    abc1 = C1(size_C1:size(C1,1)-size_C1,size_C1:size(C1,1)-size_C1,:);
                                    clear C1;
                                    temp_rgb = abc.*mask2 - (abc1./1.5).*mask2;
                                    temp_rgb = imhistmatchn(temp_rgb,abc.*mask2);
                                    temp_rgb = temp_rgb.*mask2;
                                    a1{I} = a1{I} + temp_rgb;
%                                 figure,imagesc(a1{1}),axis off,pause(1);
                                end
                           end
                        end
                end
            end
        end
        temp = round(a2{I}.*10);
        mask1 = (temp==Slope1((size(Slope1,2)+1)/2));
        mask1 = mask1.*(a1{I}==0);
        abc = Shift2{((size(Slope1,2)+1)/2)};
        temp_rgb = abc.*mask1;
        a1{I} = a1{I} + temp_rgb;
    end
    
p1 = 4;
p2 = -0;

%x = 0:abs(Slope(1,1)-Slope(1,2)):4.1;
x = 0:abs(Slope(1,1)-Slope(1,2)):Slope(1,end);
fx = double(round(p1*x + p2));

    for I=6:9
        fx1=double(round(fx/(10-I)));

        
        temp = round(a2{I}.*10);
        temp_opp = round(a2_opp{I}.*10);
        for i=(size(Slope1,2)+1)/2+1:size(Slope1,2)
%             abc=0;
            mask1 = (temp==Slope1(i));
            if (mask1==0)
                continue;
            else
                mask1 = mask1.*(a1{I}==0);
                mask_opp = temp_opp.*mask1;
                C = unique(mask_opp);
                abc = Shift2{i}(1+fx2(1,i-size(fx2,2)):512,1:512-fx1(1,i-size(fx1,2)),:);
                abc =(padarray(abc,[0 fx1(1,i-size(fx1,2))], 'pre'));
                abc =(padarray(abc,[fx2(1,i-size(fx2,2)) 0], 'post'));
%                 X = find(C==0);
%                 C(X)=[];
                if isempty(C)
                    continue;
                else
                    for j=1:size(C,1)
                        if (C(j)==Slope1(i))
%                             continue;
                            mask1 = mask1.*(a1{I}==0);
                            temp_rgb = abc.*mask1;
                            a1{I} = a1{I} + temp_rgb;
                        else
                            mask2 = (temp_opp==C(j)).*mask1;
                            if (mask2==0)
                                continue;
                            else
                                X = find(C(j)==Slope1);
                                X1 = (find(abs(C(j))==Slope1)-1)/2;
                                abc1 = Shift2{X}(1+fx2(1,i-size(fx2,2)):512,1:512-fx1(1,i-size(fx1,2)),:);
                                abc1 =(padarray(abc1,[0 fx1(1,i-size(fx1,2))], 'pre'));
                                abc1 =(padarray(abc1,[fx2(1,i-size(fx2,2)) 0], 'post'));
                                C1(:,:,1) = conv2(abc1(:,:,1),PSF{X1*2});
                                C1(:,:,2) = conv2(abc1(:,:,2),PSF{X1*2});
                                C1(:,:,3) = conv2(abc1(:,:,3),PSF{X1*2});
                                size_C1 = ((size(C1,1) - size(abc1,1))+1)/2;
                                abc1 = C1(size_C1:size(C1,1)-size_C1,size_C1:size(C1,1)-size_C1,:);
                                clear C1;                   
                                temp_rgb = abc.*mask2 - (abc1./1.5).*mask2;
                                temp_rgb = imhistmatchn(temp_rgb,abc.*mask2);
                                temp_rgb = temp_rgb.*mask2;
                                a1{I} = a1{I} + temp_rgb;
                            end
                        end
                    end
                end
            end
        end
        %figure,imagesc(a1{I}),colormap(gray),axis off

        temp = round(a2{I}.*10);
        temp_opp = round(a2_opp{I}.*10);
        for i=1:(size(Slope1,2)+1)/2-1
%             abc=0;
            mask1 = (temp==Slope1(i));
            if (mask1==0)
                continue;
            else
                mask1 = mask1.*(a1{I}==0);
                mask_opp = temp_opp.*mask1;
                C = unique(mask_opp);
                abc = Shift2{i}(1:512-fx2(1,size(fx2,2)-i),1+fx1(1,size(fx1,2)-i):512,:);
                abc =(padarray(abc,[0 fx1(1,size(fx1,2)-i)], 'post'));
                abc =(padarray(abc,[fx2(1,size(fx2,2)-i) 0], 'pre'));
%                 X = find(C==0);
%                 C(X)=[];
                if isempty(C)
                    continue;
                else
                    for j=1:size(C,1)
                        if (C(j)==Slope1(i))
                            %continue;
                            mask1 = mask1.*(a1{I}==0);
                            temp_rgb = abc.*mask1;
                            a1{I} = a1{I} + temp_rgb;
                        else
                            mask2 = (temp_opp==C(j)).*mask1;
                            if (mask2==0)
                                continue;
                            else
                                X = find(C(j)==Slope1);
                                X1 = (find(abs(C(j))==Slope1)-1)/2;
                                abc1 = Shift2{X}(1:512-fx2(1,size(fx2,2)-i),1+fx1(1,size(fx1,2)-i):512,:);
                                abc1 =(padarray(abc1,[0 fx1(1,size(fx1,2)-i)], 'post'));
                                abc1 =(padarray(abc1,[fx2(1,size(fx2,2)-i) 0], 'pre'));
                                C1(:,:,1) = conv2(abc1(:,:,1),PSF{X1*2});
                                C1(:,:,2) = conv2(abc1(:,:,2),PSF{X1*2});
                                C1(:,:,3) = conv2(abc1(:,:,3),PSF{X1*2});
                                size_C1 = ((size(C1,1) - size(abc1,1))+1)/2;
                                abc1 = C1(size_C1:size(C1,1)-size_C1,size_C1:size(C1,1)-size_C1,:);
                                clear C1;
                                temp_rgb = abc.*mask2 - (abc1./1.5).*mask2;
                                temp_rgb = imhistmatchn(temp_rgb,abc.*mask2);
                                temp_rgb = temp_rgb.*mask2;
                                a1{I} = a1{I} + temp_rgb;
                            end
                        end
                    end
                end
            end
        end
        temp = round(a2{I}.*10);
        mask1 = (temp==Slope1((size(Slope1,2)+1)/2));
        mask1 = mask1.*(a1{I}==0);
        abc = Shift2{((size(Slope1,2)+1)/2)};
        temp_rgb = abc.*mask1;
        a1{I} = a1{I} + temp_rgb;
    end 
%% for 1:4 a1_1{} - a1_4{} sub-aperture images as the focal stack images need to be shifted in vertical direction      
    elseif (var_fill >0)

        fx2=double(round(fx/var_fill));
        for I=1:4
        fx1=double(round(fx/I));
        
        
        temp = round(a2{I}.*10);
        temp_opp = round(a2_opp{I}.*10);        
        for i=1:(size(Slope1,2)+1)/2-1
%             abc=0;
            mask1 = (temp==Slope1(i));
            if (mask1==0)
                continue;
            else
                mask1 = mask1.*(a1{I}==0);
                mask_opp = temp_opp.*mask1;
                C = unique(mask_opp);
                abc = Shift2{i}(1:512-fx2(1,size(fx2,2)-i),1:512-fx1(1,size(fx1,2)-i),:);
                abc =(padarray(abc,[fx2(1,size(fx2,2)-i) fx1(1,size(fx1,2)-i)], 'pre'));
                if isempty(C)
                    continue;
                else
                    for j=1:size(C,1)
                        if (C(j)==Slope1(i))
                            %continue;
                            mask1 = mask1.*(a1{I}==0);
                            temp_rgb = abc.*mask1;
                            a1{I} = a1{I} + temp_rgb.*2;
                        else
                            mask2 = (temp_opp==C(j)).*mask1;
                            if (mask2==0)
                                continue;
                            else
                                X = find(C(j)==Slope1);
                                X1 = find(abs(C(j))==Slope1)-1;
                                abc1 = Shift2{X}(1:512-fx2(1,size(fx2,2)-i),1:512-fx1(1,size(fx1,2)-i),:);
                                abc1 =(padarray(abc1,[fx2(1,size(fx2,2)-i) fx1(1,size(fx1,2)-i)], 'pre'));
                                C1(:,:,1) = conv2(abc1(:,:,1),PSF{X1*2});
                                C1(:,:,2) = conv2(abc1(:,:,2),PSF{X1*2});
                                C1(:,:,3) = conv2(abc1(:,:,3),PSF{X1*2});
                                size_C1 = ((size(C1,1) - size(abc1,1))+1)/2;
                                abc1 = C1(size_C1:size(C1,1)-size_C1,size_C1:size(C1,2)-size_C1,:);
                                clear C1;
                                temp_rgb = abc.*mask2 - (abc1./1.5).*mask2;
                                temp_rgb = imhistmatchn(temp_rgb,abc.*mask2);
                                temp_rgb = temp_rgb.*mask2;
                                a1{I} = a1{I} + temp_rgb;
                            end
                        end
                    end
                end
            end
        end
        %figure,imagesc(a1{I}),colormap(gray),axis off

        temp = round(a2{I}.*10);
        temp_opp = round(a2_opp{I}.*10);
        for i=(size(Slope1,2)+1)/2+1:size(Slope1,2)
%             abc=0;
            mask1 = (temp==Slope1(i));
            if (mask1==0)
                continue;
            else
                mask1 = mask1.*(a1{I}==0);
                mask_opp = temp_opp.*mask1;
                C = unique(mask_opp);
                abc = Shift2{i}(1+fx2(1,i-size(fx2,2)):512,1+fx1(1,i-size(fx1,2)):512,:);
                abc =(padarray(abc,[fx2(1,i-size(fx2,2)) fx1(1,i-size(fx1,2))], 'post'));
                if isempty(C)
                    continue;
                else
                        for j=1:size(C,1)
                            if (C(j)==Slope1(i))
                            %continue;
                            mask1 = mask1.*(a1{I}==0);
                            temp_rgb = abc.*mask1;
                            a1{I} = a1{I} + temp_rgb;
                            else
                                mask2 = (temp_opp==C(j)).*mask1;
                                if (mask2==0)
                                    continue;
                                else
                                    X = find(C(j)==Slope1);
                                    X1 = find(abs(C(j))==Slope1)-1;
                                    abc1 = Shift2{X}(1+fx2(1,i-size(fx2,2)):512,1+fx1(1,i-size(fx1,2)):512,:);
                                    abc1 =(padarray(abc1,[fx2(1,i-size(fx2,2)) fx1(1,i-size(fx1,2))], 'post'));
                                    C1(:,:,1) = conv2(abc1(:,:,1),PSF{X1*2});
                                    C1(:,:,2) = conv2(abc1(:,:,2),PSF{X1*2});
                                    C1(:,:,3) = conv2(abc1(:,:,3),PSF{X1*2});
                                    size_C1 = ((size(C1,1) - size(abc1,1))+1)/2;
                                    abc1 = C1(size_C1:size(C1,1)-size_C1,size_C1:size(C1,1)-size_C1,:);
                                    clear C1;
                                    temp_rgb = abc.*mask2 - (abc1./1.5).*mask2;
                                    temp_rgb = imhistmatchn(temp_rgb,abc.*mask2);
                                    temp_rgb = temp_rgb.*mask2;
                                    a1{I} = a1{I} + temp_rgb;
                                end
%                                 figure,imagesc(a1{1}),axis off,pause(1);
                            end
                        end
                end

            end
        end
        temp = round(a2{I}.*10);
        mask1 = (temp==Slope1((size(Slope1,2)+1)/2));
        mask1 = mask1.*(a1{I}==0);
        abc = Shift2{((size(Slope1,2)+1)/2)};
        temp_rgb = abc.*mask1;
        a1{I} = a1{I} + temp_rgb;
    end
    
p1 = 4;
p2 = -0;

% x = 0:abs(Slope(1,1)-Slope(1,2)):4.1;
x = 0:abs(Slope(1,1)-Slope(1,2)):Slope(1,end);
fx = double(round(p1*x + p2));

    for I=6:9
        fx1=double(round(fx/(10-I)));

        
        temp = round(a2{I}.*10);
        temp_opp = round(a2_opp{I}.*10);
        for i=(size(Slope1,2)+1)/2+1:size(Slope1,2)
%             abc=0;
            mask1 = (temp==Slope1(i));
            if (mask1==0)
                continue;
            else
                mask1 = mask1.*(a1{I}==0);
                mask_opp = temp_opp.*mask1;
                C = unique(mask_opp);
                abc = Shift2{i}(1+fx2(1,i-size(fx2,2)):512,1:512-fx1(1,i-size(fx1,2)),:);
                abc =(padarray(abc,[0 fx1(1,i-size(fx1,2))], 'pre'));
                abc =(padarray(abc,[fx2(1,i-size(fx2,2)) 0], 'post'));
                if isempty(C)
                    continue;
                else
                    for j=1:size(C,1)
                        if (C(j)==Slope1(i))
                            %continue;
                            mask1 = mask1.*(a1{I}==0);
                            temp_rgb = abc.*mask1;
                            a1{I} = a1{I} + temp_rgb;
                        else
                            mask2 = (temp_opp==C(j)).*mask1;
                            if (mask2==0)
                                continue;
                            else
                                X = find(C(j)==Slope1);
                                X1 = find(abs(C(j))==Slope1)-1;
                                abc1 = Shift2{X}(1+fx2(1,i-size(fx2,2)):512,1:512-fx1(1,i-size(fx1,2)),:);
                                abc1 =(padarray(abc1,[0 fx1(1,i-size(fx1,2))], 'pre'));
                                abc1 =(padarray(abc1,[fx2(1,i-size(fx2,2)) 0], 'post'));
                                C1(:,:,1) = conv2(abc1(:,:,1),PSF{X1*2});
                                C1(:,:,2) = conv2(abc1(:,:,2),PSF{X1*2});
                                C1(:,:,3) = conv2(abc1(:,:,3),PSF{X1*2});
                                size_C1 = ((size(C1,1) - size(abc1,1))+1)/2;
                                abc1 = C1(size_C1:size(C1,1)-size_C1,size_C1:size(C1,1)-size_C1,:);
                                clear C1;                   
                                temp_rgb = abc.*mask2 - (abc1./1.5).*mask2;
                                temp_rgb = imhistmatchn(temp_rgb,abc.*mask2);
                                temp_rgb = temp_rgb.*mask2;
                                a1{I} = a1{I} + temp_rgb;
                            end
                        end
                    end
                end

            end
        end
        %figure,imagesc(a1{I}),colormap(gray),axis off

        temp = round(a2{I}.*10);
        temp_opp = round(a2_opp{I}.*10);
        for i=1:(size(Slope1,2)+1)/2-1
%             abc=0;
            mask1 = (temp==Slope1(i));
            if (mask1==0)
                continue;
            else
                mask1 = mask1.*(a1{I}==0);
                mask_opp = temp_opp.*mask1;
                C = unique(mask_opp);
                abc = Shift2{i}(1:512-fx2(1,size(fx2,2)-i),1+fx1(1,size(fx1,2)-i):512,:);
                abc =(padarray(abc,[0 fx1(1,size(fx1,2)-i)], 'post'));
                abc =(padarray(abc,[fx2(1,size(fx2,2)-i) 0], 'pre'));
                if isempty(C)
                    continue;
                else
                    for j=1:size(C,1)
                        if (C(j)==Slope1(i))
                            %continue;
                            mask1 = mask1.*(a1{I}==0);
                            temp_rgb = abc.*mask1;
                            a1{I} = a1{I} + temp_rgb;
                        else
                            mask2 = (temp_opp==C(j)).*mask1;
                            if (mask2==0)
                                continue;
                            else
                                X = find(C(j)==Slope1);
                                X1 = find(abs(C(j))==Slope1)-1;
                                abc1 = Shift2{X}(1:512-fx2(1,size(fx2,2)-i),1+fx1(1,size(fx1,2)-i):512,:);
                                abc1 =(padarray(abc1,[0 fx1(1,size(fx1,2)-i)], 'post'));
                                abc1 =(padarray(abc1,[fx2(1,size(fx2,2)-i) 0], 'pre'));
                                C1(:,:,1) = conv2(abc1(:,:,1),PSF{X1*2});
                                C1(:,:,2) = conv2(abc1(:,:,2),PSF{X1*2});
                                C1(:,:,3) = conv2(abc1(:,:,3),PSF{X1*2});
                                size_C1 = ((size(C1,1) - size(abc1,1))+1)/2;
                                abc1 = C1(size_C1:size(C1,1)-size_C1,size_C1:size(C1,1)-size_C1,:);
                                clear C1;
                                temp_rgb = abc.*mask2 - (abc1./1.5).*mask2;
                                temp_rgb = imhistmatchn(temp_rgb,abc.*mask2);
                                temp_rgb = temp_rgb.*mask2;
                                a1{I} = a1{I} + temp_rgb;
                            end
                        end
                    end
                end
            end
        end
        temp = round(a2{I}.*10);
        mask1 = (temp==Slope1((size(Slope1,2)+1)/2));
        mask1 = mask1.*(a1{I}==0);
        abc = Shift2{((size(Slope1,2)+1)/2)};
        temp_rgb = abc.*mask1;
        a1{I} = a1{I} + temp_rgb;
    end

    
%% for 6:9 a1_6{} - a1_9{} sub-aperture images as the focal stack images need to be shifted in vertical direction
    elseif(var_fill <0)
        var_fill = -var_fill;
        fx2=double(round(fx/var_fill));
        for I=1:4
        fx1=double(round(fx/I));
        
        
        temp = round(a2{I}.*10);
        temp_opp = round(a2_opp{I}.*10);
        for i=1:(size(Slope1,2)+1)/2-1
%             abc=0;
            mask1 = (temp==Slope1(i));
            if (mask1==0)
                continue;
            else
                mask1 = mask1.*(a1{I}==0);
                mask_opp = temp_opp.*mask1;
                C = unique(mask_opp);
                abc = Shift2{i}(1+fx2(1,size(fx2,2)-i):512,1:512-fx1(1,size(fx1,2)-i),:);
                abc =(padarray(abc,[0 fx1(1,size(fx1,2)-i)], 'pre'));
                abc =(padarray(abc,[fx2(1,size(fx2,2)-i) 0], 'post'));
                if isempty(C)
                    continue;
                else
                    
                    for j=1:size(C,1)
                        if (C(j)==Slope1(i))
                            %continue;
                            mask1 = mask1.*(a1{I}==0);
                            temp_rgb = abc.*mask1;
                            a1{I} = a1{I} + temp_rgb;
                        else
                            mask2 = (temp_opp==C(j)).*mask1;
                            if (mask2==0)
                                continue;
                            else
                                X = find(C(j)==Slope1);
                                X1 = find(abs(C(j))==Slope1)-1;
                                abc1 = Shift2{X}(1+fx2(1,size(fx2,2)-i):512,1:512-fx1(1,size(fx1,2)-i),:);
                                abc1 =(padarray(abc1,[0 fx1(1,size(fx1,2)-i)], 'pre'));
                                abc1 =(padarray(abc1,[fx2(1,size(fx2,2)-i) 0], 'post'));
                                C1(:,:,1) = conv2(abc1(:,:,1),PSF{X1*2});
                                C1(:,:,2) = conv2(abc1(:,:,2),PSF{X1*2});
                                C1(:,:,3) = conv2(abc1(:,:,3),PSF{X1*2});
                                size_C1 = ((size(C1,1) - size(abc1,1))+1)/2;
                                abc1 = C1(size_C1:size(C1,1)-size_C1,size_C1:size(C1,2)-size_C1,:);
                                clear C1;
                                temp_rgb = abc.*mask2 - (abc1./1.5).*mask2;
                                temp_rgb = imhistmatchn(temp_rgb,abc.*mask2);
                                temp_rgb = temp_rgb.*mask2;
                                a1{I} = a1{I} + temp_rgb;
                            end
                        end
                    end
                end
            end
        end
        %figure,imagesc(a1{I}),colormap(gray),axis off

        temp = round(a2{I}.*10);
        temp_opp = round(a2_opp{I}.*10);
        for i=(size(Slope1,2)+1)/2+1:size(Slope1,2)
%             abc=0;
            mask1 = (temp==Slope1(i));
            if (mask1==0)
                continue;
            else
                mask1 = mask1.*(a1{I}==0);
                mask_opp = temp_opp.*mask1;
                C = unique(mask_opp);
                abc = Shift2{i}(1:512-fx2(1,i-size(fx2,2)),1+fx1(1,i-size(fx1,2)):512,:);
                abc =(padarray(abc,[0 fx1(1,i-size(fx1,2))], 'post'));
                abc =(padarray(abc,[fx2(1,i-size(fx2,2)) 0], 'pre'));
                if isempty(C)
                    continue;
                else
                        for j=1:size(C,1)
                            if (C(j)==Slope1(i))
                            %continue;
                            mask1 = mask1.*(a1{I}==0);
                            temp_rgb = abc.*mask1;
                            a1{I} = a1{I} + temp_rgb;
                            else
                                mask2 = (temp_opp==C(j)).*mask1;
                                if (mask2==0)
                                    continue;
                                else
                                    X = find(C(j)==Slope1);
                                    X1 = find(abs(C(j))==Slope1)-1;
                                    abc1 = Shift2{X}(1:512-fx2(1,i-size(fx2,2)),1+fx1(1,i-size(fx1,2)):512,:);
                                    abc1 =(padarray(abc1,[0 fx1(1,i-size(fx1,2))], 'post'));
                                    abc1 =(padarray(abc1,[fx2(1,i-size(fx2,2)) 0], 'pre'));
                                    C1(:,:,1) = conv2(abc1(:,:,1),PSF{X1*2});
                                    C1(:,:,2) = conv2(abc1(:,:,2),PSF{X1*2});
                                    C1(:,:,3) = conv2(abc1(:,:,3),PSF{X1*2});
                                    size_C1 = ((size(C1,1) - size(abc1,1))+1)/2;
                                    abc1 = C1(size_C1:size(C1,1)-size_C1,size_C1:size(C1,1)-size_C1,:);
                                    clear C1;
                                    temp_rgb = abc.*mask2 - (abc1./1.5).*mask2;
                                    temp_rgb = imhistmatchn(temp_rgb,abc.*mask2);
                                    temp_rgb = temp_rgb.*mask2;
                                    a1{I} = a1{I} + temp_rgb;
                                end
%                                 figure,imagesc(a1{1}),axis off,pause(1);
                            end
                        end
                end
            end
        end
        temp = round(a2{I}.*10);
        mask1 = (temp==Slope1((size(Slope1,2)+1)/2));
        mask1 = mask1.*(a1{I}==0);
        abc = Shift2{((size(Slope1,2)+1)/2)};
        temp_rgb = abc.*mask1;
        a1{I} = a1{I} + temp_rgb;
    end
    
p1 = 4;
p2 = -0;

% x = 0:abs(Slope(1,1)-Slope(1,2)):4.1;
x = 0:abs(Slope(1,1)-Slope(1,2)):Slope(1,end);
fx = double(round(p1*x + p2));

    for I=6:9
        fx1=double(round(fx/(10-I)));

        
        temp = round(a2{I}.*10);
        temp_opp = round(a2_opp{I}.*10);
        for i=(size(Slope1,2)+1)/2+1:size(Slope1,2)
%             abc=0;
            mask1 = (temp==Slope1(i));
            if (mask1==0)
                continue;
            else
                mask1 = mask1.*(a1{I}==0);
                mask_opp = temp_opp.*mask1;
                C = unique(mask_opp);
                abc = Shift2{i}(1:512-fx2(1,i-size(fx2,2)),1:512-fx1(1,i-size(fx1,2)),:);
                abc =(padarray(abc,[fx2(1,i-size(fx2,2)) fx1(1,i-size(fx1,2))], 'pre'));
                if isempty(C)
                    continue;
                else
                    for j=1:size(C,1)
                        if (C(j)==Slope1(i))
                            %continue;
                            mask1 = mask1.*(a1{I}==0);
                            temp_rgb = abc.*mask1;
                            a1{I} = a1{I} + temp_rgb;
                        else
                            mask2 = (temp_opp==C(j)).*mask1;
                            if (mask2==0)
                                continue;
                            else
                                X = find(C(j)==Slope1);
                                X1 = find(abs(C(j))==Slope1)-1;
                                abc1 = Shift2{X}(1:512-fx2(1,i-size(fx2,2)),1:512-fx1(1,i-size(fx1,2)),:);
                                abc1 =(padarray(abc1,[fx2(1,i-size(fx2,2)) fx1(1,i-size(fx1,2))], 'pre'));
                                C1(:,:,1) = conv2(abc1(:,:,1),PSF{X1*2});
                                C1(:,:,2) = conv2(abc1(:,:,2),PSF{X1*2});
                                C1(:,:,3) = conv2(abc1(:,:,3),PSF{X1*2});
                                size_C1 = ((size(C1,1) - size(abc1,1))+1)/2;
                                abc1 = C1(size_C1:size(C1,1)-size_C1,size_C1:size(C1,1)-size_C1,:);
                                clear C1;                   
                                temp_rgb = abc.*mask2 - (abc1./1.5).*mask2;
                                temp_rgb = imhistmatchn(temp_rgb,abc.*mask2);
                                temp_rgb = temp_rgb.*mask2;
                                a1{I} = a1{I} + temp_rgb;
                            end
                        end
                    end
                end
            end
        end
        %figure,imagesc(a1{I}),colormap(gray),axis off

        temp = round(a2{I}.*10);
        temp_opp = round(a2_opp{I}.*10);
        for i=1:(size(Slope1,2)+1)/2-1
%             abc=0;
            mask1 = (temp==Slope1(i));
            if (mask1==0)
                continue;
            else
                mask1 = mask1.*(a1{I}==0);
                mask_opp = temp_opp.*mask1;
                C = unique(mask_opp);
                abc = Shift2{i}(1+fx1(1,size(fx1,2)-i):512,1+fx1(1,size(fx1,2)-i):512,:);
                abc =(padarray(abc,[fx1(1,size(fx1,2)-i) fx1(1,size(fx1,2)-i)], 'post'));
                if isempty(C)
                    continue;
                else
                    for j=1:size(C,1)
                        if (C(j)==Slope1(i))
                            %continue;
                            mask1 = mask1.*(a1{I}==0);
                            temp_rgb = abc.*mask1;
                            a1{I} = a1{I} + temp_rgb;
                        else
                            mask2 = (temp_opp==C(j)).*mask1;
                            if (mask2==0)
                                continue;
                            else
                                X = find(C(j)==Slope1);
                                X1 = find(abs(C(j))==Slope1)-1;
                                abc1 = Shift2{X}(1+fx1(1,size(fx1,2)-i):512,1+fx1(1,size(fx1,2)-i):512,:);
                                abc1 =(padarray(abc1,[fx1(1,size(fx1,2)-i) fx1(1,size(fx1,2)-i)], 'post'));
                                C1(:,:,1) = conv2(abc1(:,:,1),PSF{X1*2});
                                C1(:,:,2) = conv2(abc1(:,:,2),PSF{X1*2});
                                C1(:,:,3) = conv2(abc1(:,:,3),PSF{X1*2});
                                size_C1 = ((size(C1,1) - size(abc1,1))+1)/2;
                                abc1 = C1(size_C1:size(C1,1)-size_C1,size_C1:size(C1,1)-size_C1,:);
                                clear C1;
                                temp_rgb = abc.*mask2 - (abc1./1.5).*mask2;
                                temp_rgb = imhistmatchn(temp_rgb,abc.*mask2);
                                temp_rgb = temp_rgb.*mask2;
                                a1{I} = a1{I} + temp_rgb;
                            end
                        end
                    end
                end
            end
        end
        temp = round(a2{I}.*10);
        mask1 = (temp==Slope1((size(Slope1,2)+1)/2));
        mask1 = mask1.*(a1{I}==0);
        abc = Shift2{((size(Slope1,2)+1)/2)};
        temp_rgb = abc.*mask1;
        a1{I} = a1{I} + temp_rgb;
        end
    end
    
end