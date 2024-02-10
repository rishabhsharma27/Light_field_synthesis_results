function [g,g1, w_size2 ,w_width2]= octagon_gen(w_size1,w_width1)

g=zeros(w_size1,w_size1);
% w_size1=15;
% w_width1=1;             %vertical ansd horizontal width of the octagon on either sides 
                    %for the total width 5 , w_width1 should be +2 , -2 i.e.2
                    
     mid1=round(w_size1/2);
     w_size2= w_size1-(w_width1*2+1);
     mid2=round(w_size2/2);
     w_width2=1;
 if (mod(w_size1,2)==1) 
     
%      mid1=round(w_size1/2);
%      w_size2= w_size1-w_width1*2+1;
%      mid2=round(w_size2/2);
%      w_width2=2;

     for i=1:w_size1
        for j=1:w_size1
          if( i <=mid1+w_width1)  
             if( (j>=mid1+w_width1+i || j<=mid1-w_width1-i))  %if window is even use (j<=mid1-w_width1-i+1))
                    g(i,j)=0;
             else
                 g(i,j)=1;
             end
          elseif (i>mid1+w_width1)
              if(  (j>=w_size1-i+mid1+w_width1+1 || j<=i-mid1-w_width1))
                   g(i,j)=0;
              else 
                  g(i,j)=1;
              end
          end
        end
     end
     
     for i=1:w_size2
        for j=1:w_size2
          if( i <=mid2+w_width2)  
             if( (j>=mid2+w_width2+i || j<=mid2-w_width2-i+1))  %if window is even use (j<=mid2-w_width2-i+1))
                    g1(i,j)=0;
             else
                 g1(i,j)=1;
             end
          elseif (i>mid2+w_width2)
              if(  (j>=w_size2-i+mid2+w_width2+1 || j<=i-mid2-w_width2))
                   g1(i,j)=0;
              else 
                  g1(i,j)=1;
              end
          end
        end
     end

    % g=[g g;g g];
    % figure
    % imshow(g)
 else 
     
     mid1=round(w_size1/2); 
     w_size2= w_size1-w_width1*2;
     mid2=round(w_size2/2);
     for i=1:w_size1
        for j=1:w_size1
          if( i <=mid1+w_width1)  
             if( (j>=mid1+w_width1+i || j<=mid1-w_width1-i+1))  %if window is even use (j<=mid1-w_width1-i+1))
                    g(i,j)=0;
             else
                 g(i,j)=1;
             end
          elseif (i>mid1+w_width1)
              if(  (j>=w_size1-i+mid1+w_width1+1 || j<=i-mid1-w_width1))
                   g(i,j)=0;
              else 
                  g(i,j)=1;
              end
          end
        end
     end
     for i=1:w_size2
        for j=1:w_size2
          if( i <=mid2+w_width2)  
             if( (j>=mid2+w_width2+i || j<=mid2-w_width2-i+1))  %if window is even use (j<=mid2-w_width2-i+1))
                    g1(i,j)=0;
             else
                 g1(i,j)=1;
             end
          elseif (i>mid2+w_width2)
              if(  (j>=w_size2-i+mid2+w_width2+1 || j<=i-mid2-w_width2))
                   g1(i,j)=0;
              else 
                  g1(i,j)=1;
              end
          end
        end
     end
     
 end
 