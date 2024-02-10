function [a2] = Depth_fill(a2)


for k=1:4
        for i=2:size(a2{k},1)-1
            for j=2:size(a2{k},2)-1
%             for j=2:size(a2{k},2)-1
                if(isnan(a2{k}(i,j)))

                    temp1_d=a2{k}(i,j-1);

                    temp3_d=a2{k}(i-1,j);
                    temp4_d=a2{k}(i+1,j);
                    
                    if(~isnan(temp1_d) && ~isnan(temp3_d) && ~isnan(temp4_d))  
                            
                        temp=[temp1_d,temp3_d,temp4_d];

                        [temp_max,temp_index]=max(temp,[],'omitnan');

                            if(temp_index==1)
                                a2{k}(i,j)=temp1_d;

        
                            elseif(temp_index==2)
                                a2{k}(i,j)=temp3_d;

        
                            elseif(temp_index==3)
                                a2{k}(i,j)=temp4_d;


                            end
                         
                    else
                       temp_var=NaN;
                          
                       if(isnan(temp1_d) )
                           count3=1;
                          for nan=j-1:-1:2
                              
                              if(~isnan(a2{k}(i,nan)) || count3>9 )
                                  temp1_d=a2{k}(i,nan);
                                  break;
                              end
                              count3=count3+1;
                          end
                       end
                     
                       if(isnan(temp3_d) )
                           count4=1;
                           for nan=i-1:-1:2
                              if(~isnan(a2{k}(nan,j)) || count4>9 ) 
                                  temp3_d=a2{k}(nan,j);
                                  break;
                              end
                              count4=count4+1;
                          end
                       end 
                       
                       if(isnan(temp4_d))
                           count5=1;
                           for nan=i+1:1:size(a2{k},1)
                              if(~isnan(a2{k}(nan,j)) || count5>9 ) 
                                  temp4_d=a2{k}(nan,j);
                                  break;
                              end
                           count5=count5+1;
                          end
                           
                       end
                       
                           
                           temp=[temp1_d,temp3_d,temp4_d];

                           [temp_max,temp_index]=max(temp,[],'omitnan');

                           if(temp_index==1)
                                    a2{k}(i,j)=temp1_d;


                           elseif(temp_index==2)
                                    a2{k}(i,j)=temp3_d;


                           elseif(temp_index==3)
                                    a2{k}(i,j)=temp4_d;


                           end
                        
                    
                    end
                        
                 end
             end
        end
        %make main for loop go from 4:1
        %add a for loop here that adds these depth values from {K} to the next depth
        %map{k-1}
end
% end


for k=6:9
        for i=2:size(a2{k},1)-1
            for j=size(a2{k},2)-1:-1:2

                if(isnan(a2{k}(i,j)))

                    temp1_d=a2{k}(i,j+1);

                    temp3_d=a2{k}(i-1,j);
                    temp4_d=a2{k}(i+1,j);
                    
                    if(~isnan(temp1_d) && ~isnan(temp3_d) && ~isnan(temp4_d))  
                            
                        temp=[temp1_d,temp3_d,temp4_d];

                        [temp_max,temp_index]=max(temp,[],'omitnan');

    %                     temp1_rgb=a1{k}(i,j+1,:);
    %                     temp2_rgb=a1{k}(i,j-1,:);
                        
%                             a2{k}(i,j)=temp1_d;
                            if(temp_index==1)
                                a2{k}(i,j)=temp1_d;

        
                            elseif(temp_index==2)
                                a2{k}(i,j)=temp3_d;

        
                            elseif(temp_index==3)
                                a2{k}(i,j)=temp4_d;


                            end
                         
                    else
                       temp_var=NaN;
                          
                       if(isnan(temp1_d) )
                           count3=1;
                          for nan=j+1:1:size(a2{k},2)
                              
                              if(~isnan(a2{k}(i,nan)) || count3>9 )
                                  temp1_d=a2{k}(i,nan);
                                  break;
                              end
                              count3=count3+1;
                          end
                       end
                     
                       if(isnan(temp3_d) )
                           count4=1;
                           for nan=i-1:-1:2
                              if(~isnan(a2{k}(nan,j)) || count4>9 ) 
                                  temp3_d=a2{k}(nan,j);
                                  break;
                              end
                              count4=count4+1;
                          end
                       end 
                       
                       if(isnan(temp4_d))
                           count5=1;
                           for nan=i+1:1:size(a2{k},1)
                              if(~isnan(a2{k}(nan,j)) || count5>9 ) 
                                  temp4_d=a2{k}(nan,j);
                                  break;
                              end
                           count5=count5+1;
                          end
                           
                       end
                       
                           
                           temp=[temp1_d,temp3_d,temp4_d];

                           [temp_max,temp_index]=max(temp,[],'omitnan');

                           if(temp_index==1)
                                    a2{k}(i,j)=temp1_d;


                           elseif(temp_index==2)
                                    a2{k}(i,j)=temp3_d;


                           elseif(temp_index==3)
                                    a2{k}(i,j)=temp4_d;


                           end
                        
                    
                    end
                        
                 end
             end
        end
        %make main for loop go from 6:9
        %add a for loop here that adds these depth values from {K} to the next depth
        %map{k+1}
        
end


% for i=1:9
%     a2{i}=(medfilt2(a2{i}, [ 3 3 ]));
%     a2{i}=round(a2{i});
%     mask=isnan(a2{i}); 
%     a2{i}=round(regionfill(a2{i},mask));
% end
    
    
end

