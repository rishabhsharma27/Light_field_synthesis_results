%if mistake, check https://github.com/rishabhsharma27/Light-field-refocusing

function [ImgOut_med] = LFFiltShiftSum_rish_mod_freq_domain_individ_subap( LF,LF_sep, Slope )
%#codegen
% 
% coder.extrinsic('fft2');
% coder.extrinsic('squeeze');

% coder.varsize('LFOut',[9 9 512 512 3]);


% coder.varsize('v',[1 512]);
% coder.varsize('u',[1 512]);
% 
% coder.varsize('IF_image',[512 512 3]);

%%

% FiltOptions = struct('Precision', 'single','Aspect4D', [1,1,1,1],'FlattenMethod', 'sum','InterpMethod', 'linear',...
%     'ExtrapVal', 0,'MaskThresh', 0.5,'Mask', false,'UpsampRate', 1,...
%     'MinWeight', 1.192092895507813e-06,'Normalize', 'true', 'ExtrapMethod', 'nearest');

FiltOptions = struct('UpsampRate', 1);



LFSize = size(LF);
NColChans = size(LF,5);
% HasWeight = ( NColChans == 4 || NColChans == 2 );
% if( HasWeight )
% 	NColChans = NColChans-1;
% end


% TVSlope = Slope * FiltOptions.Aspect4D(3) / FiltOptions.Aspect4D(1);
% SUSlope = Slope * FiltOptions.Aspect4D(4) / FiltOptions.Aspect4D(2);

TVSlope = Slope ;
SUSlope = Slope ;

%%
v = linspace(1,LFSize(3), round(LFSize(3)*FiltOptions.UpsampRate));
u = linspace(1,LFSize(4), round(LFSize(4)*FiltOptions.UpsampRate));


%%


NewLFSize = LFSize;
NewLFSize(3:4) = [length(v), length(u)];



%%

VOffsetVec=zeros(1,LFSize(1));

UOffsetVec=zeros(1,LFSize(2));

for i=1:LFSize(1)
    VOffsetVec(i)=(i-1-floor(LFSize(1)/2)).*TVSlope;
    UOffsetVec(i)=(i-1-floor(LFSize(2)/2)).*SUSlope;
end


%%


%     LFOut = zeros(81,512,512,3);
    LFOut = zeros(NewLFSize, 'like', LF);


    [a,b,~]=size(LF_sep{1,1});
    Nr = ifftshift([-fix(a/2):ceil(a/2)-1]);
    Nc = ifftshift([-fix(b/2):ceil(b/2)-1]);
    [xF,yF] = meshgrid(Nr,Nc);
    
   
          
for( TIdx = 1:LFSize(1) )
	VOffset = VOffsetVec(TIdx)*FiltOptions.UpsampRate;
    
    for( SIdx = 1:LFSize(2) )
		UOffset = UOffsetVec(SIdx)*FiltOptions.UpsampRate;
        
%           in= imresize(CurSlice,FiltOptions.UpsampRate);  
          in=LF_sep{TIdx,SIdx}; %// Define input signal
          
%           in= imresize(in,FiltOptions.UpsampRate);  %uncomment if using UpsampRate

            x0=-UOffset; %// Define shifts
            y0=-VOffset;
            phase = 2;
            

%             H=exp(1i*2*pi*(x0*xF/a+y0*yF/b)).*exp(-1i*phase);
            H=exp(1i*2*pi*(x0*xF/a+y0*yF/b));
            IF_image = real(ifft2(in.*H));
            
            x0=round(x0);
            y0=round(y0);
            if (x0>0)
               IF_image(:,end-abs(x0)+1:end,:)=0;
            elseif(x0<0)
                IF_image(:,1:abs(x0),:)=0;
            end
            
            if (y0>0)
               IF_image(end-abs(y0)+1:end,:,:)=0;
            elseif(y0<0)
                IF_image(1:abs(y0),:,:)=0;
            end
%             LFOut(TIdx,SIdx, :,:, :)=IF_image*-2.4033;
                LFOut(TIdx,SIdx, :,:, :)=IF_image;
%             r1=IF_image;
%            r=(SIdx-1)*9+1+(TIdx-1);
% 		r2(r, :,:, :) = (IF_image);
        
    end
end

        LF = LFOut;

		t = reshape(LF(:,:,:,:,1:NColChans), [prod(LFSize(1:2)), NewLFSize(3:4), NColChans]);
		ImgOut_med = squeeze(median(t,"omit")); %if issue check github for code
        


end

% TimeStamp = datestr(now,'ddmmmyyyy_HHMMSS');
% FiltOptions.FilterInfo = struct('mfilename', mfilename, 'time', TimeStamp, 'VersionStr', LFToolboxVersion);


%% previous could, we don't need the weight dimension in LF if we use the median instead of the summ method

% function [ImgOut_med] = LFFiltShiftSum_rish_mod_freq_domain_individ_subap( LF,LF_sep, Slope )
% %#codegen
% 
% LFSize = size(LF);
% 
% FiltOptions = struct('Precision', 'single','Aspect4D', [1,1,1,1],'FlattenMethod', 'sum','InterpMethod', 'linear',...
%     'ExtrapVal', 0,'MaskThresh', 0.5,'Mask', false,'UpsampRate', 1,...
%     'MinWeight', 1.192092895507813e-06,'Normalize', 'true', 'ExtrapMethod', 'nearest');
% 
% 
% 
% LFSize = size(LF);
% NColChans = size(LF,5);
% HasWeight = ( NColChans == 4 || NColChans == 2 );
% if( HasWeight )
% 	NColChans = NColChans-1;
% end
% 
% 
% TVSlope = Slope * FiltOptions.Aspect4D(3) / FiltOptions.Aspect4D(1);
% SUSlope = Slope * FiltOptions.Aspect4D(4) / FiltOptions.Aspect4D(2);
% 
% %%
% v = linspace(1,LFSize(3), round(LFSize(3)*FiltOptions.UpsampRate));
% u = linspace(1,LFSize(4), round(LFSize(4)*FiltOptions.UpsampRate));
% 
% 
% %%
% 
% 
% NewLFSize = LFSize;
% NewLFSize(3:4) = [length(v), length(u)];
% 
% 
% 
% %%
% 
% VOffsetVec=zeros(1,LFSize(1));
% 
% UOffsetVec=zeros(1,LFSize(2));
% 
% for i=1:LFSize(1)
%     VOffsetVec(i)=(i-1-floor(LFSize(1)/2)).*TVSlope;
%     UOffsetVec(i)=(i-1-floor(LFSize(2)/2)).*SUSlope;
% end
% 
% 
% %%
% 
% 
% LFOut = zeros(NewLFSize, 'like', LF);
% LFOut(:,:,:,:,4)=LF(:,:,:,:,4);
% IF_image = complex(zeros(512,512,3));
% 
% for( TIdx = 1:LFSize(1) )
% 	VOffset = VOffsetVec(TIdx)*FiltOptions.UpsampRate;
%     for( SIdx = 1:LFSize(2) )
% 		UOffset = UOffsetVec(SIdx)*FiltOptions.UpsampRate;
%         
% %           in= imresize(CurSlice,FiltOptions.UpsampRate);  
%           in=LF_sep{TIdx,SIdx}; %// Define input signal
%           in= imresize(in,FiltOptions.UpsampRate); 
% %             H=fftshift(fft2(in)); %// Compute 2D Fourier Transform
%             x0=-UOffset; %// Define shifts
%             y0=-VOffset;
%             phase = 2;
%             [a,b,c]=size(in);
%             Nr = ifftshift([-fix(a/2):ceil(a/2)-1]);
%             Nc = ifftshift([-fix(b/2):ceil(b/2)-1]);
%             [xF,yF] = meshgrid(Nr,Nc);
%             IF_image = ifft2(fft2(in).*exp(1i*2*pi*(x0*xF/a+y0*yF/b))).*exp(-1i*phase);
% 
%             x0=round(x0);
%             y0=round(y0);
%             if (x0>0)
%                IF_image(:,end-abs(x0)+1:end,:)=0;
%             elseif(x0<0)
%                 IF_image(:,1:abs(x0),:)=0;
%             end
%             
%             if (y0>0)
%                IF_image(end-abs(y0)+1:end,:,:)=0;
%             elseif(y0<0)
%                 IF_image(1:abs(y0),:,:)=0;
%             end
%             CurSlice=real(IF_image)*-2.4033;
%         
%         
%         
% 		LFOut(TIdx,SIdx, :,:, 1:3) = CurSlice;
% 	end
% % 	if( mod(TIdx, ceil(LFSize(1)/10)) == 0 )
% % 		fprintf('.');
% % 	end
% end
% 
% LF = LFOut;
% 
% % if( FiltOptions.Normalize )
% % 	if( HasWeight )
% % 		for( iColChan = 1:NColChans )
% % 			LF(:,:,:,:,iColChan) = LF(:,:,:,:,iColChan) .* LF(:,:,:,:,end);
% % 		end
% % 	else % add a weight channel
% % 		LF(:,:,:,:,4) = ones(size(LF(:,:,:,:,1)), FiltOptions.Precision);
% % 		LFSize = size(LF);
% % 	end
% % end
% 
% clear LFOut
% %%
% 
% 
% if( FiltOptions.Normalize )
% 	W = LF(:,:,:,:,end);
% 	W(isnan(W)) = 0;
% 	LF(:,:,:,:,end) = W;
% end
% 
% 
% % 		ImgOut_sum = squeeze(nansum(nansum(LF,1),2));
% 
% 		t = reshape(LF(:,:,:,:,1:NColChans), [prod(LFSize(1:2)), NewLFSize(3:4), NColChans]);
% 		ImgOut_med = squeeze(nanmedian(t));
% % ImgOut_sum=ImgOut_med;
% 
% % %---
% % if( FiltOptions.Normalize )
% % 	WeightChan = ImgOut_sum(:,:,end);
% % 	InvalidIdx = find(WeightChan < FiltOptions.MinWeight);
% % 	ChanSize = numel(ImgOut_sum(:,:,1));
% % 	for( iColChan = 1:NColChans )
% % 		ImgOut_sum(:,:,iColChan) = ImgOut_sum(:,:,iColChan) ./ WeightChan;
% % 		ImgOut_sum( InvalidIdx + ChanSize.*(iColChan-1) ) = 0;
% % 	end
% % end
% % 
% %     ImgOut_sum(1:16,:,1:3)= ImgOut_med(1:16,:,1:3);
% %     ImgOut_sum(:,1:16,1:3)= ImgOut_med(:,1:16,1:3);
% %     ImgOut_sum(end-16+1:end,:,1:3)= ImgOut_med(end-16+1:end,:,1:3);
% %     ImgOut_sum(:,end-16+1:end,1:3)= ImgOut_med(:,end-16+1:end,1:3);
% 
% end
% 
% % TimeStamp = datestr(now,'ddmmmyyyy_HHMMSS');
% % FiltOptions.FilterInfo = struct('mfilename', mfilename, 'time', TimeStamp, 'VersionStr', LFToolboxVersion);
% 
% 
% 
% 
% 
