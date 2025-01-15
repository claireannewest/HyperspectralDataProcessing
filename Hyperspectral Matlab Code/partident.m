function [ptu mm mark2] = partident(specrgbnorm,lower,upper,nhood)

srgb=sum(specrgbnorm,3,'double');
[x,y]=size(srgb);
pts=[];

for i=lower:0.001:upper
rbs=im2bw(specrgbnorm,i); %figure; imshow(rbs,'initialmagnification','fit'); title(num2str(i))
rbs=imclearborder(rbs);
mask1=ones(x-(nhood-1),y-(nhood-1));
mask2=padarray(mask1,[(nhood-1)/2,(nhood-1)/2],0,'both');
rbs=rbs.*mask2;
CC = bwconncomp(rbs);
numPixels = cellfun(@numel,CC.PixelIdxList);
[biggest,idx] = max(numPixels);
while biggest > 1
rbs(CC.PixelIdxList{idx}) = 0;
CC = bwconncomp(rbs);
numPixels = cellfun(@numel,CC.PixelIdxList);
[biggest,idx] = max(numPixels);
end
% going to multiply by a picture frame here so that I can't select a
% particle on the edge ever

pts=[pts;find(rbs)];

end

%%
ptu=unique(pts);
[mark2(:,2),mark2(:,1)]=ind2sub(size(srgb),ptu);
[mm,~]=size(mark2);
%{
figure; subplot(1,2,1)
imshow(specrgbnorm,'initialmagnification','fit')
subplot(1,2,2)
imshow(specrgbnorm,'initialmagnification','fit')

hold on
    
    for ki = 1:mm
        
        plot(mark2(ki,1),mark2(ki,2),'go','MarkerSize',5,'LineWidth',1)
        hold on
    end
    hold off
    title(num2str(mm))

%}
   
