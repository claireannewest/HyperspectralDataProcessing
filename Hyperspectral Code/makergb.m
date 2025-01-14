function rgb=makergb(min,max,pix,gamma)
%clear all
%close all

intmax=1; %intensity max
%gamma=0.8;
%min=380; %min wvlth
%max=780; %max wvlth
%pix=1340; %total pixels


nmstep=(max-min)/(pix-1);
wvlth=min:nmstep:max;

for i=1:pix

    if((wvlth(i) >= 380) && (wvlth(i)<=440))
        red(i) = -(wvlth(i) - 440) / (440 - 380);
        green(i) = 0.0;
        blue(i) = 1.0;
    elseif((wvlth(i) >= 440) && (wvlth(i)<=490))
        red(i) = 0.0;
        green(i) = (wvlth(i) - 440) / (490 - 440);
        blue(i) = 1.0;
    elseif((wvlth(i) >= 490) && (wvlth(i)<=510))
        red(i) = 0.0;
        green(i) = 1.0;
        blue(i) = -(wvlth(i) - 510) / (510 - 490);
    elseif((wvlth(i) >= 510) && (wvlth(i)<=580))
        red(i) = (wvlth(i) - 510) / (580 - 510);
        green(i) = 1.0;
        blue(i) = 0.0;
    elseif((wvlth(i) >= 580) && (wvlth(i)<=645))
        red(i) = 1.0;
        green(i) = -(wvlth(i) - 645) / (645 - 580);
        blue(i) = 0.0;
    elseif((wvlth(i) >= 645) && (wvlth(i)<=780))
        red(i) = 1.0;
        green(i) = 0.0;
        blue(i) = 0.0;
        else
        red(i) = 0.0;
        green(i) = 0.0;
        blue(i) = 0.0;
    end
end

% figure
% plot(wvlth,red,'r',wvlth,green,'g',wvlth,blue,'b')    

for j=1:pix
    if((wvlth(j) >= 380) && (wvlth(j)<420))
        factor = 0.3 + 0.7*(wvlth(j) - 380) / (420 - 380);
    elseif((wvlth(j) >= 420) && (wvlth(j)<700))
        factor = 1.0;
    elseif((wvlth(j) >= 700) && (wvlth(j)<=780))
        factor = 0.3 + 0.7*(780 - wvlth(j)) / (780 - 700);
        else
        factor = 0.0;
    end
    if red(j)==0
        rgb(1,j)=0;
    else
    rgb(1,j)=(intmax*(red(j)*factor)^gamma);
    end
    if green(j)==0
        rgb(2,j)=0;
    else
    rgb(2,j)=(intmax*(green(j)*factor)^gamma);
    end
    if blue(j)==0
        rgb(3,j)=0;
    else
    rgb(3,j)=(intmax*(blue(j)*factor)^gamma);
    end
    
end

% figure
% plot(wvlth,rgb(1,:),'r',wvlth,rgb(2,:),'g',wvlth,rgb(3,:),'b')

for r=1:400
rgbim(r,:,1)=rgb(1,:);
rgbim(r,:,2)=rgb(2,:);
rgbim(r,:,3)=rgb(3,:);
end
% {
% figure
% subplot(2,1,1)
% imshow(rgbim)
% subplot(2,1,2)
% plot(wvlth,rgb(1,:),'r',wvlth,rgb(2,:),'g',wvlth,rgb(3,:),'b')
% 
% figure
% subaxis(2,1,1)
% imshow(rgbim)
% subaxis(2,1,2)
% plot(wvlth,rgb(1,:),'r',wvlth,rgb(2,:),'g',wvlth,rgb(3,:),'b')
% axis tight
% ylabel('Relative Intensity','fontsize',12,'fontweight','b')
% xlabel('Wavelength (nm)','fontsize',12,'fontweight','b')
% set(gca,'fontsize',12,'fontweight','b')
% }
    