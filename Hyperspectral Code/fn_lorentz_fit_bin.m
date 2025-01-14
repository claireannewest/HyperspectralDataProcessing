function [c,g,xbin]=fn_lorentz_fit_bin(x,y,fitnum,binfac)
%Lorentz fit
%Fits the spectrum to a sum of one, two, or three Lorentzian peaks
%b1 should be the peak maximum, c1 should be the FWHM

%One Lorentzian
if fitnum ==1
f_t=fittype('(2*a1/pi).*(c1./(4*(x-b1).^2+c1.^2))',...
         'dependent',{'y'},'independent',{'x'},...
     'coefficients',{'a1', 'b1', 'c1'});

 switch binfac
     case 'one'
         xbin=x;
         ybin=y;
     case 'two'
 
         for i = 1:numel(x)/2
             xbin(i,1)=(x(i*2-1)+x(i*2))/2;
             ybin(i,1)=y(i*2-1)+y(i*2);
         end
         
     case 'four'
         
         for i = 1:numel(x)/4
             xbin(i,1)=(x(i*4-1)+x(i*4-2))/2;
             ybin(i,1)=y(i*4-3)+y(i*4-2)+y(i*4-1)+y(i*4);
         end
 end
%parameters
fo_ = fitoptions('method','NonlinearLeastSquares',...
    'Lower',[-Inf -Inf -Inf],...
    'Upper',[Inf Inf Inf],'MaxIter',8000,'tolfun',1e-8,'tolx',1e-8);
st_ = [0.7 540 60];  %starting points
set(fo_,'Startpoint',st_);

[c,g] = fit(xbin,ybin,f_t,fo_);
end

%%%%% Change the limits if you use these multi-Lorentzian fits in the 
%%%%% future - LSS 021813 

%Two Lorentzian
if fitnum==2
    
    switch binfac
     case 'one'
         xbin=x;
         ybin=y;
     case 'two'
 
         for i = 1:numel(x)/2
             xbin(i,1)=(x(i*2-1)+x(i*2))/2;
             ybin(i,1)=y(i*2-1)+y(i*2);
         end
         
     case 'four'
         
         for i = 1:numel(x)/4
             xbin(i,1)=(x(i*4-1)+x(i*4-2))/2;
             ybin(i,1)=y(i*4-3)+y(i*4-2)+y(i*4-1)+y(i*4);
         end
 end
    
f_t=fittype('(2*a1/pi).*(c1./(4*(x-b1).^2+c1.^2))+(2*a2/pi).*(c2./(4*(x-b2).^2+c2.^2))');

%parameters
fo_ = fitoptions('method','NonlinearLeastSquares','Lower',[0 0 400 400 0 0],'Upper',[Inf Inf 800 1000 Inf Inf],'MaxIter',8000);
st_ = [4631 2491 565 692 58.54 62.5];  %starting points
set(fo_,'Startpoint',st_);

[c,g] = fit(x,y,f_t,fo_);
end

%Three Lorentzian
if fitnum ==3
f_t=fittype('(2*a1/pi).*(c1./(4*(x-b1).^2+c1.^2))+(2*a2/pi).*(c2./(4*(x-b2).^2+c2.^2))+(2*a3/pi).*(c3./(4*(x-b3).^2+c3.^2))');

%parameters
fo_ = fitoptions('method','NonlinearLeastSquares','Lower',[0 0 0 400 400 600 0 0 0],'Upper',[Inf Inf Inf 800 800 1200 Inf Inf Inf],'MaxIter',2000);
st_ = [0.5 0.5 0.5 600 650 850 111.6 100 100];  %starting points
set(fo_,'Startpoint',st_);

[c,g] = fit(x,y,f_t,fo_);
end