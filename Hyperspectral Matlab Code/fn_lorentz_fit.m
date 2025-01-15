function [c,g]=fn_lorentz_fit(x,y,guess,fitnum,lambda_1,lambda_2)
%Lorentz fit
%Fits the spectrum to a sum of one, two, or three Lorentzian peaks
%b1 should be the peak maximum, c1 should be the FWHM

%One Lorentzian
if fitnum ==1
f_t=fittype('(2*a1/pi).*(c1./(4*(x-b1).^2+c1.^2))',...
         'dependent',{'y'},'independent',{'x'},...
     'coefficients',{'a1', 'b1', 'c1'});

%parameters
fo_ = fitoptions('method','NonlinearLeastSquares',...
    'Lower',[0 min(x) 0],...
    'Upper',[Inf max(x) Inf], ...
    'StartPoint',[100 guess 70],...
    'MaxIter',8000,'tolfun',1e-6,'tolx',1e-6);
% 
% st_ = [100 guess 70];  %starting points
% set(fo_,'Startpoint',st_);

[c,g] = fit(x,y,f_t,fo_);
end

