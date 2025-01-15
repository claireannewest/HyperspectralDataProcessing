
function stan = standardreadindark(tdmslamp,tdmsdark,dataloc)

%generates a standard to 
%written by Chad Byers 3-26-13


tdmsfile=tdmslamp;
[ConvertedData,~,~,~,~] = convertTDMS(0,tdmsfile,dataloc);    % TDMS to Matlab converter
prow = ConvertedData.Data.MeasuredData(1,1).Property(1,10).Value-ConvertedData.Data.MeasuredData(1,1).Property(1,9).Value+1;
pcol = ConvertedData.Data.MeasuredData(1,1).Property(1,8).Value;
ncol = length(ConvertedData.Data.MeasuredData(4).Data);     % number of columns in future array (CCD spanning wavelengths)
stanim = zeros(prow,pcol,ncol);


for c4 = 1:pcol         % for all columns in the image
    for c5 = 1:prow     % for all rows in the image
        stanim(c5,c4,:) = ConvertedData.Data.MeasuredData(1,c5+(c4-1)*prow+2).Data';
    end
end

stan =sum(stanim,2)/pcol;
%size(stan)

clear ConvertedData stanim

tdmsfile=tdmsdark;
[ConvertedData,~,~,~,~] = convertTDMS(0,tdmsfile,dataloc);    % TDMS to Matlab converter
prow = ConvertedData.Data.MeasuredData(1,1).Property(1,10).Value-ConvertedData.Data.MeasuredData(1,1).Property(1,9).Value+1;
pcol = ConvertedData.Data.MeasuredData(1,1).Property(1,8).Value;
ncol = length(ConvertedData.Data.MeasuredData(4).Data);     % number of columns in future array (CCD spanning wavelengths)
stanim = zeros(prow,pcol,ncol);


for c4 = 1:pcol         % for all columns in the image
    for c5 = 1:prow     % for all rows in the image
        stanim(c5,c4,:) = ConvertedData.Data.MeasuredData(1,c5+(c4-1)*prow+2).Data';
    end
end

dark =sum(stanim,2)/pcol;
a = squeeze(stan);
b = squeeze(dark);
c = a-b;
imshow(c)
size(dark)
stan = stan - dark;
fname=[dataloc '\standark.mat'];
save(fname,'stan')