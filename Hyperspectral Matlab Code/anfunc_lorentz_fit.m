function anfunc_lorentz_fit(dataloc,names,c3,stanbig,backper,lowercut,uppercut,lower,upper,nhood,rsquarelim,binfac)

%anfunc_lorentz_fit is a subfunction of hyper_analysis which processes
%hyperspectral images. This function is used primarily for memory
%management issue and passes no output back to the its parent function. It
%instead writes data to files named as the original TDMS filename followed
%by analysis.mat
%function dependencies: convertTDMS.m makergb.m partident.m
%fn_lorentz_fit.m plot_lorentz_fit.m 
%As of 08-21-13, framework for local background subtraction is provided but
%not verified. Some level of assurance must be reached to identify and
%exclude other particles, perhaps 2D gaussian fitting to a sum image, or a
%local version of the existing global subtraction.
    che=strcat(dataloc,'/',names{1,1}(c3+2),'.mat');
    che2=strcat(dataloc,'/',names{1,1}(c3+2),'analysis','.mat');
    fch=exist(che{1},'file');
    if fch==2
        load(che{1})
    else
        tdmsfile = strcat(names{1,1}(c3+2),'.tdms');
        [ConvertedData,~,~,~,~] = convertTDMS(1,tdmsfile,dataloc);   % TDMS to Matlab converter
    end
    grpnum = ConvertedData.Data.Root.Property(1,7).Value;   %wavelength grouping
    top = ConvertedData.Data.Root.Property(1,9).Value;      %top pixel used for binning
    bottom = ConvertedData.Data.Root.Property(1,10).Value;  %bottom pixel used for binning
    prow = ConvertedData.Data.MeasuredData(1,1).Property(1,10).Value-ConvertedData.Data.MeasuredData(1,1).Property(1,9).Value+1;
    pcol = ConvertedData.Data.MeasuredData(1,1).Property(1,8).Value;
    ncol = length(ConvertedData.Data.MeasuredData(4).Data); % number of columns in future array (CCD spanning wavelengths)
    npdata = zeros(prow*pcol,ncol);                        
    for c2 = 3:prow*pcol+2   % TDMS includes two empty structures for names/properties; data starts in third structure
        npdata(c2-2,:) = (ConvertedData.Data.MeasuredData(c2).Data');     % storing spectral data for one pixel as a row in matrix (turned horizontally)
    end
    wvlths=ConvertedData.Data.MeasuredData(end-2).Data';
    rgbspec=makergb(wvlths(1),wvlths(end),ncol,0.8);
    hypr=rgbspec(1,:);
    hypg=rgbspec(2,:);
    hypb=rgbspec(3,:);
    clear rgbspec
    sumimg = (zeros(prow,pcol));            % sum image matrix
    specim = (zeros(prow,pcol,ncol));       % 3D data cube
    bkavg = zeros(1,ncol);                  % background
    stanim = zeros(prow,pcol,ncol);         % 3D standard
    redim = zeros(prow,pcol,ncol);          % calculate R
    greenim = zeros(prow,pcol,ncol);        % calculate G
    blueim = zeros(prow,pcol,ncol);         % calculate B
    bkgim = zeros(prow,pcol,ncol);          % 3D global background
      if size(stanbig,3) == 1340
        grpnum = 2;
    else
        grpnum = 1;
    end
    for c4 = 1:pcol         % for all columns in the image
        for c5 = 1:prow     % for all rows in the image
            sumimg(c5,c4) = (sum(npdata(c5 + (c4-1)*prow,:)));  % summed data of pixel to produce image
            redim(c5,c4,:)= (hypr);
            greenim(c5,c4,:)= (hypg);
            blueim(c5,c4,:)= (hypb);
            specim(c5,c4,:) = ConvertedData.Data.MeasuredData(1,c5+(c4-1)*prow+2).Data';
        end
        stanim(:,c4,:)=stanbig(top:bottom,1,1:grpnum:end); %selects the local standard based on position on CCD
    end
    clear ConvertedData npdata
    
    smln=ceil(backper*prow*pcol);   % number of pixels to be included in background
%     imnorm=specim./stanim;          % Flatfield correct image first for summing
    imnorm=specim./1;          % Flatfield correct image first for summing
    sumnorm=sum(imnorm,3);          % creates sum image over all wavelengths
    sumnorm=(sumnorm-(min(min(sumnorm))))/(max(max(sumnorm-(min(min(sumnorm))))));  %normalize over this image
    sorted=unique(sumnorm);         % find the lowest index pixels
    nthsmlst=sorted(smln+1);        % find the cutoff value for background
    [a,b]=ind2sub(size(sumnorm), find(sumnorm<nthsmlst,smln)); %Record the backgound positions
    
    clear imnorm
    
    for d1 = 1:smln
        bkavg=bkavg+reshape(specim(a(d1),b(d1),:),1,ncol); % calculate global background sum
    end
    bkavg = bkavg/smln;             %global background average
    
    
    for c4 = 1:pcol                 % for all columns in the image
        for c5 = 1:prow             % for all rows in the image
            bkgim(c5,c4,:) = bkavg; % makes a background cube
        end
    end
    
    specfin = (specim - bkgim) ./ stanim;       % spectral correction
    % specfin = stanim;       % spectral correction


    clear bkgim bkgim stanim
  
    
   save(che2{1},'specfin','wvlths');
   
 
end