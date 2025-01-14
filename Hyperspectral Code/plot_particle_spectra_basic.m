%% Plots the spectra of particles at specified pixels
close all 
clear all
clc
%% Parameters

% Specify path to folder with analysis data
path = 'data_from_2024-12-21/';
datafile = 'mg110_glyc_33';
sample = [datafile, 'analysis'];
% cd(strcat(path));


% Specify file name of data you want to plot
load(sample, 'specfin', 'wvlths')

place2save = [path, datafile, '/'];
% if ~exist(place2save, 'dir')
%    mkdir(place2save)
% end


% Specify pixel coordinates
part_coord = readmatrix([place2save, 'positions.txt']);
part_coord_flip = circshift(part_coord,[0,1]);

 
% Plot hyperspectral data cube
% figure1 = figure;
% specfin_final = sum(specfin(:,:,100:end),3); % Plot the spectrum image, summing over wavelengths excluding higher energy
% imshow(specfin_final,[0 100])
% colormap('turbo');
% 
% hold all
% for n = 1:size(part_coord_flip,1)
%     text(part_coord_flip(n,2),part_coord_flip(n,1),num2str(n),'Color','White', 'FontSize', 24)
% end
% 
% figure1.Position(3:4)=[800,1500];
% 
% hgsave([place2save, 'hyperspec_image'],'-v6')
% close;

%% Specify cutoff bounds in units of pixels, not nm. 

lower_bound = 100;
upper_bound = 670;
rawwvlths = wvlths(100:end);
wvlths = wvlths(lower_bound:upper_bound);

% Describes the number of pixels you want to integrate/particle. Typically 3 or 5.
int_size = 3;                                              
int_var_low = -1*(int_size-1)./2;
int_var_high = (int_size-1)./2;


%% Integrates over the specified square and then plots each specified pixel.
cd(place2save)

numPart = size(part_coord,1);
numWave_raw = 670;
allSpec = zeros(571, numPart);

for n = 1:numPart
    % Full data background correcting
    x_part = part_coord(n,2);
    y_part = part_coord(n,1);

    spec_part_n = zeros(1, 1, numWave_raw);
    spec_part_n = reshape(spec_part_n, [numWave_raw,1]);

    for m = int_var_low:int_var_high
        for l = int_var_low:int_var_high
            a = x_part + m;
            b = y_part + l;
            spec_part_n = spec_part_n + reshape(specfin(a,b,:), [numWave_raw,1]);    

         end
    end

    part_spec = spec_part_n(lower_bound:upper_bound);
       
    allSpec(:, n) = part_spec;
    
    % [param_1, param_2] = fn_lorentz_fit(wvlths', part_spec, 1, 1);
    % a1 = param_1.a1;
    % b1 = param_1.b1;
    % c1 = param_1.c1;
    % resonance = b1;
    % FWHM = c1;
    % r_list = param_2.rsquare;
    % lorentz_fit =(2*a1/pi).*(c1./(4*(wvlths'-b1).^2+c1.^2));
    % 
    % Noi = std(part_spec-lorentz_fit);
    % [Notneeded, IndiMax]=min(abs(wvlths-resonance));
    % Sigy = part_spec(IndiMax);
    % SnN= Sigy/Noi;    
    % 
    % figure1 = figure;
    % hold all
    % plot(rawwvlths, part_spec,'b','linewidth',3)
    % plot(wvlths, lorentz_fit,'k--','linewidth',3)
    % 
    % 
    % xlabel('Wavelength (nm)','fontsize',32)
    % ylabel('Scattering (unitless)','fontsize',32)
    % set(gca,'FontSize',22,'box','on')
    % text(0.1,0.9,['\lambda_m_a_x = ',num2str(round(resonance)), ' nm'],'fontsize',20,'Units','normalized')
    % text(0.1,0.78,['\Gamma = ',num2str(round(FWHM)), ' nm'],'fontsize',20,'Units','normalized')
    % text(0.1,0.66,['S/N = ',num2str(round(SnN))],'fontsize',20,'Units','normalized')
    % xlim([450 950])
    % ylim([0 1.1])
    % 
    % title(['NP: ', num2str(n)])
    % 
    % saveas(figure1,[datafile, '_', num2str(n),'.tif'])
    % close;
    
end
save('all_spectra_mg110_glyc_33','rawwvlths','allSpec')

