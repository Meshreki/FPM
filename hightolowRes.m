function hightolowRes(O,P,dirac_cen,idx_led,Np,opts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file is to convert the high resolution image
% into low resolution images with the knowledge of their
% corresponding dirac peak positions.
%
% by John Meshreki, john.meshreki@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% helper functions
%Fourier and inverse Transforms
F = @(x) fftshift(fft2(ifftshift(x)));
Ft = @(x) fftshift(ifft2(ifftshift(x)));

% operator to crop region of O
downsamp = @(x,cen,Np) x(cen(1)-floor(Np(1)/2):cen(1)-floor(Np(1)/2)+Np(1)-1,...
    cen(2)-floor(Np(2)/2):cen(2)-floor(Np(2)/2)+Np(2)-1);

%% Specify dependencies

% add path to dependencies
addpath('../dependencies/natsortfiles');
addpath('../dependencies/export_fig');
addpath('../dependencies/labelpoints');
addpath('../dependencies/min_max_elements_index_n_values');

% add path for functions files
addpath('FP_Func/');

out_dir = opts.out_dir;


%% read in all images into the memory first
fprintf(['loading the high resolution image...\n']);
tic;

% high res complex field
hig_res_O = O; 

%% ROI, pupil, etc.
Np = Np;
P = P;
% No. images
Nimg=293;
% read led indices
idx_led = idx_led;
%read dirac centers
dirac_cen = dirac_cen;

toc;

%% saving the high res image for testing it is loaded correctly
imwrite(uint16(hig_res_O), strcat(out_dir,'high_res_O_image.tif'));
f1 = figure('visible','off');imshow(hig_res_O,[]);
title('(high res O)');
export_fig(f1,strcat(out_dir,'high_res_O_figure.tif'),'-m4');


%% go to Fourier space
hig_res_O_fourier = F(hig_res_O);


%% Estimate low res from high res

% all intensities for estimated low res
I_est_stack = zeros(Np(1),Np(2),Nimg);

for m = 1:Nimg

    % find Dirac peak position for image m (1 LED per image)
    index_m = find(idx_led==m);
    dirac_cen_pos = dirac_cen(index_m, :);
    

    % crope and apply the circular filter(P) at the dirac peak position 
    O_cropped = downsamp(hig_res_O_fourier,dirac_cen_pos,Np).*P;

    % go to real space
    O_est_low_res = Ft(O_cropped);

    % get estimated low res intensity from Object
    I_est = abs(O_est_low_res).^2;
    
    % save all estimated low res intensities in one array
    I_est_stack(:,:,m) = I_est;
    
    % saving them as TIFF images
    imwrite2tif(I_est, [], convertStringsToChars(strcat(opts.out_dir,'I_est_',num2str(m),'_image.tiff')), 'single','Compression',1);

    %% saving variables to matlab files (save every 10 images together!) 
    %    save(low_res_est_I_matfile, 'I_all_est_low_res', '-v7.3');
    %if rem(m,10)==0 | m==293
    %    fprintf('\nfinished m = %2d images\n', m);
    %    low_res_est_I_matfile = fullfile(strcat(out_dir,'low_res_est_I_till_m_',num2str(m),'.mat'));
    %    save(low_res_est_I_matfile, 'I_all_est_low_res', '-v7.3');
    %end
end    

I_est_stack_matfile = fullfile(strcat(out_dir,'I_est_stack','.mat'));
save(I_est_stack_matfile, 'I_est_stack', '-v7.3');
end
