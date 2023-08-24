%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% main file to implement Fourier Ptychography reconstruction algorithm
% ref
% Lei Tian, et.al, Biomedical Optics Express 5, 2376-2389 (2014).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All main.m, AtlerMin.m, GDUpdate_Multiplication_rank1.m, USAF_Parameter.m, 
% hela_Parameter.m files are modified and extended to:
%
% 1- allowing a FULL reconstruction of the input images for 
% any dimesnions (i.e. not only square crops of them!).
% 2- saving images (and figures) for multiple variables.
% 3- saving dirac peaks positions.
% 4- adding more descriptive comments (in addition to those added by Ivo Ihrke).
% 4- investigating the huge error of the Algorithm to be due to pixels w/ 
% very large (intensity) values in the input stacks (and not due to
% the Algorithm itself!).
%
% last modified on 27.05.2022
% by John Meshreki, john.meshreki@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% last modified on 10/07/2015
% by Lei Tian, lei_tian@alum.mit.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% To do list for the user: (marked by 'TODO#')
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1) specify where the data located in 'filedir'
% 2) specify where you want to store the results in 'out_dir'
% 3) Find coordinates for estimating background levels and input into 'bk1'
% and 'bk2'.
% 4) specify a threshold value, above which the estimated background value
% is signal rather than noise.
% 5) make sure the LED index (used for taking the images) are properly defined
% in 'lit'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% use today's date to create new output directory
todaysdatetime = string(datetime('now','Format','yyyy_MM_dd_HH_mm_ss'));
todaysdate = string(datetime('today','Format','yyyy_MM_dd'));



%% Reconstruction library locates here
%clear all;
addpath('../dependencies/natsortfiles');
addpath('../dependencies/export_fig');
addpath('../dependencies/labelpoints');
addpath('../dependencies/min_max_elements_index_n_values');
addpath('../dependencies/normalize_img');
addpath('../dependencies/imwrite2tif');


% add path for functions files
addpath('FP_Func/');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TODO 1: specify the file directory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% specify which sample to use as input
sample_name= 'USAF'; %'stained', 'USAF', 'hela'

% map container of input directory path to
% multiplex reading of images
input_dir_name = containers.Map({'stained'; 'USAF'; 'hela'},...
{'../data/Tian14/1LED/tif/';...
   %'../data/Tian14_ResTarget/1LED/';...
   '../data/CSE/2022_11_07/2022_11_07_13_34_41/1LED/';...
   '../data/Tian15_inVitroHeLa/data/'});

filedir = input_dir_name(sample_name);

% Generate the image list, in 'tif' image format (depending on your image format)
imglist = dir([filedir,'*.tif']);
%imglist = dir([filedir,'*.tiff']);
nstart = [100, 100];
N = natsortfiles({imglist.name});%sorting the images in 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TODO 2: specify output folder
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% map container of output directory path to
out_dir_name = containers.Map({'stained'; 'USAF'; 'hela'},...
{strcat('../out_dir/Tian14_StainedHistologySlide/',todaysdate,'/',todaysdatetime,'/');...
    strcat('/tmp/fpm/FPM/out_dir/Tian14_ResTarget/',todaysdate,'/',todaysdatetime,'/');...
    strcat('../out_dir/Out_Tian15_inVitroHeLa/',todaysdate,'/',todaysdatetime,'/')});

out_dir = out_dir_name(sample_name);
mkdir(out_dir);

% keep a log
diary(strcat(out_dir,'/','log_',todaysdatetime,'.txt'));


%% define # of LEDs used to capture each image
numlit = 1;
% raw image size
n1 = 2160; n2 = 2560;
%n1 = 3753; n2 = 5634;


%% read in all images into the memory first
fprintf(['loading the images...\n']);
tic;
Nimg = length(imglist);
fprintf('no. of images loaded:%2d \n', Nimg);
Iall = zeros(n1,n2,Nimg,'uint16'); %3d array; each pixel in the 2d image for all images
Ibk = zeros(Nimg,1);
for m = 1:Nimg

    %     fn = [filedir,imglist(m).name]; %removed for sorting purpose
    fn = [filedir,N{m}];  %added for sorting purpose
    disp(fn);
    % all image data
    I = double(imread(fn));
    %I = normalize_img(I); % normalize measured image
    Iall(:,:,m) = I  ; %filling 3d array by replacing zeros in 1st&2nd elements by the image values and 3rd element is index of image
    
    % creating cropped input files for Onofre
    I_input_stack_cropped(:,:, m) = double(I(1080:1080+100-1,1280:1280+150-1));
    imwrite(I_input_stack_cropped(:,:, m), strcat(out_dir,'I_input_stack_cropped','_image_',sprintf('%04d',m),'.png'));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% TODO 3: specify background region coordinates
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % backgr4ound noise esimtation 
    bk1 = mean2(double(Iall(1:100,1:100,m)));
    bk2 = mean2(double(Iall(492:600,2380:2520,m)));

    Ibk(m) = mean([bk1,bk2]);
    % brightfield background processing
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % TODO 4: if Ibk is larger than some threshold, it is not noise! (chnage 
    % this value correpondingly)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     if Ibk(m)>300*(2^8-1)/(2^16-1)&& m>1
%         Ibk(m) = Ibk(m-1);
%     end
    if Ibk(m)>300 
        Ibk(m) = Ibk(m-1);
    end
end
fprintf(['\nfinish loading images\n']);
toc;

%saving cropped input stacks
all_vars_matfile = fullfile(out_dir, ['I_input_stack_cropped','.mat']);
save(all_vars_matfile, 'I_input_stack_cropped');


%% define processing ROI
%Np = [2160, 2560];
%Np = [344, 344];
Np = [902, 902];
%Np = [100, 150];


%% read system parameters
if strcmp(sample_name,'USAF') | strcmp(sample_name,'stained')
    USAF_Parameter();
else
    hela_Parameter();
end

%% load in data: read in the patch from the memory
Imea = double(Iall(nstart(1):nstart(1)+Np(1)-1,nstart(2):nstart(2)+Np(2)-1,:));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TODO 5: Define the LED index
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% in this example, the LEDs are scanned sequentially within illumination NA
% defined in the systemsetup();
% In other cases, make the correponding changes. In the end, all we need is
% the matrix 'lit' that specifies the LED index corresponding image.

ledidx = 1:Nled;
ledidx = reshape(ledidx,numlit,Nimg);
lit = Litidx(ledidx); %index of the sequence of LEDs that were used in experiment corresponding to the 293 LEDs in the 32x32 grid!
lit = reshape(lit,numlit,Nimg);
% reorder LED indices based on illumination NA
[dis_lit2,idx_led] = sort(reshape(illumination_na_used,1,Nled)); %reshape NA of the 293 LEDs to be a vector from 1:293 and sort them so that
                                                                 %they store SEPARATELY the index and values of the NA after sorting them!
                                                                 %i.e., idx_led is the SORTED NA-indicies
                                                                 % while dis_lit2 is
                                                                 % the SORTED NA-values
                                                                 
Nsh_lit = zeros(numlit,Nimg); % a vector (in case of numlit==1) of size Nimg==293
Nsv_lit = zeros(numlit,Nimg);

for m = 1:Nimg
    % corresponding index of spatial freq for the LEDs are lit
    lit0 = lit(:,m); %read index of sequence of LEDs that were used in experiment for each image
                     %in case of 1 LED sequence; lit0 just stores the
                     %index of that LED that was on for the corresponding
                     %image!
                     %Therefore, for that image, which LED was on!!
                     
    % this is a translation from the LED plane represented by lit0 or lit
    % to the pupil plane represented by idx_u and idx_v
    Nsh_lit(:,m) = idx_u(lit0); % for the LED that was on, get dirac peak position in u-coordinate (of the horizontal) relative to the center of image
    Nsv_lit(:,m) = idx_v(lit0); % for the LED that was on, get dirac peak position in v-coordinate (of the vertical) relative to the center of image
end

% reorder the LED indices and intensity measurements according the previous
% dis_lit or rather their indices: idx_led
Ns = []; % 3D array where for each LED, the corresponding index of k_u and index of k_v are to be stored
Ns(:,:,1) = Nsv_lit;
Ns(:,:,2) = Nsh_lit;

Imea_reorder = Imea(:,:,idx_led);
Ibk_reorder = Ibk(idx_led);

% gain = a(idx_led);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pre-processing the data to DENOISING is IMPORTANT
% background subtraction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ithresh_reorder = Imea_reorder;
for m = 1:Nimg
    Itmp = Ithresh_reorder(:,:,m);
    Itmp = Itmp-Ibk_reorder(m);
%     Itmp = awgn(Itmp,0,'measured');
    Itmp(Itmp<0) = 0;
    Ithresh_reorder(:,:,m) = Itmp;
    
end

%load('C:\Users\Muneeb\Desktop\ajmal fpm\laura\multiplexed fpm\Ns_cal289.mat\');
% Ns_reorder = Ns(:,idx_led,:);
% Ithresh_reorder = Ithresh_reorder(:,:,1:89);
Ns_reorder = Ns(:,idx_led,:); % 1st element is empty!
                              % 3rd element has the indicies of the k_u and k_v
                              % now store in 2nd element the NA-indicies
                              % idx_led is the NA-indices, e.g.
                              % for USAF samples, the following:
                              % idx_led(1) = 147,
                              % illumination_na_used(147)=0.0012,
                              % dis_lit2(1)=0.0012
                              % Therefore, for each NA-index, store indicies of the corresponding k_u and k_v
clear Imea;
%% reconstruction algorithm
% select the index of images that will be used in the processing
%Nused = 293;
Nused = 1;
idx_used = 1:Nused;
I = Ithresh_reorder(:,:,idx_used);
Ns2 = Ns_reorder(:,idx_used,:);

%% reconstruction algorithm options: opts
%   tol: maximum change of error allowed in two consecutive iterations
    %   maxIter: maximum iterations 
    %   minIter: minimum iterations
    %   monotone (1, default): if monotone, error has to monotonically dropping
    %   when iters>minIter
%   display: display results (0: no (default) 1: yes)
    %   saveIterResult: save results at each step as images (0: no (default) 1: yes)
    %   mode: display in 'real' space or 'fourier' space.
    %   out_dir: saving directory
%   O0, P0: initial guesses for O and P
    %   OP_alpha: regularization parameter for O
    %   OP_beta: regularization parameter for P
%   scale: LED brightness map
%   H0: known portion of the aberration function, 
        % e.g. sample with a known defocus induce a quardratic aberration
        % function can be defined here
%   poscalibrate: flag for LED position correction using
    % '0': no correction
    % 'sa': simulated annealing method
        % calbratetol: parameter in controlling error tolence in sa
    % 'ga': genetic algorithm
    % caution: takes consierably much longer time to compute a single iteration
%   F, Ft: operators of Fourier transform and inverse
opts.tol = 1;
opts.maxIter = 10; 
opts.minIter = 1;
opts.monotone = 1;
% 'full', display every subroutin,
% 'iter', display only results from outer loop
% 0, no display
opts.display = 'full';%0;%'iter';
opts.out_dir = out_dir;

upsamp = @(x) padarray(x,[(N_obj(1)-Np(1))/2,(N_obj(2)-Np(2))/2]);
opts.O0 = F(sqrt(I(:,:,1)));
opts.O0 = upsamp(opts.O0);
opts.P0 = w_NA;
opts.Ps = w_NA;
opts.iters = 1;
opts.mode = 'fourier';
opts.scale = ones(Nused,1);
opts.OP_alpha = 1;
opts.OP_beta = 1e3 ;
opts.poscalibrate =0; % using correction or not
opts.calbratetol = 1e-1;
opts.F = F;
opts.Ft = Ft;
opts.StepSize = 0.1;
opts.saveIterResult = 1;
opts.idx_led = idx_led(:);

disp(opts);
disp(['nstart: ', num2str(nstart)]);
disp(['Np: ', num2str(Np)]);
f88 = [];
%% algorithm starts
[O,P,dirac_cen,err_pc,c,Ns_cal] = AlterMin(I,[N_obj(1),N_obj(2)],round(Ns2),opts);

%% save results for Onofre
recons_img = O;
recons_img_matfile = fullfile(opts.out_dir, ['recons_img','.mat']);
save(recons_img_matfile, "recons_img");

proc_abs_O = abs(O);
proc_abs_O(abs(O)>25) = 25;
recons_img = proc_abs_O;
recons_img_proc_matfile = fullfile(opts.out_dir, ['recons_img_proc','.mat']);
save(recons_img_proc_matfile, "recons_img");

%% save results

%saving high res (here: after going to real) and dirac peak positions as matfile & figure
hig_res_O_matfile = fullfile(opts.out_dir, ['high_res_O','.mat']);
save(hig_res_O_matfile, 'O', 'P', 'dirac_cen', 'idx_led', 'Np', '-v7.3');

% add tag for Np + maxIter
Np_Iter = ['Np_',num2str(Np(1)),'_',num2str(Np(2)),'_minIter_',num2str(opts.minIter),'_maxIter_',num2str(opts.maxIter)];

%real_O = real(O);
%imwrite(uint16(real_O), strcat(out_dir,'O_',Np_Iter,'_image.png'),'BitDepth',16);
%f1 = figure('visible','off');imshow(real_O);
%title('(O)');
%export_fig(f1,strcat(out_dir,'O_',Np_Iter,'_figure.png'),'-m4');

angle_O = angle(O);
imwrite(uint16(angle_O), strcat(out_dir,'angle_O_',Np_Iter,'_image.png'),'BitDepth',16);
f2 = figure('visible','off');imshow(angle_O,[]);
title('angle (O)');
export_fig(f2,strcat(out_dir,'angle_O_',Np_Iter,'_figure.png'),'-m4');

abs_O = abs(O);
imwrite(uint16(abs_O), strcat(out_dir,'abs_O_',Np_Iter,'_image.png'),'BitDepth',16);
f3 = figure('visible','off');imshow(abs_O,[]);
title('abs (O)');
export_fig(f3,strcat(out_dir,'abs_O_',Np_Iter,'_figure.png'),'-m4');

f4 = figure('visible','off'); 
plot(1:length(err_pc), (err_pc))
title('err/iter');
xlabel('iter');
ylabel('err');
export_fig(f4,strcat(out_dir,'err_',Np_Iter,'.png'),'-m4');

f_err_log = figure('visible','off');
plot(1:length(err_pc), (log(err_pc)))
title('err/iter');
xlabel('iter');
ylabel('log err');
export_fig(f_err_log,strcat(out_dir,'err_log_',Np_Iter,'.png'),'-m4');


%% saving variables to matlab files
%all_vars_matfile = fullfile(out_dir, ['all_vars_',Np_Iter,'.mat']);
%%save(all_vars_matfile);
%save(all_vars_matfile, '-v7.3');
% To save all but some large variables,e.g. I, Iall,etc (uncomment next line)
%%save(all_vars_matfile, '-regexp','^(?!(I|Iall|Imea_reorder|Ithresh_reorder)$).');


%post-processing of the white pixels w/ high intensity
%proc_abs_O = abs(O);
%proc_abs_O(abs(O)>25) = 25;
%proc_abs_O(abs(O)>18) = 18;
%imwrite(uint16(proc_abs_O), strcat(out_dir,'proc_abs_O_',Np_Iter,'_image.png'),'BitDepth',16);

f5 = figure('visible','off');imshow(proc_abs_O,[]);
title(['processed abs (O)']);
export_fig(f5,strcat(out_dir,'proc_abs_O_',Np_Iter,'_figure.png'),'-m4');

fprintf('\nfinished Altermin.m\n');
fprintf('\nGetting/saving low res intensities\n');

hightolowRes(O,P,dirac_cen,idx_led,Np,opts);

fprintf('processing completes\n');

% writing README
readmef = fopen( strcat(out_dir,'/','README.txt'), 'wt' );

fprintf('\nwriting the README file.\n');
fprintf(readmef, '\nThis directory has the measured & estimated images of size %2d x %2d.\n', Np(1), Np(2));
for j=1:20, fprintf(readmef,'-'); end
fprintf(readmef,'\n');
fprintf(readmef,'\nParamters used for these computed images:\n');
fprintf(readmef,'\nNp: %2d x %2d.\n', Np(1), Np(2));
fprintf(readmef,'\nminIter: %2d.\n', opts.minIter);
fprintf(readmef,'\nmaxIter: %2d.\n', opts.maxIter);

fclose(readmef);
