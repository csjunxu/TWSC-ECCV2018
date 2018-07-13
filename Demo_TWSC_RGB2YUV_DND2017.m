clear;

Original_image_dir = 'C:\Users\csjunxu\Desktop\CVPR2018 Denoising\dnd_2017\images_srgb\';
fpath = fullfile(Original_image_dir, '*.mat');
im_dir  = dir(fpath);
im_num = length(im_dir);
load 'C:\Users\csjunxu\Desktop\CVPR2018 Denoising\dnd_2017\info.mat';

method = 'TWSCrgb2yuv_y';
% write image directory
write_MAT_dir = ['C:/Users/csjunxu/Desktop/CVPR2018 Denoising/dnd_2017Results/'];
write_sRGB_dir = [write_MAT_dir method];
if ~isdir(write_sRGB_dir)
    mkdir(write_sRGB_dir)
end
dataset = 'DND';
% Parameters
Par.ps = 6;        % patch size
Par.step = 3;       % the step of two neighbor patches
Par.win = 20;   % size of window around the patch
Par.Outerloop = 8;
Par.Innerloop = 2;
Par.maxIter = 10;
Par.maxrho = 100;
Par.nlspini = 70;
Par.model = 1;
Par.display = 0;
Par.delta = 0;
Par.nlspgap = 0; %10

for mu = [1.1]
    Par.mu = mu;
    for rho = [0.5]
        Par.rho = rho;
        for lambda1 = [0]
            Par.lambda1 = lambda1;
            for lambda2 = [3]
                Par.lambda2 = lambda2;
                for nlspini = [70]
                    Par.nlspini = nlspini;
                    % set Parameters
                    % record all the results in each iteration
                    Par.PSNR = zeros(Par.Outerloop, im_num*20, 'double');
                    Par.SSIM = zeros(Par.Outerloop, im_num*20, 'double');
                    CCPSNR = [];
                    CCSSIM = [];
                    for i = 1 :im_num
                        Par.image = i;
                        load(fullfile(Original_image_dir, im_dir(i).name));
                        S = regexp(im_dir(i).name, '\.', 'split');
                        [h,w,ch] = size(InoisySRGB);
                        for j = 1:size(info(1).boundingboxes,1)
                            Par.nlsp = Par.nlspini;  % number of non-local patches
                            IMinname = [S{1} '_' num2str(j)];
                            bb = info(i).boundingboxes(j,:);
                            Par.nim = InoisySRGB(bb(1):bb(3), bb(2):bb(4),:);
                            Par.I = Par.nim;
                            Par.nim = rgb2ycbcr(Par.nim);
                            Par.nimcb = Par.nim(:,:,2);
                            Par.nimcr = Par.nim(:,:,3);
                            Par.nim = Par.nim(:,:,1);
                            [h,w,ch] = size(Par.nim);
                            % noise estimation
                            for c = 1:ch
                                % Par.nSig0(c) = NoiseLevel(Par.nim(:, :, c));
                                Par.nSig(c) = NoiseEstimation(Par.nim(:, :, c)*255, Par.ps)/255;
                            end
                            [IMout, Par]  =  TWSC_Sigma_WAR(Par);
                            IMout(:,:,2) = Par.nimcb;
                            IMout(:,:,3) = Par.nimcr;
                            IMout = ycbcr2rgb(IMout);
                            % initial PSNR and SSIM
                            fprintf('%s: \n', IMinname);
                            
                            % calculate the PSNR
                            Par.PSNR(Par.Outerloop, Par.image)  =   csnr( IMout*255, Par.I*255, 0, 0 );
                            Par.SSIM(Par.Outerloop, Par.image)      =  cal_ssim( IMout*255, Par.I*255, 0, 0 );
                            %% output
                            imwrite(IMout, [write_sRGB_dir '/' method '_' dataset '_l1_' num2str(lambda1) '_l2_' num2str(lambda2) '_N' num2str(nlspini) '_' IMinname '.png']);
                        end
                    end
                end
            end
        end
    end
end

Original_image_dir = 'C:\Users\csjunxu\Desktop\CVPR2018 Denoising\dnd_2017\images_srgb\';
fpath = fullfile(Original_image_dir, '*.mat');
im_dir  = dir(fpath);
im_num = length(im_dir);
load 'C:\Users\csjunxu\Desktop\CVPR2018 Denoising\dnd_2017\info.mat';

method = 'TWSCrgb2yuv_y_u_v';
% write image directory
write_MAT_dir = ['C:/Users/csjunxu/Desktop/CVPR2018 Denoising/dnd_2017Results/'];
write_sRGB_dir = [write_MAT_dir method];
if ~isdir(write_sRGB_dir)
    mkdir(write_sRGB_dir)
end
dataset = 'DND';
% Parameters
Par.ps = 6;        % patch size
Par.step = 3;       % the step of two neighbor patches
Par.win = 20;   % size of window around the patch
Par.Outerloop = 8;
Par.Innerloop = 2;
Par.maxIter = 10;
Par.maxrho = 100;
Par.nlspini = 70;
Par.model = 1;
Par.display = 0;
Par.delta = 0;
Par.nlspgap = 0; %10

for mu = [1.1]
    Par.mu = mu;
    for rho = [0.5]
        Par.rho = rho;
        for lambda1 = [0]
            Par.lambda1 = lambda1;
            for lambda2 = [3]
                Par.lambda2 = lambda2;
                for nlspini = [70]
                    Par.nlspini = nlspini;
                    % set Parameters
                    % record all the results in each iteration
                    Par.PSNR = zeros(Par.Outerloop, im_num*20, 'double');
                    Par.SSIM = zeros(Par.Outerloop, im_num*20, 'double');
                    CCPSNR = [];
                    CCSSIM = [];
                    for i = 1 :im_num
                        Par.image = i;
                        load(fullfile(Original_image_dir, im_dir(i).name));
                        S = regexp(im_dir(i).name, '\.', 'split');
                        [h,w,ch] = size(InoisySRGB);
                        for j = 1:size(info(1).boundingboxes,1)
                            Par.nlsp = Par.nlspini;  % number of non-local patches
                            IMinname = [S{1} '_' num2str(j)];
                            bb = info(i).boundingboxes(j,:);
                            Par.nim = InoisySRGB(bb(1):bb(3), bb(2):bb(4),:);
                            Par.I = Par.nim;
                            Par.nimycbcr = rgb2ycbcr(Par.nim);
                            [h,w,ch] = size(Par.nim);
                            % noise estimation
                            for c = 1:ch
                                Par.nim = Par.nimycbcr(:, :, c);
                                Par.nSig(c) = NoiseEstimation(Par.nim*255, Par.ps)/255;
                                % denoising
                                [IMout(:,:,c), Par]  =  TWSC_Sigma_WAR(Par);
                            end
                            IMout = ycbcr2rgb(IMout);
                            % calculate the PSNR
                            Par.PSNR(Par.Outerloop, Par.image)  =   csnr( IMout*255, Par.I*255, 0, 0 );
                            Par.SSIM(Par.Outerloop, Par.image)      =  cal_ssim( IMout*255, Par.I*255, 0, 0 );
                            %% output
                            imwrite(IMout, [write_sRGB_dir '/' method '_' dataset '_l1_' num2str(lambda1) '_l2_' num2str(lambda2) '_N' num2str(nlspini) '_' IMinname '.png']);
                        end
                    end
                end
            end
        end
    end
end

clear;

Original_image_dir = 'C:\Users\csjunxu\Desktop\CVPR2018 Denoising\dnd_2017\images_srgb\';
fpath = fullfile(Original_image_dir, '*.mat');
im_dir  = dir(fpath);
im_num = length(im_dir);
load 'C:\Users\csjunxu\Desktop\CVPR2018 Denoising\dnd_2017\info.mat';

method = 'TWSCrgb2yuv';
% write image directory
write_MAT_dir = ['C:/Users/csjunxu/Desktop/CVPR2018 Denoising/dnd_2017Results/'];
write_sRGB_dir = [write_MAT_dir method];
if ~isdir(write_sRGB_dir)
    mkdir(write_sRGB_dir)
end
dataset = 'DND';
% Parameters
Par.ps = 6;        % patch size
Par.step = 3;       % the step of two neighbor patches
Par.win = 20;   % size of window around the patch
Par.Outerloop = 8;
Par.Innerloop = 2;
Par.maxIter = 10;
Par.maxrho = 100;
Par.nlspini = 70;
Par.model = 1;
Par.display = 0;
Par.delta = 0;
Par.nlspgap = 0; %10

for mu = [1.1]
    Par.mu = mu;
    for rho = [0.5]
        Par.rho = rho;
        for lambda1 = [0]
            Par.lambda1 = lambda1;
            for lambda2 = [3]
                Par.lambda2 = lambda2;
                for nlspini = [70]
                    Par.nlspini = nlspini;
                    % set Parameters
                    % record all the results in each iteration
                    Par.PSNR = zeros(Par.Outerloop, im_num*20, 'double');
                    Par.SSIM = zeros(Par.Outerloop, im_num*20, 'double');
                    CCPSNR = [];
                    CCSSIM = [];
                    for i = 1 :im_num
                        Par.image = i;
                        load(fullfile(Original_image_dir, im_dir(i).name));
                        S = regexp(im_dir(i).name, '\.', 'split');
                        [h,w,ch] = size(InoisySRGB);
                        for j = 1:size(info(1).boundingboxes,1)
                            Par.nlsp = Par.nlspini;  % number of non-local patches
                            IMinname = [S{1} '_' num2str(j)];
                            bb = info(i).boundingboxes(j,:);
                            Par.nim = InoisySRGB(bb(1):bb(3), bb(2):bb(4),:);
                            Par.I = Par.nim;
                            Par.nim = rgb2ycbcr(Par.nim);
                            [h,w,ch] = size(Par.nim);
                            % noise estimation
                            for c = 1:ch
                                % Par.nSig0(c) = NoiseLevel(Par.nim(:, :, c));
                                Par.nSig(c) = NoiseEstimation(Par.nim(:, :, c)*255, Par.ps)/255;
                            end
                            [IMout, Par]  =  TWSC_Sigma_WAR(Par);
                            IMout = ycbcr2rgb(IMout);
                            % initial PSNR and SSIM
                            fprintf('%s: \n', IMinname);
                            % calculate the PSNR
                            Par.PSNR(Par.Outerloop, Par.image)  =   csnr( IMout*255, Par.I*255, 0, 0 );
                            Par.SSIM(Par.Outerloop, Par.image)      =  cal_ssim( IMout*255, Par.I*255, 0, 0 );
                            %% output
                            imwrite(IMout, [write_sRGB_dir '/' method '_' dataset '_l1_' num2str(lambda1) '_l2_' num2str(lambda2) '_N' num2str(nlspini) '_' IMinname '.png']);
                        end
                    end
                end
            end
        end
    end
end