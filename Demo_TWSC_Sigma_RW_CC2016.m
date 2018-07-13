%-------------------------------------------------------------------------------------------------------------
% This is an implementation of the TWSC algorithm for real-world image denoising
% 
% Author:  Jun Xu, csjunxu@comp.polyu.edu.hk / nankaimathxujun@gmail.com
%          The Hong Kong Polytechnic University
% 
% Please refer to the following paper if you find this code helps:
%
% @article{TWSC_ECCV2018,
% 	author = {Jun Xu and Lei Zhang and David Zhang},
% 	title = {A Trilateral Weighted Sparse Coding Scheme for Real-World Image Denoising},
% 	journal = {ECCV},
% 	year = {2018}
% }
%
% Please see the file License.txt for the license governing this code.
%-------------------------------------------------------------------------------------------------------------
clear;
GT_Original_image_dir = 'Real_ccnoise_denoised_Part/';
GT_fpath = fullfile(GT_Original_image_dir, '*mean.png');
TT_Original_image_dir = 'Real_ccnoise_denoised_Part';
TT_fpath = fullfile(TT_Original_image_dir, '*real.png');
% GT_Original_image_dir = 'C:\Users\csjunxu\Desktop\CVPR2018 Denoising\PolyU\';
% GT_fpath = fullfile(GT_Original_image_dir, '*mean.JPG');
% TT_Original_image_dir = 'C:\Users\csjunxu\Desktop\CVPR2018 Denoising\PolyU\';
% TT_fpath = fullfile(TT_Original_image_dir, '*real.JPG');
GT_im_dir  = dir(GT_fpath);
TT_im_dir  = dir(TT_fpath);
im_num = length(TT_im_dir);

method = 'TWSC';
dataset = 'cc';
write_MAT_dir = ['C:/Users/csjunxu/Desktop/CVPR2018 Denoising/' dataset '_Results/'];
write_sRGB_dir = [write_MAT_dir method];
if ~isdir(write_sRGB_dir)
    mkdir(write_sRGB_dir)
end

% Parameters
Par.ps = 6;        % patch size
Par.step = 3;       % the step of two neighbor patches
Par.win = 20;   % size of window around the patch
Par.Outerloop = 7;
Par.Innerloop = 2;
Par.maxIter = 10;
Par.maxrho = 100;
Par.nlspini = 70;
Par.model = 1;
Par.display = 1;
Par.delta = 0;
Par.nlspgap = 0;
for mu = [1.1]
    Par.mu = mu;
    for rho = [0.5]
        Par.rho = rho;
        for lambda1 = [0]
            Par.lambda1 = lambda1;
            for lambda2 = [4.9]
                Par.lambda2 = lambda2;
                for nlspini = [40:10:100]
                    Par.nlspini = nlspini;
                    % set Parameters
                    % record all the results in each iteration
                    Par.PSNR = zeros(Par.Outerloop, im_num, 'double');
                    Par.SSIM = zeros(Par.Outerloop, im_num, 'double');
                    for i = 1 : im_num
                        Par.nlsp = Par.nlspini;  % number of non-local patches
                        Par.image = i;
                        Par.nim = im2double(imread(fullfile(TT_Original_image_dir, TT_im_dir(i).name) ));
                        Par.I = im2double(imread(fullfile(GT_Original_image_dir, GT_im_dir(i).name)));
                        S = regexp(TT_im_dir(i).name, '\.', 'split');
                        IMname = S{1};
                        [h,w,ch] = size(Par.nim);
                        % noise estimation
                        for c = 1:ch
                            Par.nSig(c) = NoiseEstimation(Par.nim(:, :, c)*255, Par.ps)/255;
                        end
                        % initial PSNR and SSIM
                        fprintf('%s: \n', TT_im_dir(i).name);
                        fprintf('The initial PSNR = %2.4f, SSIM = %2.4f. \n', csnr( Par.nim*255, Par.I*255, 0, 0 ), cal_ssim( Par.nim*255, Par.I*255, 0, 0 ));
                        % denoising
                        [IMout, Par]  =  TWSC_Sigma_WAR(Par);
                        % calculate the PSNR
                        Par.PSNR(Par.Outerloop, Par.image)  =   csnr( IMout*255, Par.I*255, 0, 0 );
                        Par.SSIM(Par.Outerloop, Par.image)      =  cal_ssim( IMout*255, Par.I*255, 0, 0 );
                        %% output
                        imwrite(IMout, [write_sRGB_dir '/' method '_' dataset '_l1_' num2str(lambda1) '_l2_' num2str(lambda2) '_N' num2str(nlspini) '_' IMname '.png']);
                        fprintf('%s : PSNR = %2.4f, SSIM = %2.4f \n',TT_im_dir(i).name, Par.PSNR(Par.Outerloop, Par.image),Par.SSIM(Par.Outerloop, Par.image)     );
                    end
                    mPSNR=mean(Par.PSNR,2);
                    [~, idx] = max(mPSNR);
                    PSNR =Par.PSNR(idx,:);
                    SSIM = Par.SSIM(idx,:);
                    mSSIM=mean(SSIM,2);
                    matname = sprintf([write_MAT_dir method '_' dataset '_l1_' num2str(lambda1) '_l2_' num2str(lambda2) '_N' num2str(nlspini) '.mat']);
                    save(matname,'PSNR','SSIM','mPSNR','mSSIM');
                end
            end
        end
    end
end