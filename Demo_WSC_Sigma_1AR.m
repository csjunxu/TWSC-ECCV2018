clear;
GT_Original_image_dir = 'C:\Users\csjunxu\Desktop\CVPR2018 Denoising\d1_Results\Real_NoisyImage\';
GT_fpath = fullfile(GT_Original_image_dir, '*.png');
TT_Original_image_dir = 'C:\Users\csjunxu\Desktop\CVPR2018 Denoising\d1_Results\Real_NoisyImage\';
TT_fpath = fullfile(TT_Original_image_dir, '*.png');
% GT_Original_image_dir = 'C:\Users\csjunxu\Desktop\CVPR2018 Denoising\cc_Results\Real_ccnoise_denoised_Part\';
% GT_fpath = fullfile(GT_Original_image_dir, '*mean.png');
% TT_Original_image_dir = 'C:\Users\csjunxu\Desktop\CVPR2018 Denoising\cc_Results\Real_ccnoise_denoised_Part\';
% TT_fpath = fullfile(TT_Original_image_dir, '*real.png');
GT_im_dir  = dir(GT_fpath);
TT_im_dir  = dir(TT_fpath);
im_num = length(TT_im_dir);

method = 'WSC';
% write image directory
write_MAT_dir = ['C:/Users/csjunxu/Desktop/CVPR2018 Denoising/d1_Results/'];
write_sRGB_dir = [write_MAT_dir method '/'];
if ~isdir(write_sRGB_dir)
    mkdir(write_sRGB_dir)
end
dataset = 'd1'; % cc

% Parameters
Par.ps = 6;        % patch size
Par.step = 3;       % the step of two neighbor patches
Par.win = 20;   % size of window around the patch

Par.Outerloop = 2;
Par.Innerloop = 2;
Par.nlspini = 70;
Par.delta = 0;
Par.nlspgap = 10;

for lambda = [100]
    Par.lambda = lambda;
    % set Parameters
    % record all the results in each iteration
    Par.PSNR = zeros(Par.Outerloop, im_num, 'double');
    Par.SSIM = zeros(Par.Outerloop, im_num, 'double');
    CCPSNR = [];
    CCSSIM = [];
    for i = [1 4 8] %1 : im_num
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
        Par.nSig = sqrt(sum(Par.nSig.^2)/3);
        % initial PSNR and SSIM
        fprintf('%s: \n', TT_im_dir(i).name);
        CCPSNR = [CCPSNR csnr( Par.nim*255, Par.I*255, 0, 0 )];
        CCSSIM = [CCSSIM cal_ssim( Par.nim*255, Par.I*255, 0, 0 )];
        fprintf('The initial PSNR = %2.4f, SSIM = %2.4f. \n', CCPSNR(end), CCSSIM(end));
        % denoising
        t1=clock;
        [IMout, Par]  =  WSC_Sigma_1AR(Par);
        t2=clock;
        etime(t2,t1)
        alltime(Par.image)  = etime(t2, t1);
        % calculate the PSNR
        Par.PSNR(Par.Outerloop, Par.image)  =   csnr( IMout*255, Par.I*255, 0, 0 );
        Par.SSIM(Par.Outerloop, Par.image)      =  cal_ssim( IMout*255, Par.I*255, 0, 0 );
        %% output
        fprintf('The final PSNR = %2.4f, SSIM = %2.4f. \n', Par.PSNR(Par.Outerloop, Par.image), Par.SSIM(Par.Outerloop, Par.image));
        %% output
        imname = sprintf([write_sRGB_dir 'WSC_' dataset '_'  TT_im_dir(i).name]);
        imwrite(IMout, imname);
    end
    %     mPSNR=mean(Par.PSNR,2);
    %     [~, idx] = max(mPSNR);
    %     PSNR =Par.PSNR(idx,:);
    %     SSIM = Par.SSIM(idx,:);
    %     mSSIM=mean(SSIM,2);
    %     mtime  = mean(alltime);
    %     mCCPSNR = mean(CCPSNR);
    %     mCCSSIM = mean(CCSSIM);
    %     save([write_MAT_dir 'WSC_1AR_delta' num2str(Par.delta) '_ps' num2str(Par.ps) '_step' num2str(Par.step) '_nlspini' num2str(Par.nlspini) '_nlspgap' num2str(Par.nlspgap) '_l' num2str(Par.lambda) '.mat'],'alltime','mtime','PSNR','mPSNR','SSIM','mSSIM','CCPSNR','mCCPSNR','CCSSIM','mCCSSIM');
end