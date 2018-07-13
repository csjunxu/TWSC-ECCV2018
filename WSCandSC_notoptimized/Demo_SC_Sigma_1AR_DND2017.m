clear;
Original_image_dir = 'C:\Users\csjunxu\Desktop\CVPR2018 Denoising\dnd_2017\images_srgb\';
fpath = fullfile(Original_image_dir, '*.mat');
im_dir  = dir(fpath);
im_num = length(im_dir);
load 'C:\Users\csjunxu\Desktop\CVPR2018 Denoising\dnd_2017\info.mat';

method = 'SC';
% write image directory
write_MAT_dir = ['C:/Users/csjunxu/Desktop/CVPR2018 Denoising/dnd_2017Results/'];
write_sRGB_dir = [write_MAT_dir method];
if ~isdir(write_sRGB_dir)
    mkdir(write_sRGB_dir)
end

% Parameters
Par.ps = 6;        % patch size
Par.step = 3;     % the step of two neighbor patches
Par.win = 20;    % size of window around the patch

Par.Outerloop = 2; % 6
Par.Innerloop = 2;
Par.nlspini = 70;
Par.delta = 0;
Par.nlspgap = 10;

for lambda = [100]
    Par.lambda = lambda;
    % set Parameters
    % record all the results in each iteration
    Par.PSNR = zeros(Par.Outerloop, im_num*20, 'double');
    Par.SSIM = zeros(Par.Outerloop, im_num*20, 'double');
    CCPSNR = [];
    CCSSIM = [];
    for i = 1:im_num
        Par.image = i;
        load(fullfile(Original_image_dir, im_dir(i).name));
        S = regexp(im_dir(i).name, '\.', 'split');
        [h,w,ch] = size(InoisySRGB);
        for j = 1:size(info(1).boundingboxes,1)
            Par.nlsp = Par.nlspini;  % number of non-local patches
            IMinname = [S{1} '_' num2str(j)];
            bb = info(i).boundingboxes(j,:);
            Par.nim = InoisySRGB(bb(1):bb(3), bb(2):bb(4),1:3);
            Par.I = Par.nim;
            % noise estimation
            for c = 1:ch
                % Par.nSig0(c) = NoiseLevel(Par.nim(:, :, c));
                Par.nSig(c) = NoiseEstimation(Par.nim(:, :, c)*255, Par.ps)/255;
            end
            % initial PSNR and SSIM
            fprintf('%s: \n', IMinname);
            % denoising
%             t1=clock;
            [IMout, Par]  =  SC_Sigma_1AR(Par);
%             t2=clock;
%             etime(t2,t1)
%             alltime(Par.image)  = etime(t2, t1);
            % calculate the PSNR
            Par.PSNR(Par.Outerloop, Par.image)  =   csnr( IMout*255, Par.I*255, 0, 0 );
            Par.SSIM(Par.Outerloop, Par.image)      =  cal_ssim( IMout*255, Par.I*255, 0, 0 );
            %% output
            fprintf('The final PSNR = %2.4f, SSIM = %2.4f. \n', Par.PSNR(Par.Outerloop, Par.image), Par.SSIM(Par.Outerloop, Par.image));
            %% output
            IMoutname = sprintf([write_sRGB_dir '/' method '_DND_' num2str(lambda) '_' IMinname '.png']);
            imwrite(IMout, IMoutname);
        end
    end
end

