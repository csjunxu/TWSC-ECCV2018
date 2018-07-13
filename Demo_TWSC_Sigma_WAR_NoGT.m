clear;
GT_Original_image_dir = 'C:\Users\csjunxu\Desktop\CVPR2018 Denoising\d1_Results\Real_NoisyImage\';
GT_fpath = fullfile(GT_Original_image_dir, '*.png');
TT_Original_image_dir = 'C:\Users\csjunxu\Desktop\CVPR2018 Denoising\d1_Results\Real_NoisyImage\';
TT_fpath = fullfile(TT_Original_image_dir, '*.png');
GT_im_dir  = dir(GT_fpath);
TT_im_dir  = dir(TT_fpath);
im_num = length(TT_im_dir);

method = 'TWSC';
dataset = 'd1';
write_MAT_dir = ['C:/Users/csjunxu/Desktop/CVPR2018 Denoising/' dataset '_Results/'];
write_sRGB_dir = [write_MAT_dir method];
if ~isdir(write_sRGB_dir)
    mkdir(write_sRGB_dir)
end

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
Par.display = 1;
Par.delta = 0;
Par.nlspgap = 10;
for mu = [1.1]
    Par.mu = mu;
    for rho = [0.5]
        Par.rho = rho;
        for lambda1 = [0]
            Par.lambda1 = lambda1;
            for lambda2 = [1.4 1.8 1.6 1.2 1]
                Par.lambda2 = lambda2;
                % set Parameters
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
                    fprintf('%s: \n', TT_im_dir(i).name);
                    % denoising
                    t1=clock;
                    [IMout, Par]  =  WLSWSC_Sigma_WAR(Par);
                    t2=clock;
                    etime(t2,t1)
                    alltime(Par.image)  = etime(t2, t1);
                    %% output
                    imwrite(IMout, [write_sRGB_dir '/' method '_' dataset '_l1_' num2str(lambda1) '_l2_' num2str(lambda2) '_' IMname '.png']);
                end
            end
        end
    end
end