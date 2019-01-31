% Utility function for denoising all bounding boxes in all sRGB images of
% the DND dataset.
%
% denoiser      Function handle
%               It is called as Idenoised = denoiser(Inoisy, nlf) where Inoisy is the noisy image patch
%               and nlf is a struct containing the  mean noise strength (nlf.sigma)
% data_folder   Folder where the DND dataset resides
% out_folder    Folder where denoised output should be written to
%
% You can parallelize by having a parfor over the bounding boxes.
%
% Author: Tobias Plötz, TU Darmstadt (tobias.ploetz@visinf.tu-darmstadt.de)
%
% This file is part of the implementation as described in the CVPR 2017 paper:
% Tobias Plötz and Stefan Roth, Benchmarking Denoising Algorithms with Real Photographs.
% Please see the file LICENSE.txt for the license governing this code.
%
% modified @ 2019-01-31 by Jun Xu, nankaimathxujun@gmail.com
clear;
data_folder = '../dnd_2017/';
load '../dnd_2017/info.mat';

% iterate over images
for i=1:50
    img = load(fullfile(data_folder, 'images_srgb', sprintf('%04d.mat', i)));
    Inoisy = img.InoisySRGB;
    
    % iterate over bounding boxes
    Idenoised_crop_bbs = cell(1,20);
    for b=1:20
        bb = info(i).boundingboxes(b,:);
        Inoisy_crop = Inoisy(bb(1):bb(3), bb(2):bb(4), :);
        nlf = info(i).nlf;
        
        nlf.sigma = info(i).sigma_srgb(b);
        Idenoised_crop = denoiser(Inoisy_crop, nlf);
        
        
        Idenoised_crop_bbs{b} = single(Idenoised_crop);
    end
    for b=1:20
        Idenoised_crop = Idenoised_crop_bbs{b};
        save(fullfile(output_folder, sprintf('%04d_%02d.mat', i, b)), 'Idenoised_crop');
    end
    fprintf('Image %d/%d done\n', i,50);
end



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
Original_image_dir = '../dnd_2017/images_srgb/';
fpath = fullfile(Original_image_dir, '*.mat');
im_dir  = dir(fpath);
im_num = length(im_dir);
load '../dnd_2017/info.mat';

method = 'TWSC';
dataset = 'dnd_2017';
% write image directory
write_MAT_dir = ['../' dataset '_Results/'];
write_sRGB_dir = [write_MAT_dir method];
if ~isdir(write_sRGB_dir)
    mkdir(write_sRGB_dir)
end

% set parameters
Par.ps   = 6;        % patch size
Par.step = 3;       % the step of two neighbor patches
Par.win  = 20;   % size of window around the patch
Par.Outerloop = 8;
Par.Innerloop = 2;
Par.nlspini = 70;
Par.display = 0;
Par.delta   = 0;
Par.nlspgap = 0; %10
Par.lambda1 = 0;
Par.lambda2 = 3;

alltime  = zeros(im_num, 1, 'double');
for i = 1 :im_num
    Par.image = i;
    load(fullfile(Original_image_dir, im_dir(i).name));
    S = regexp(im_dir(i).name, '\.', 'split');
    [h,w,ch] = size(InoisySRGB);
    % iterate over bounding boxes
    Idenoised_crop_bbs = cell(1,20);
    for j = 1:size(info(1).boundingboxes,1)
        Par.nlsp = Par.nlspini;  % number of non-local patches
        IMinname = [S{1} '_' num2str(j)];
        bb = info(i).boundingboxes(j,:);
        Par.nim = InoisySRGB(bb(1):bb(3), bb(2):bb(4),:);
        Par.I = Par.nim;
        % noise estimation
        for c = 1:ch
            Par.nSig(c) = NoiseEstimation(Par.nim(:, :, c)*255, Par.ps)/255;
        end
        % initial PSNR and SSIM
        fprintf('%s: \n', IMinname);
        % denoising
        t1=clock;
        [IMout, Par]  =  TWSC_Sigma_RW(Par);
        t2=clock;
        etime(t2,t1)
        alltime(Par.image)  = etime(t2, t1);
        %% output
        IMoutname = sprintf([write_sRGB_dir '/' method '_dnd_' IMinname '.png']);
        imwrite(IMout, IMoutname);
        Idenoised_crop_bbs{j} = single(IMout);
    end
    for j = 1:size(info(1).boundingboxes,1)
        Idenoised_crop = Idenoised_crop_bbs{j};
        save(fullfile(write_MAT_dir, sprintf('%04d_%02d.mat', i, b)), 'Idenoised_crop');
    end
    fprintf('Image %d/%d done\n', i,50);
end


