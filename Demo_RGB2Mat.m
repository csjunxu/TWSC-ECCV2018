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
%-------------------------------------------------------------------------------------------------------------
%%% This file will be deleted for clean coding!!!
%-------------------------------------------------------------------------------------------------------------
clear;
Original_image_dir = '../dnd_2017/images_srgb/';
fpath = fullfile(Original_image_dir, '*.mat');
im_dir  = dir(fpath);
im_num = length(im_dir);
load '../dnd_2017/info.mat';

method = 'TWSCrgb2yuv';
dataset = 'dnd_2017';
% write image directory
%write_MAT_dir = ['../' dataset '_Results/'];
write_MAT_dir = ['/home/csjunxu/Paper/DeNoise/Results_' dataset '/'];
write_sRGB_dir = [write_MAT_dir method];
if ~isdir(write_sRGB_dir)
    mkdir(write_sRGB_dir)
end
Par.lambda2 = 3.2;

for i = 1 :im_num
    Par.image = i;
    load(fullfile(Original_image_dir, im_dir(i).name));
    S = regexp(im_dir(i).name, '\.', 'split');
    [h,w,ch] = size(InoisySRGB);
    % iterate over bounding boxes
    Idenoised_crop_bbs = cell(1,20);
    for j = 1:size(info(1).boundingboxes,1)
        IMinname = [S{1} '_' num2str(j)];
        bb = info(i).boundingboxes(j,:);
        Par.nim = InoisySRGB(bb(1):bb(3), bb(2):bb(4),:);
        IMout = im2double(imread([write_sRGB_dir '/' method '_' dataset '_' num2str(Par.lambda2) '_' IMinname '.png']));
        Idenoised_crop_bbs{j} = single(IMout);
    end
    for j = 1:size(info(1).boundingboxes,1)
        Idenoised_crop = Idenoised_crop_bbs{j};
        save(fullfile(write_MAT_dir, sprintf('%04d_%02d.mat', i, j)), 'Idenoised_crop');
    end
    fprintf('Image %d/%d done\n', i,50);
end
% generate submission files
bundle_submission_srgb( write_MAT_dir, Par.lambda2 );
