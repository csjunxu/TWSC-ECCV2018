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
TT_Original_image_dir = 'nc/';
TT_fpath = fullfile(TT_Original_image_dir, '*.png');
TT_im_dir  = dir(TT_fpath);
im_num = length(TT_im_dir);

method = 'TWSC';
dataset = 'nc';
write_MAT_dir = [dataset '_Results/'];
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
Par.nlspini = 70;
Par.display = 0;
Par.delta = 0;
Par.nlspgap = 10;
Par.lambda1 = 0;
Par.lambda2 = 1; % set randomly as 1, different for each image
% set Parameters
alltime  = zeros(im_num, 1, 'double');
for i = 1 : im_num
    Par.nlsp = Par.nlspini;  % number of non-local patches
    Par.image = i;
    Par.nim = im2double(imread(fullfile(TT_Original_image_dir, TT_im_dir(i).name) ));
    Par.I = im2double(imread(fullfile(TT_Original_image_dir, TT_im_dir(i).name))); % just for consistency in 'TWSC_Sigma_RW.m' function
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
    [IMout, Par]  =  TWSC_Sigma_RW(Par);
    t2=clock;
    etime(t2,t1)
    alltime(Par.image)  = etime(t2, t1);
    %% output
    imwrite(IMout, [write_sRGB_dir '/' method '_' dataset '_' IMname '.png']);
end