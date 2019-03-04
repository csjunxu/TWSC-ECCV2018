%-------------------------------------------------------------------------------------------------------------
% This is an implementation of the TWSC algorithm for real-world image denoising
% by processing the YUV channels insead of RGB channels, it would improve
% TWSC by 0.1dB on PSNR.
% Author:  Jun Xu, csjunxu@comp.polyu.edu.hk / nankaimathxujun@gmail.com
%          The Hong Kong Polytechnic University
%
% Please refer to the following paper if you find this code helps:
%
% @article{TWSC_ECCV2018,
% 	author = {Jun Xu and Lei Zhang and David Zhang},
% 	title = {A Trilateral Weighted Sparse Coding Scheme for Real-World Image Denoising},
% 	journal = {ECCV},
%   year = {2018}
% }
% Please see the file License.txt for the license governing this code.
%-------------------------------------------------------------------------------------------------------------
clear,clc;
load '/home/csjunxu/Github/data/NTIRE2019/ValidationNoisyBlocksSrgb.mat';

method = 'TWSCrgb2yuv';
dataset = 'SIDD_2018';
% write image directory
write_MAT_dir = ['/home/csjunxu/Github/data/NTIRE2019/'];
% Parameters
Par.ps = 6;        % patch size
Par.step = 3;      % the step of two neighbor patches
Par.win = 20;      % size of window around the patch
Par.Outerloop = 4;
Par.Innerloop = 2;
Par.nlspini = 70;
Par.display = 0;
Par.delta = 0;
Par.nlspgap = 0; %10
Par.lambda1 = 0;
% Par.lambda2 = 3;
for lambda2 = [4]
    results = zeros(size(ValidationNoisyBlocksSrgb));
    Par.lambda2 =lambda2; 
    save([write_MAT_dir 'TWSC' num2str(lambda2) '.mat'], 'results');
    for i = 1:size(ValidationNoisyBlocksSrgb,1)
        for j = 1:size(ValidationNoisyBlocksSrgb,2)
            Par.nlsp = Par.nlspini;  % number of non-local patches
            % iterate over bounding boxes
            Par.nim = im2double(squeeze(ValidationNoisyBlocksSrgb(i,j,:,:,:)));
            Par.I = Par.nim;
            Par.nim = rgb2ycbcr(Par.nim);
            [h,w,ch] = size(Par.nim);
            % noise estimation
            for c = 1:ch
                Par.nSig(c) = NoiseEstimation(Par.nim(:, :, c)*255, Par.ps)/255;
            end
            [IMout, Par]  =  TWSC_RW(Par);
            IMout = ycbcr2rgb(IMout);
            fprintf('%s/%s: \n', num2str(i), num2str(j));
            %% output
            results(i,j,:,:,:) = uint8(IMout*255);
        end
        fprintf('Image %d/%d done\n', i,j);
    end
end