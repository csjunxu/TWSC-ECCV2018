% Bundles submission data for sRGB denoising

% submission_folder Folder where denoised images reside

% Output is written to <submission_folder>/bundled/. Please submit
% the content of this folder.
%
% Author: Tobias Plötz, TU Darmstadt (tobias.ploetz@visinf.tu-darmstadt.de)
%
% This file is part of the implementation as described in the CVPR 2017 paper:
% Tobias Plötz and Stefan Roth, Benchmarking Denoising Algorithms with Real Photographs.
% Please see the file LICENSE.txt for the license governing this code.

function bundle_submission_srgb( submission_folder, para)

if nargin < 2
    out_folder = fullfile(submission_folder, 'bundled/');
else
    out_folder = [submission_folder, 'bundled' num2str(para) '/'];
end
mkdir(out_folder);

israw = false;
eval_version = '1.0';

for i=1:50
    Idenoised = cell(1,20);
    parfor b=1:20
        filename = sprintf('%04d_%02d.mat', i,b);
        s = load(fullfile(submission_folder, filename));
        
        Idenoised{b} = s.Idenoised_crop;
    end
    filename = sprintf('%04d.mat', i);
    save(fullfile(out_folder, filename), 'Idenoised', 'israw', 'eval_version', '-v7.3');
end
end

