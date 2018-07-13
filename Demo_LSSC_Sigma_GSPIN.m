clear;

Original_image_dir  =    'C:\Users\csjunxu\Desktop\TWSCGIN\cleanimages\';
Sdir = regexp(Original_image_dir, '\', 'split');
fpath = fullfile(Original_image_dir, '*.png');
im_dir  = dir(fpath);
im_num = length(im_dir);

method = 'LSSC';
write_MAT_dir = ['C:/Users/csjunxu/Desktop/TWSCGIN/'];
write_sRGB_dir = [write_MAT_dir method];
if ~isdir(write_sRGB_dir)
    mkdir(write_sRGB_dir)
end

for nSig                      =  [10 20 30]         % The standard variance of the additive Gaussian noise;
    for sp                        =  [.1 .3 .5]        % salt and pepper
        Type = 0;
        Par.ps = 8;
        Par.step = 4;
        Par.Outerloop = 10;
        Par.nlspini = 60;
        Par.win = 30;
        Par.model = 1;
        Par.display = 1;
        Par.Innerloop = 2;
        Par.maxIter = 100;
        Par.maxrho = 100;
        Par.nlspgap = 5;
        Par.delta = 0;
        for lambda = [1.1:.1:1.5]
            Par.lambda = lambda;
            % record all the results in each iteration
            Par.PSNR = zeros(Par.Outerloop, im_num, 'single');
            Par.SSIM = zeros(Par.Outerloop, im_num, 'single');
            T512 = [];
            T256 = [];
            for i = 1:im_num
                Par.nlsp = Par.nlspini;  % number of non-local patches
                Par.image = i;
                Par.nSig = nSig;
                Par.I =  single( imread(fullfile(Original_image_dir, im_dir(i).name)) );
                S = regexp(im_dir(i).name, '\.', 'split');
                %                                 %% add Gaussian noise
                %                                 randn('seed',0);
                %                                 Par.nim =   Par.I + Par.nSig*randn(size(Par.I));
                %                                 %% add "salt and pepper" noise
                %                                 rand('seed', 0)
                %                                 nI = imnoise(nI, 'salt & pepper', sp); %"salt and pepper" noise
                %                                 nI = nI*255;
                %                                 %% add "salt and pepper" noise 0 or RVIN noise 1
                %                                 randn('seed',0);
                %                                 [Par.nim,Narr]          =   impulsenoise(Par.nim,sp,Type);
                if Type == 0
                    imname = sprintf([write_MAT_dir 'noisyimages/G' num2str(nSig) '_SPIN' num2str(sp) '_' im_dir(i).name]);
                    %                                     imwrite(Par.nim/255,imname);
                    Par.nim = double( imread(imname));
                elseif Type == 1
                    imname = sprintf([write_MAT_dir 'noisyimages/G' num2str(nSig) '_RVIN' num2str(sp) '_' im_dir(i).name]);
                    %                                     imwrite(Par.nim/255,imname);
                    Par.nim = double( imread(imname));
                else
                    break;
                end
                [Par.pim,ind]           =   adpmedft(Par.nim,19);
                %                                 ind=(Par.pim~=Par.nim)&((Par.nim==255)|(Par.nim==0));
                %                                 Par.pim(~ind)=Par.nim(~ind);
                %
                fprintf('%s :\n',im_dir(i).name);
                PSNR =   csnr( Par.nim, Par.I, 0, 0 );
                SSIM      =  cal_ssim( Par.nim, Par.I, 0, 0 );
                fprintf('The initial value of PSNR = %2.4f, SSIM = %2.4f \n', PSNR,SSIM);
                
                PSNR =   csnr( Par.pim, Par.I, 0, 0 );
                SSIM      =  cal_ssim( Par.pim, Par.I, 0, 0 );
                fprintf('The initial value of PSNR = %2.4f, SSIM = %2.4f \n', PSNR,SSIM);
                Par.nSig = NoiseEstimation(Par.pim, Par.ps);
                %
                time0 = clock;
                [im_out, Par]  =  LSSC_Sigma_1GIN(Par);
                if size(Par.I,1) == 512
                    T512 = [T512 etime(clock,time0)];
                    fprintf('Total elapsed time = %f s\n', (etime(clock,time0)) );
                elseif size(Par.I,1) ==256
                    T256 = [T256 etime(clock,time0)];
                    fprintf('Total elapsed time = %f s\n', (etime(clock,time0)) );
                end
                im_out(im_out>255)=255;
                im_out(im_out<0)=0;
                % calculate the PSNR
                Par.PSNR(Par.Outerloop, Par.image)  =   csnr( im_out, Par.I, 0, 0 );
                Par.SSIM(Par.Outerloop, Par.image)      =  cal_ssim( im_out, Par.I, 0, 0 );
                %                 imname = sprintf('C:/Users/csjunxu/Desktop/NIPS2017/W3Results/SC/LSSC_GSPIN_nSig%d_%1.1f_%s', nSig, sp, im_dir(i).name);
                %                 imwrite(im_out/255,imname);
                fprintf('%s : PSNR = %2.4f, SSIM = %2.4f \n',im_dir(i).name, Par.PSNR(Par.Outerloop, Par.image),Par.SSIM(Par.Outerloop, Par.image)     );
            end
            mPSNR=mean(Par.PSNR,2);
            [~, idx] = max(mPSNR);
            PSNR =Par.PSNR(idx,:);
            SSIM = Par.SSIM(idx,:);
            mSSIM=mean(SSIM,2);
            mT512 = mean(T512);
            sT512 = std(T512);
            mT256 = mean(T256);
            sT256 = std(T256);
            fprintf('The best PSNR result is at %d iteration. \n',idx);
            fprintf('The average PSNR = %2.4f, SSIM = %2.4f. \n', mPSNR(idx),mSSIM);
            name = sprintf([write_MAT_dir '/' method '_GSPIN_p_' Sdir{end-1} '_nSig' num2str(nSig) '_sp' num2str(sp) '_ps' num2str(Par.ps) '_step' num2str(Par.step) '_nlspini' num2str(Par.nlspini) '_nlspgap' num2str(Par.nlspgap) '_delta' num2str(Par.delta) '_l' num2str(Par.lambda) '.mat']);
            save(name,'nSig','sp','PSNR','SSIM','mPSNR','mSSIM','mT512','sT512','mT256','sT256');
        end
    end
end
