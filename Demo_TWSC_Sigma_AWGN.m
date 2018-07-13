clear;
Original_image_dir  =    'Path to clean images';
Sdir = regexp(Original_image_dir, '\', 'split');
fpath = fullfile(Original_image_dir, '*.png');
im_dir  = dir(fpath);
im_num = length(im_dir);

method = 'TWSC';
writematpath = '/Users/xujun/Desktop/TWSC/Results_AWGN/';
writefilepath  = [writematpath method '/'];
if ~isdir(writefilepath)
    mkdir(writefilepath);
end
for nSig = [15 25 35 50 75]
    %% Parameters
    Par.model = 1;
    Par.innerIter = 2;
    Par.win = 30;
    Par.lambda1 = 0;
    Par.ps = 8;
    Par.outerIter = 10;
    Par.step = 3;
    Par.nlspini = 90;
    Par.nlspgap = 10;
    if 0 < nSig <= 20
        Par.outerIter = 8;
        Par.delta = .07;
        Par.nlspini = 70;
        Par.lambda2 = .9;
    elseif 20 < nSig <= 30
        Par.delta = .06;
        Par.lambda2 = .76;
    elseif 30 < nSig <= 40
        Par.delta = .07;
        Par.lambda2 = .78;
    elseif 40 < nSig <= 60
        Par.nlspini = 120;
        Par.nlspgap = 15;
        Par.delta = .05;
        Par.lambda2 = .72;
    elseif 60 < nSig <= 80
        Par.ps = 9;
        Par.outerIter = 14;
        Par.step = 4;
        Par.nlspini = 140;
        Par.delta = .05;
        Par.lambda2 = .68; % .66
    else
        disp('Please tune the above parameters by yourself, thanks!');
    end
    % record all the results in each iteration
    Par.PSNR = zeros(Par.outerIter, im_num, 'double');
    Par.SSIM = zeros(Par.outerIter, im_num, 'double');
    T512 = [];
    T256 = [];
    for i = 1:im_num
        Par.nlsp = Par.nlspini;  % number of non-local patches
        Par.image = i;
        Par.nSig = nSig/255;
        Par.I =  im2double( imread(fullfile(Original_image_dir, im_dir(i).name)) );
        S = regexp(im_dir(i).name, '\.', 'split');
        randn('seed',0);
        Par.nim =   Par.I + Par.nSig*randn(size(Par.I));
        fprintf('%s :\n',im_dir(i).name);
        PSNR =   csnr( Par.nim*255, Par.I*255, 0, 0 );
        SSIM      =  cal_ssim( Par.nim*255, Par.I*255, 0, 0 );
        fprintf('The initial value of PSNR = %2.4f, SSIM = %2.4f \n', PSNR,SSIM);
        time0 = clock;
        [im_out, Par]  =  TWSC_Sigma_WA(Par);
        if size(Par.I,1) == 512
            T512 = [T512 etime(clock,time0)];
            fprintf('Total elapsed time = %f s\n', (etime(clock,time0)) );
        elseif size(Par.I,1) ==256
            T256 = [T256 etime(clock,time0)];
            fprintf('Total elapsed time = %f s\n', (etime(clock,time0)) );
        end
        im_out(im_out>1)=1;
        im_out(im_out<0)=0;
        % calculate the PSNR
        Par.PSNR(Par.outerIter, Par.image)  =   csnr( im_out*255, Par.I*255, 0, 0 );
        Par.SSIM(Par.outerIter, Par.image)      =  cal_ssim( im_out*255, Par.I*255, 0, 0 );
        imname = sprintf([writefilepath method '_nSig' num2str(nSig) '_oIte' num2str(Par.outerIter) '_iIte' num2str(Par.innerIter) '_ps' num2str(Par.ps) '_step' num2str(Par.step) '_nlspini' num2str(Par.nlspini) '_nlspgap' num2str(Par.nlspgap) '_delta' num2str(delta) '_l1' num2str(lambda1) '_l2' num2str(lambda2) '_' im_dir(i).name]);
        imwrite(im_out,imname);
        fprintf('%s : PSNR = %2.4f, SSIM = %2.4f \n',im_dir(i).name, Par.PSNR(Par.outerIter, Par.image),Par.SSIM(Par.outerIter, Par.image)     );
    end
    mPSNR=mean(Par.PSNR(end,:),2);
    mSSIM=mean(Par.SSIM(end,:),2);
    mT512 = mean(T512);
    sT512 = std(T512);
    mT256 = mean(T256);
    sT256 = std(T256);
    fprintf('The average PSNR = %2.4f, SSIM = %2.4f. \n', mPSNR,mSSIM);
    name = sprintf([writematpath method '_nSig' num2str(nSig) '_oIte' num2str(Par.outerIter) '_iIte' num2str(Par.innerIter) '_ps' num2str(Par.ps) '_step' num2str(Par.step) '_nlspini' num2str(Par.nlspini) '_nlspgap' num2str(Par.nlspgap) '_delta' num2str(delta) '_l1' num2str(lambda1) '_l2' num2str(lambda2) '.mat']);
    save(name,'nSig','PSNR','SSIM','mPSNR','mSSIM','mT512','sT512','mT256','sT256');
end