function  [im_out, par] = LSSC_Sigma_1AG(par)
im_out    =   par.nim;
% parameters for noisy image
[h,  w, ch]      =  size(im_out);
par.h = h;
par.w = w;
par.ch = ch;
par = SearchNeighborIndex( par );
for ite  =  1 : par.Outerloop
    %     % iterative regularization
    im_out = im_out+par.delta*(par.nim - im_out);
    % image to patches and estimate local noise variance
    Y = Image2PatchNew( im_out, par );
    % estimation of noise variance
    dif = mean( mean( mean( (par.nim-im_out).^2 ) ) );
    par.sigma = sqrt( abs( par.nSig^2 - dif ) );
    
    % estimation of noise variance
    if mod(ite-1, par.Innerloop)==0
        par.nlsp = par.nlsp - par.nlspgap;
        % searching  non-local patches
        blk_arr = Block_Matching( Y, par );
    end
    % Weighted Sparse Coding
    Y_hat = zeros(par.ps2ch, par.maxrc, 'single');
    W_hat = zeros(par.ps2ch, par.maxrc, 'single');
    for i = 1:par.lenrc
        index = blk_arr(:, i);
        nlY = Y( : , index );
        DC = mean(nlY, 2);
        nDCnlY = bsxfun(@minus, nlY, DC);
        % update D and S
        [D, ~, ~] = svd( full(nDCnlY), 'econ' );
        % update C by soft thresholding
        B = D' * nDCnlY;
        C = sign(B) .* max( abs(B) - par.lambda*par.sigma^2, 0 );
        % update Y
        nDCnlYhat = D * C;
        % add back DC components
        nlYhat = bsxfun(@plus, nDCnlYhat, DC);
        % aggregation
        Y_hat(:, index) = Y_hat(:, index) + nlYhat;
        W_hat(:, index) = W_hat(:, index) + ones(par.ps2ch, par.nlsp);
    end
    % Reconstruction
    im_out = PGs2Image(Y_hat, W_hat, par);
    % calculate the PSNR and SSIM
    PSNR =   csnr( im_out, par.I, 0, 0 );
    SSIM      =  cal_ssim( im_out, par.I, 0, 0 );
    fprintf('Iter %d : PSNR = %2.4f, SSIM = %2.4f\n', ite, PSNR, SSIM);
    par.PSNR(ite, par.image) = PSNR;
    par.SSIM(ite, par.image) = SSIM;
end
return;

