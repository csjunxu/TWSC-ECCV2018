function  [im_out, Par] = SC_Sigma_1AR(Par)
im_out    =   Par.nim;
% parameters for noisy image
[h,  w, ch]      =  size(im_out);
Par.h = h;
Par.w = w;
Par.ch = ch;
Par = SearchNeighborIndex( Par );
Par.nSig = sqrt(mean(Par.nSig.^2));
for ite  =  1 : Par.Outerloop
    %     % iterative regularization
    im_out = im_out+Par.delta*(Par.nim - im_out);
    % image to patches and estimate local noise variance
    Y = Image2PatchNew( im_out, Par );
    % estimation of noise variance
    dif = mean( mean( mean( (Par.nim-im_out).^2 ) ) );
    Par.sigma = sqrt( abs( Par.nSig^2 - dif ) );
    % estimation of noise variance
    if mod(ite-1, Par.Innerloop)==0
        Par.nlsp = Par.nlsp - Par.nlspgap;
        % searching  non-local patches
        blk_arr = Block_Matching( Y, Par );
    end
    % Weighted Sparse Coding
    Y_hat = zeros(Par.ps2ch, Par.maxrc, 'double');
    W_hat = zeros(Par.ps2ch, Par.maxrc, 'double');
    for i = 1:Par.lenrc
        index = blk_arr(:, i);
        nlY = Y( : , index );
        DC = mean(nlY, 2);
        nDCnlY = bsxfun(@minus, nlY, DC);
        % update D and S
        [D, ~, ~] = svd( full(nDCnlY), 'econ' );
        % update C by soft thresholding
        B = D' * nDCnlY;
        C = sign(B) .* max( abs(B) - Par.lambda*Par.sigma^2, 0 );
        % update Y
        nDCnlYhat = D * C;
        % add back DC components
        nlYhat = bsxfun(@plus, nDCnlYhat, DC);
        % aggregation
        Y_hat(:, index) = Y_hat(:, index) + nlYhat;
        W_hat(:, index) = W_hat(:, index) + ones(Par.ps2ch, Par.nlsp);
    end
    % Reconstruction
    im_out = PGs2Image(Y_hat, W_hat, Par);
    % calculate the PSNR and SSIM
    PSNR =   csnr( im_out*255, Par.I*255, 0, 0 );
    SSIM      =  cal_ssim( im_out*255, Par.I*255, 0, 0 );
    fprintf('Iter %d : PSNR = %2.4f, SSIM = %2.4f\n', ite, PSNR, SSIM);
    Par.PSNR(ite, Par.image) = PSNR;
    Par.SSIM(ite, Par.image) = SSIM;
end
return;

