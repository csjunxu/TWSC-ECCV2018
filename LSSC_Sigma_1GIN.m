function  [im_out, Par] = LSSC_Sigma_1GIN(Par)
im_out    =   Par.pim;
% parameters for noisy image
[h,  w, ch]      =  size(im_out);
Par.h = h;
Par.w = w;
Par.ch = ch;
Par = SearchNeighborIndex( Par );
% % original noisy image to patches
% [NX, NY] = Image2PG( Par.nim, Par );
% [PX, PY]  = Image2PG( Par.pim, Par );
for ite  =  1 : Par.Outerloop
    %     % iterative regularization
    im_out = im_out+Par.delta*(Par.nim - im_out);
    % image to patches and estimate local noise variance
    [X, Y] = Image2PG( im_out, Par );
    % estimation of noise variance
    if mod(ite-1, Par.Innerloop)==0
        Par.nlsp = Par.nlsp - Par.nlspgap;
        % searching  non-local patches
        blk_arr = Block_Matching( X, Par );
    end
    % Weighted Sparse Coding
    Y_hat = zeros(Par.ps2ch, Par.lenrc, 'double');
    W_hat = zeros(Par.ps2ch, Par.lenrc, 'double');
    for i = 1:Par.lenrc
        index = blk_arr(:, i);
        nlY = X( : , index );
        y        = Y(:, i);
        DC = mean(nlY, 2);
        nDCnlY = bsxfun(@minus, nlY, DC);
        nDCy = y - DC;
        % update D and S
        [D, ~, ~] = svd( full(nDCnlY), 'econ' );
        % update C by soft thresholding
        B = D' * nDCy;
        C = sign(B) .* max( abs(B) - Par.lambda, 0 );
        % update Y
        nDCyhat = D * C;
        % add back DC components
        nlYhat = bsxfun(@plus, nDCyhat, DC);
        % aggregation
        Y_hat(:, i) = Y_hat(:, i) + nlYhat;
        W_hat(:, i) = W_hat(:, i) + ones(Par.ps2ch, 1);
    end
    % Reconstruction
    im_out = Patch2Image(Y_hat, W_hat, Par);
    % calculate the PSNR and SSIM
    PSNR =   csnr( im_out, Par.I, 0, 0 );
    SSIM      =  cal_ssim( im_out, Par.I, 0, 0 );
    fprintf('Iter %d : PSNR = %2.4f, SSIM = %2.4f\n', ite, PSNR, SSIM);
    Par.PSNR(ite, Par.image) = PSNR;
    Par.SSIM(ite, Par.image) = SSIM;
end
return;

