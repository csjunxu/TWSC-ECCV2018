function  [im_out, Par] = TWSC_Sigma_AWGN(Par)
im_out    =   Par.nim;
% parameters for noisy image
[h,  w, ch]      =  size(im_out);
Par.h = h;
Par.w = w;
Par.ch = ch;
Par = SearchNeighborIndex( Par );
% original noisy image to patches
NY = Image2Patch( Par.nim, Par );
for ite  =  1 : Par.outerIter
    % iterative regularization
    im_out = im_out + Par.delta * (Par.nim - im_out);
    % image to patches
    Y = Image2Patch( im_out, Par );
    % estimate local noise variance, par.lambdals is put here since the MAP
    % and Bayesian rules
    if Par.lambda1 ~= 0
        SigmaRow = (NY-Y).^2;
    end
    SigmaCol = Par.lambda2 * sqrt(abs(repmat(Par.nSig^2, 1, size(Y,2)) - mean((NY - Y).^2))); %Estimated Local Noise Level
    % estimation of noise variance
    if mod(ite-1, Par.innerIter)==0
        Par.nlsp = max(Par.nlspgap, Par.nlsp - Par.nlspgap);
        % searching  non-local patches
        blk_arr = Block_Matching( Y, Par );
        if ite == 1
            SigmaCol = Par.nSig * ones(size(SigmaCol));
        end
    end
    % Trilateral Weighted Sparse Coding
    Y_hat = zeros(Par.ps2ch, Par.maxrc, 'double');
    W_hat = zeros(Par.ps2ch, Par.maxrc, 'double');
    for i = 1:Par.lenrc
        index = blk_arr(:, i);
        nlY = Y( : , index );
        DC = mean(nlY, 2);
        nDCnlY = bsxfun(@minus, nlY, DC);
        
        % Compute W2
        W2 = 1 ./ (SigmaCol(index) + eps);
        % update D
        [D, S, ~] = svd( full(nDCnlY), 'econ' );
        % update S
        S = sqrt(max( diag(S).^2 - length(index) * SigmaCol(index(1))^2, 0 ));
        if Par.lambda1 == 0
            % update weight for sparse coding
            Wsc = bsxfun( @rdivide, SigmaCol(index).^2, S + eps );
            % update C by soft thresholding
            B = D' * nDCnlY;
            C = sign(B) .* max( abs(B) - Wsc, 0 );
            % update Y
            nDCnlYhat = D * C;
        else
            W1 = exp( - Par.lambda1*mean(SigmaRow(:, index), 2)); % mean or max?
            S = diag(S);
            % min |Z|_1 + |W1(Y-DSC)W2|_F,2  s.t.  C=Z
            C = TWSC_ADMM( nDCnlY, D, S, W1, W2, Par );
            % update Y
            nDCnlYhat = D * S * C;
        end
        % add back DC components
        nlYhat = bsxfun(@plus, nDCnlYhat, DC);
        % aggregation
        Y_hat(:, index) = Y_hat(:, index) + bsxfun(@times, nlYhat, W2);
        W_hat(:, index) = W_hat(:, index) + repmat(W2, [Par.ps2ch, 1]);
        %         Y_hat(:, index) = Y_hat(:, index) + bsxfun(@times, bsxfun(@times, W1, nlYhat), W2);
        %         W_hat(:, index) = W_hat(:, index) + W1 * W2;
    end
    % Reconstruction
    im_out = PGs2Image(Y_hat, W_hat, Par);
    % calculate the PSNR and SSIM
    PSNR = csnr( im_out * 255, Par.I * 255, 0, 0 );
    SSIM = cal_ssim( im_out * 255, Par.I * 255, 0, 0 );
    fprintf('Iter %d : PSNR = %2.4f, SSIM = %2.4f\n', ite, PSNR, SSIM);
    Par.PSNR(ite, Par.image) = PSNR;
    Par.SSIM(ite, Par.image) = SSIM;
end
return;

