function  [im_out, Par] = TWSC_RW(Par)
im_out     =   Par.nim;
% parameters for noisy image
[h,  w, ch] = size(im_out);
Par.h = h;
Par.w = w;
Par.ch = ch;
Par = SearchNeighborIndex( Par );
% original noisy image to patches
NY = Image2Patch( Par.nim, Par );
Par.Sigma = sqrt(mean(Par.nSig.^2));
for ite  =  1 : Par.Outerloop
    % iterative regularization
    % im_out = im_out + Par.delta * (Par.nim - im_out);
    % image to patches
    Y = Image2Patch( im_out, Par );
    % estimate local noise variance, par.lambdals is put here since the MAP and Bayesian rules
    if Par.lambda1 ~= 0
        if ite == 1
            SigmaRow = zeros(size(NY));
        else
            SigmaRow = (NY-Y).^2;
        end
    end
    SigmaCol = Par.lambda2 * sqrt(abs(repmat(Par.Sigma^2, 1, size(Y,2)) - mean((NY - Y).^2))); % Estimated Local Noise Level
    % estimation of noise variance
    if mod(ite-1, Par.Innerloop)==0
        Par.nlsp = max(Par.nlspgap, Par.nlsp - Par.nlspgap);
        % searching non-local patches
        blk_arr = Block_Matching_RW( Y, Par );
    end
    % Weighted Sparse Coding
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
        % update weight for sparse coding
        Wsc = bsxfun( @rdivide, SigmaCol(index).^2, S + eps );
        % update C by soft thresholding
        B = D' * nDCnlY;
        C = sign(B) .* max( abs(B) - Wsc, 0 );
        % update Y
        nDCnlYhat = D * C;
        % add back DC components
        nlYhat = bsxfun(@plus, nDCnlYhat, DC);
        % aggregation
        Y_hat(:, index) = Y_hat(:, index) + bsxfun(@times, nlYhat, W2);
        W_hat(:, index) = W_hat(:, index) + repmat(W2, [Par.ps2ch, 1]);
    end
    % Reconstruction
    im_out = PGs2Image(Y_hat, W_hat, Par);
end
return;

