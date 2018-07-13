function       [Y, Sigma] = Image2Patch( im_out, im_in, par )
% record the non-local patch set and the index of each patch in
% of seed patches in image
im_out         =  single(im_out);
Y          =  zeros(par.ps2ch, par.maxrc, 'single');
NY          =  zeros(par.ps2ch, par.maxrc, 'single');
k    =  0;
for l = 1:par.ch
    for i = 1:par.ps
        for j = 1:par.ps
            k    =  k+1;
            blk  = im_out(i:end-par.ps+i,j:end-par.ps+j, l);
            Y(k,:) = blk(:)';
            nblk  = im_in(i:end-par.ps+i,j:end-par.ps+j, l);
            NY(k,:) = nblk(:)';
        end
    end
end
Sigma = sqrt(abs(repmat(par.nSig0^2, 1, size(Y,2)) - mean((NY - Y).^2))); %Estimated Local Noise Level
