function       [X, Y] = Image2PG( im, par )
% record the non-local patch set and the index of each patch in
% of seed patches in image
im         =  single(im);
X          =  zeros(par.ps2ch, par.maxrc, 'double');
Y          =  zeros(par.ps2ch, par.lenrc, 'double');
k    =  0;
for l = 1:par.ch
    for i = 1:par.ps
        for j = 1:par.ps
            k    =  k+1;
            blk  = im(i:end-par.ps+i,j:end-par.ps+j, l);
            X(k,:) = blk(:)';
            blk     =  im(par.r-1+i, par.c-1+j, l);
            Y(k,:) =  blk(:)';
        end
    end
end