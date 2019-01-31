function       [nDCnlX,blk_arr,DC,par] = CalNonLocal( im, par)
% record the non-local patch set and the index of each patch in
% of seed patches in image
X   = zeros(par.ps2ch, par.maxrc, 'double');
k   = 0;
for l = 1:par.ch
    for i = 1:par.ps
        for j = 1:par.ps
            k = k+1;
            blk = im(i:end-par.ps+i,j:end-par.ps+j, l);
            X(k,:) = blk(:)';
        end
    end
end
% index of each patch in image
Index    =   (1:par.maxrc);
Index    =   reshape(Index,par.maxr,par.maxc);
% record the indexs of patches similar to the seed patch
blk_arr  =  zeros(par.nlsp, par.lenrc ,'double');
% Patch Group Means
DC = zeros(par.ps2ch,par.lenrc,'double');
% non-local patch groups
nDCnlX = zeros(par.ps2ch,par.lenrc*par.nlsp,'double');
% for  i  =  1 :par.lenr
for  j  =  1 : par.lenc
    for  i  =  1 :par.lenr
        row = par.r(i);
        col = par.c(j);
        off = (col-1)*par.maxr + row;
        off1 = (j-1)*par.lenr + i;
        % the range indexes of the window for searching the similar patches
        rmin = max( row - par.win, 1 );
        rmax = min( row + par.win, par.maxr );
        cmin = max( col - par.win, 1 );
        cmax = min( col + par.win, par.maxc );
        idx     =   Index(rmin:rmax, cmin:cmax);
        idx     =   idx(:);
        neighbor = X(:,idx); % the patches around the seed in X
        seed  = X(:,off);
        dis = sum(bsxfun(@minus,neighbor, seed).^2,1);
        [~,ind] = sort(dis);
        indc = idx( ind( 1:par.nlsp ) ); 
        indc(indc==off) = indc(1); % added on 08/01/2017
        indc(1) = off; % to make sure the first one of indc equals to off
        blk_arr(:,off1) = indc;
        temp = X( : , indc );
        DC(:,off1) = mean(temp,2);
        nDCnlX(:,(off1-1)*par.nlsp+1:off1*par.nlsp) = bsxfun(@minus,temp,DC(:,off1));
    end
end