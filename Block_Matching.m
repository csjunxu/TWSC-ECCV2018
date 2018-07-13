function       [blk_arr] = Block_Matching( X, Par)
% record the indexs of patches similar to the seed patch
blk_arr   =  zeros(Par.nlsp, Par.lenrc, 'double');
for  i  =  1 : Par.lenrc
    seed = X(:, Par.SelfIndex(i));
    neighbor = X( :, Par.NeighborIndex( 1:Par.NumIndex(i), i ) );
    Dist = mean(bsxfun(@minus, neighbor, seed).^2, 1);
    [~,index]   =  sort(Dist);
    blk_arr(:,i)        =  Par.NeighborIndex( index( 1:Par.nlsp ), i );
end