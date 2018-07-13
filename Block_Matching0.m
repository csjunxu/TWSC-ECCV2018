function       blk_arr = Block_Matching0( X, Par)
% record the indexs of patches similar to the seed patch
blk_arr   =  zeros(Par.nlsp, Par.lenrc, 'double');
for  i  =  1 : Par.lenrc
    seed = X(:, Par.SelfIndex(i));
    neighbor = X( :, Par.NeighborIndex( 1:Par.NumIndex(i), i ) );
    Dist = mean(bsxfun(@minus, neighbor, seed).^2, 1);
    [~,index]   =  sort(Dist);
    %% for synthetic noise
    blk_arr(:,i)        =  Par.NeighborIndex( index( 1:Par.nlsp ), i );
    %% for realistic noisy
    indc        =  Par.NeighborIndex( index( 1:Par.nlsp ), i );
    indc(indc == Par.SelfIndex(i)) = indc(1); % added on 08/01/2017
    indc(1) = Par.SelfIndex(i); % to make sure the first one of indc equals to off
    blk_arr(:, i) = indc;
end
