function      im_out = Patch2Image(X, W, par)
% Reconstruction
im_out = zeros(par.h, par.w, par.ch);
im_wei = zeros(par.h, par.w, par.ch);
k = 0;
for l = 1:1:par.ch
    for i = 1:1:par.ps
        for j = 1:1:par.ps
            k = k+1;
            im_out(par.r-1+i, par.c-1+j, l)  =  im_out(par.r-1+i, par.c-1+j, l) + reshape( X(k,:)', [par.lenr par.lenc] );
            im_wei(par.r-1+i, par.c-1+j, l)  =  im_wei(par.r-1+i, par.c-1+j, l) + reshape( W(k,:)', [par.lenr par.lenc] );
        end
    end
end
im_out  =  im_out ./ im_wei;