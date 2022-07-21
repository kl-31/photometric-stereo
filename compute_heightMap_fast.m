function Z = compute_heightMap_fast(N, mask)
%COMPUTE_HEIGHTMAP Computes the height map given images, their
%corresponding light directions, and the gray mask
%
%   Z = compute_heightMap(N, mask)
%
%computes the height map "Z", a m-by-n matrix, from the
%surface normals "N", a m-by-n-by-3 matrix whose each element is a
%3-vector and the gray mask "mask".
%
%Author: Xiuming Zhang (GitHub: xiumingzhang), National Univ. of Singapore
%

[im_h, im_w, ~] = size(N);
tic
##% 2D index to 1D object index

#no_pix = size(obj_h, 1);

# vectorized
[h, w] = find(mask);
to_mask = find((h+1)>im_h | (w+1)>im_w);
ind_to_mask = sub2ind(size(mask),h(to_mask),w(to_mask));
mask(ind_to_mask) = 0;
h(to_mask) = [];
w(to_mask) = [];
[obj_h, obj_w] = find(mask);

# matches location in image to index in ids
full2obj = zeros(size(mask));
for id = 1:size(h, 1)
    full2obj(h(id), w(id)) = id;
end

hw_ind = sub2ind(size(mask),h,w);
hw_down = sub2ind(size(mask),h+1,w);
to_remove1 = find(mask(hw_down)==0); # indices of hw_down, not of mask
Mind_remove1 = full2obj(hw_down(to_remove1)-1);
hw_right = sub2ind(size(mask),h,w+1);
to_remove2 = find(mask(hw_right)==0);
Mind_remove2 = full2obj(hw_right(to_remove2)-im_h);

hw_down(to_remove1) = [];
hw_right(to_remove2) = [];

no_pix = length(h);
ids = 1:no_pix;
row_ids1 = (ids-1)*2+1;
row_ids2 = (ids-1)*2+2;

M = sparse(2*no_pix, no_pix);
u = sparse(2*no_pix, 1);

N_x = N(:,:, 1);
N_y = N(:,:, 2);
N_z = N(:,:, 3);
n_x = N_x(hw_ind);
n_y = N_y(hw_ind);
n_z = N_z(hw_ind);

#to_remove1 = find(mask(sub2ind(size(mask),h+1,w))==0);

#to_remove2 = find(mask(sub2ind(size(mask),h,w+1))==0);

u((full2obj(hw_down-1)-1)*2+1) = N_y(hw_down-1);
Mindex1 = sub2ind(size(M), (full2obj(hw_down-1)-1)*2+1, full2obj(hw_down-1));
M(Mindex1) = -N_z(hw_down-1);
Mindex2 = sub2ind(size(M), (full2obj(hw_down-1)-1)*2+1, full2obj(hw_down));
M(Mindex2) = N_z(hw_down-1);

u((full2obj(hw_right-im_h)-1)*2+2) = -N_x(hw_right-im_h);
Mindex3 = sub2ind(size(M), (full2obj(hw_right-im_h)-1)*2+2, full2obj(hw_right-im_h));
M(Mindex3) = -N_z(hw_right-im_h);
Mindex4 = sub2ind(size(M), (full2obj(hw_right-im_h)-1)*2+2, full2obj(hw_right));
M(Mindex4) = N_z(hw_right-im_h);

to_remove = [(Mind_remove1-1)*2+1; (Mind_remove2-1)*2+2];
M(to_remove, :) = [];
u(to_remove, :) = [];

#### Original unvectorized
####------------------------ Assemble M and u
####% 2D index to 1D object index
##full2obj = zeros(size(mask));
##for id = 1:size(obj_h, 1)
##    full2obj(obj_h(id), obj_w(id)) = id;
##end
##
##failed_rows = [];
##for idx = 1:no_pix
##    % Position in 2D image
##    h = obj_h(idx);
##    w = obj_w(idx);
##    % Surface normal
##    n_x = N(h, w, 1);
##    n_y = N(h, w, 2);
##    n_z = N(h, w, 3);
##    % First row - vertical neighbors
##    row_idx = (idx-1)*2+1;
##    % Filter our potentially harmful points
##    if ((h+1)<=im_h) && ((h-1)>=0)
##      if mask(h+1, w) % check if down neighbor is in bound
##            idx_vertN = full2obj(h+1, w);
##            u(row_idx) = n_y;
##            M(row_idx, idx) = -n_z;
##            M(row_idx, idx_vertN) = n_z;
####      elseif mask(h-1, w) % check if up neighbor is in bound
####            idx_vertN = full2obj(h-1, w);
####            u(row_idx) = -n_y;
####            M(row_idx, idx) = -n_z;
####            M(row_idx, idx_vertN) = n_z;
##        else % no vertical neighbors
##            failed_rows = [failed_rows; row_idx];
##        end
##    end
##    % Second row - horizontal neighbors
##    row_idx = (idx-1)*2+2;
##    if ((w+1)<=im_w) && ((w-1)>=1)
##      if mask(h, w+1) % check if right neighbor is in bound
##          idx_horizN = full2obj(h, w+1);
##          u(row_idx) = -n_x;
##          M(row_idx, idx) = -n_z;
##          M(row_idx, idx_horizN) = n_z;
####      elseif mask(h, w-1) % check if left neighbor is in bound
####          idx_horizN = full2obj(h, w-1);
####          u(row_idx) = n_x;
####          M(row_idx, idx) = -n_z;
####          M(row_idx, idx_horizN) = n_z;
##      else % no horizontal neighbors
##          failed_rows = [failed_rows; row_idx];
##      end
##    end
##end
##
## Remove those all-zero rows
##M(failed_rows, :) = [];
##u(failed_rows, :) = [];
##
%------------------------ Solve

z = (M.'*M)\(M.'*u);
% z = qmr(M.'*M, M.'*u);
#z = lsqr(M, u);

% From sparse back to full matrix
z = full(z);

% Outliers due to singularity
outlier_ind = abs(zscore(z))>10;
z_min = min(z(~outlier_ind));
z_max = max(z(~outlier_ind));

%------------------------ Reassemble z back to 2D

Z = double(mask);
for idx = 1:no_pix
    % Position in 2D image
    h = obj_h(idx);
    w = obj_w(idx);
    % Rescale
    Z(h, w) = (z(idx)-z_min)/(z_max-z_min)*255;
end
toc
disp('Height map computed.')
