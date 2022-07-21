function N = compute_surfNorm_fast(I, L, M)
%COMPUTE_SURFNORM Computes the surface normals given images, their
%corresponding light directions, and the gray mask
%
%   N = compute_surfNorm(I, L)
%
%computes the surface normals "N", a m-by-n-by-3 matrix, from the
%images "I", a p-by-1 cell array whose each element is an image, their
%corresponding lightings "L", a p-by-3 matrix whose each row is the
%lighting direction for the corresponding image in "I", and the gray mask
%"M".
%
%Author: Xiuming Zhang (GitHub: xiumingzhang), National Univ. of Singapore
%
assert(size(I, 1)==size(L, 1), 'I and L mismatched!');
p = size(I, 1);
tic
%------------------------ Collect intensity at a point for all images, T
% Get image dimensions
im = I{1};
[im_h, im_w, ~] = size(im);

##% Initialize T, a im_h-by-im_w-by-p matrix, whose (h, w, :) holds the
##% intensities at (h, w) for all p different lightings
##T = zeros(im_h, im_w, p);
##
##% For each image
##for idx = 1:p
##    im = I{idx};
##    % Loop thru each pixel
##    for h = 1:im_h
##        for w = 1:im_w
##            % If in the mask
##            if M(h, w)
##                r = im(h, w, 1);
##                g = im(h, w, 2);
##                b = im(h, w, 3);
##                inten = norm(double([r g b]));
##                T(h, w, idx) = inten;
##            end
##        end
##    end
##end
##
##
##% Initialize N, a im_h-by-im_w-by-3 matrix, whose (h, w, :) holds the
##% surface norm at (h, w)
##N = zeros(im_h, im_w, 3);
##
##% Loop thru each location
##for h = 1:im_h
##    for w = 1:im_w
##        % If in the mask
##        if M(h, w)
##            % Intensities
##            i = reshape(T(h, w, :), [p, 1]);
##            % Solve surface normals
##            n = (L.'*L)\(L.'*i);
##            if norm(n) ~= 0
##                % Normalize n
##                n = n/norm(n);
##            else
##                n = [0; 0; 0];
##            end
##            % Save
##            N(h, w, :) = n;
##        end
##    end
##end

to_keep = find(M~=0);
to_remove = find(M==0);
for idx = 1:p
  im = double(I{idx});
  T(:,:,idx) = sqrt(im(:,:,1).^2 + im(:,:,2).^2 + im(:,:,3).^2);
end
i = reshape(T, [im_h*im_w,p]);
i(to_remove,:)=[];
i = i.';
n = (L.'*L)\(L.'*i);
n = n.';
n = reshape(n, [length(to_keep), 3]);
N = zeros(im_h*im_w,3);
N(to_keep,:) = n;
norm_N = (sqrt(N(:,1).^2 + N(:,2).^2 + N(:,3).^2));
N = bsxfun(@rdivide,N ,norm_N);
N(find(norm_N==0),:) = zeros(length(find(norm_N==0)),3);
N = reshape(N,[im_h,im_w,3]);
norm_N = reshape(norm_N,[im_h,im_w]);

toc
disp('Surface normal computed.')

