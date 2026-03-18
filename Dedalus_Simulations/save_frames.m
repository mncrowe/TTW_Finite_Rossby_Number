function save_frames(f,savename,ext,colormap,scale,NaN_col)
% saves 3D field as N frames f(:,:,i) for i = {1,..,N}
% - f: f(x,y,i), with x as horizontal axis (column)
% - savename: string, name of saved images, suffixed by '_00i', excluding file extension
% - ext: string, file extension
% - colormap: colormap, matrix of size [N 3] describing colour scale
% - scale: vector [v1 v2], fix the ends of the colorbar to these values
% - NaN_col: color of NaN values, [R G B] vector

if nargin < 3; ext = 'png'; end
if nargin < 4; colormap = cmap; end
if nargin < 5; scale = 0; end
if nargin < 6; NaN_col = 0; end

if ext == 0; ext = 'png'; end
if colormap == 0; colormap = cmap; end

s = size(f); N = s(3);
d = 1+floor(log(N)/log(10));

for i = 1:N
    fi = squeeze(f(:,:,i));
    save_image(fi,[savename '_' pad_zeros(i,d)],ext,colormap,scale,NaN_col)
end

end

