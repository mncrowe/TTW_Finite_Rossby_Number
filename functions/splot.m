function splot(f, x, y, options)
% plots f(x,y,n) as a multi-panel plot with n panels
%
% - f: f(x,y,n), array for n sub-figures
% - x: x array for each sub-figure (default: 1:Nx)
% - y: y array for each sub-figure (default: 1:Ny)
% - options:
%   - xlabel: x caption, string
%   - ylabel: y caption, string
%   - title: string array of titles for each figure
%   - type:   0 - surface plot
%             1 - pcolor plot (default)
%             n - (>1) contour plot with n contours
%   - fig_size: [s_x s_y] for width s_x and height s_y, [0 0] gives default size
%   - close_all: true - close all current figures, false - do nothing (default)
%   - panel_array: [n_x n_y] for n_x by n_y array of subfigures
%   - new_fig: true - open new figure (default), false - use current figure window
%   - font_size: size of font for axis labels (default: 12)
%
% ----------------------------------------------------------------------------
% Note: The array f should have first dimension corresponding to x, second
%       dimension to y. If f is entered with x and y flipped f will be 
%       transposed unless length(x) = length(y).
% ----------------------------------------------------------------------------

arguments
    f                       (:,:,:) double
    x                       (:,1)   double  = 1:size(f, 1)
    y                       (:,1)   double  = 1:size(f, 2)
    options.xlabel          (1,1)   string  = ""
    options.ylabel          (1,1)   string  = ""
    options.title           (1,:)   string  = ""
    options.type            (1,:)   double  = 1
    options.fig_size        (1,2)   double  = [0, 0]
    options.close_all       (1,1)   logical = false
    options.panel_array     (1,:)   int64   = 0
    options.new_fig         (1,1)   logical = true
    options.font_size       (1,1)   int64   = 12
    options.show_contours   (1,:)   double  = 0
    options.line_width      (1,1)   double  = 0.5
    options.xlim            (1,2)   double  = [min(min(x)) max(max(x))]
    options.ylim            (1,2)   double  = [min(min(y)) max(max(y))]
end

% Close existing figures and open new figure if required:

if options.close_all; close all; end
if options.new_fig; figure; end

% Get size of f:

[Nx, Ny, Np] = size(f);

% Check size of x and y match size of first two dimensions of f:

if length(x) ~= Nx; error('Length of x must match the first dimension of f'); end
if length(y) ~= Ny; error('Length of y must match the first dimension of f'); end

% Get size of subfigure array (required only for Np > 1):

if isequal(options.panel_array, 0)
    [nx, ny] = sq_fac(Np);
else
    nx = double(options.panel_array(1)); ny = double(options.panel_array(2));
end

% Set figure window size and swap nx and ny if the figure is longer vertically:

if ~isequal(options.fig_size, [0, 0])
    pos = get(gcf, 'position');
    pos(3) = options.fig_size(1);
    pos(4) = options.fig_size(2);
    set(gcf, 'position', pos)
    if pos(4) > pos(3); [nx, ny] = deal(ny, nx); end
end

% Extend options.type to all subfigures:

if length(options.type) < Np; options.type = repmat(options.type(1), [1 Np]); end

% Create figure:

if Np == 1

    % Creat surface plot, pseudocolour plot or contour plot depending on type:

    hold on

    switch options.type
        case 0
            surf(x, y, f')
        case 1
            pcolor(x, y, f')
        otherwise
            contour(x, y, f', options.type, LineWidth = options.line_width)
    end

    if options.show_contours > 0 | numel(options.show_contours) > 1
        contour(x, y, f', options.show_contours, 'k', LineWidth = options.line_width)
    end

    hold off

    % Align plot, set axis limits, add colourbar and interpolated shading:

    view(0, 90)
    shading interp
    colorbar
    xlim(options.xlim)
    ylim(options.ylim)
    colormap(cmap())

    % Add x and y axis labels and title:

    if options.xlabel ~= ""; xlabel(options.xlabel); end
    if options.ylabel ~= ""; ylabel(options.ylabel); end
    if options.title ~= ""; title(options.title); end

    % Set font and line size, add bounding rectangle:

     box off
     set(gca, 'FontSize', options.font_size);
     rectangle('Position', [options.xlim(1) options.ylim(1) options.xlim(2)-options.xlim(1) ...
        options.ylim(2)-options.ylim(1)], 'LineWidth', options.line_width)

else

    hold on;

    % Create Np subplots:

    for n = 1:Np

        subplot(ny, nx, n);

        % Creat surface plot, pseudocolour plot or contour plot depending on type:

        hold on

        switch options.type(n)
            case 0
                surf(x, y, f(:, :, n)')
            case 1
                pcolor(x, y, f(:, :, n)')
            otherwise
                contour(x, y, f(:, :, n)', options.type(n), LineWidth = options.line_width)
        end

        if options.show_contours > 0 | numel(options.show_contours) > 1
            contour(x, y, f(:, :, n)', options.show_contours, 'k', LineWidth = options.line_width)
        end

        hold off

        % Align plot, set axis limits, add colourbar and interpolated shading:

        view(0, 90)
        shading interp
        colorbar
        xlim(options.xlim)
        ylim(options.ylim)
        colormap(cmap())

        % Add x and y axis labels and title:

        if options.xlabel ~= ""; xlabel(options.xlabel); end
        if options.ylabel ~= ""; ylabel(options.ylabel); end
        if length(options.title) == Np; title(options.title(n)); end
        
        % Set font and line size, add bounding rectangle:
        
        box off
        set(gca, 'FontSize', options.font_size);
        rectangle('Position', [options.xlim(1) options.ylim(1) options.xlim(2)-options.xlim(1) ...
            options.ylim(2)-options.ylim(1)], 'LineWidth', options.line_width)


    end

    % Add title to whole figure:

    if length(options.title) == 1; sgtitle(options.title); end

    hold off;

end

end

function [a,b] = sq_fac(n)
% returns the two factors, [a b], closest to sqrt(n), where a*b=n
%
% - n: number to factorise
%
% ----------------------------------------------------------------------------
% Note: cleverer methods should be used when dealing with large numbers
% ----------------------------------------------------------------------------

i=floor(sqrt(n));
while n/i ~= floor(n/i)
    i = i-1; 
end
a = n/i;
b=i;

end

