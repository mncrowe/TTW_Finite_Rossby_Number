function map = cmap(options)
% outputs a custom colour map for use with colormap(map)
%
% - R: plot limits, integer or array;
%    vector: 2 element vector of [Min Max] (default: values of current figure)
%    field: F, calculate limits as min and max of F array
% - centre: value correspnding to the colorbar centre (default: (Min+Max)/2)
% - f: function f:[0,1]->[0,1], nonlinear rescaling around centre (default: f = @(x) x)
% - N: size of colorbar (default: 256)
% - rescale: true - rescale ends to max/min colour bar values, false - no scaling (default: true)
% - C: colours, integer or array;
%    integer: 0 - blue/white/red (default)
%             1 - greyscale
%    array: list of N colours (rgb), N x 3, values in [0,1], can use Matlab presets (e.g. jet(8), parula(3))
% - debug: true - outputs limits and centre value, false - no output (default:false)

arguments
    options.R (:,:) double        = get(gca, 'CLim')
    options.centre (1,:) double   = []
    options.f function_handle     = @(x) x
    options.N (1,1) int64         = 256
    options.rescale (1,1) logical = true
    options.C (:,:) double        = 0
    options.debug (1,1) logical   = false
end

% set ends of colour bar:
if numel(options.R) == 2
    Min = options.R(1); Max = options.R(2);
else
    Min = min(options.R, [], "all"); Max = max(options.R, [], "all");
end

% if a centre is not given, calculate as average of Min and Max:
if numel(options.centre) == 0
    centre = (Min + Max) / 2;
else
    centre = options.centre;
end

% Define colours when using presets:
C = options.C;
if C == 0; C = [0 0 0.5; 0 0.5 1; 1 1 1; 1 0 0; 0.5 0 0];  end
if C == 1; C = [1 1 1; 0.5 0.5 0.5; 0 0 0]; end

% Create normalised array of colour map points:
Tick_lims = (linspace(Min, Max, options.N) - centre) / max(Max - centre, centre - Min);

% Rescale map points to stretch colour bar over full colour range:
if options.rescale
    if Tick_lims(1) > -1; Tick_lims(Tick_lims < 0) = -Tick_lims(Tick_lims < 0) / Tick_lims(1); end
    if Tick_lims(end) < 1; Tick_lims(Tick_lims > 0) = Tick_lims(Tick_lims > 0) / Tick_lims(end); end
end

% Determine colours at each point using interpolation of specified colours:
map = interp1(linspace(-1, 1, length(C)), C, sign(Tick_lims) .* options.f(abs(Tick_lims)));

% Print debug info to screen:
if options.debug
    disp(['Limits = [' num2str([Min Max]) ']'])
    disp(['Centre = ' num2str(centre)])
end

end

