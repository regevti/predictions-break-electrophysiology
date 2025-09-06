function h = scalebar(ax, lenX, labelX, lenY, labelY, location, varargin)
% SCALEBAR Draw a scale bar in data units on axes AX.
%   h = scalebar(ax, lenX, labelX)                    % horizontal only
%   h = scalebar(ax, lenX, labelX, lenY, labelY)      % + vertical
%   h = scalebar(ax, ..., location)                   % 'sw','se','nw','ne'
%   h = scalebar(ax, ..., 'Color','k','LineWidth',2)  % styling
%
% Inputs:
%   ax       : target axes handle (e.g., gca)
%   lenX     : length of horizontal bar (in data units)
%   labelX   : text under/over horizontal bar (e.g., '50 \mum')
%   lenY     : (optional) length of vertical bar (data units), [] to skip
%   labelY   : (optional) text beside vertical bar (e.g., '200 ms')
%   location : (optional) 'sw' (default), 'se', 'nw', or 'ne'
%
% Returns:
%   h: struct with fields .hx, .tx, .hy, .ty (line/text handles)

if nargin < 6 || isempty(location), location = 'sw'; end
if nargin < 4 || isempty(lenY),    lenY = []; labelY = ''; end

% Parse styling
p = inputParser;
addParameter(p,'Color','k');
addParameter(p,'LineWidth',2);
parse(p, varargin{:});
col = p.Results.Color;
lw  = p.Results.LineWidth;

hold(ax, 'on');
xl = xlim(ax); yl = ylim(ax);
dx = diff(xl); dy = diff(yl);
m  = 0.05;                 % margin from edges (5% of axis span)

% Anchor position (x0,y0) at chosen corner
switch lower(location)
    case 'sw', x0 = xl(1)+m*dx;       y0 = yl(1)+m*dy;           ytxtOff = -0.02*dy; xlblVA = 'top';
    case 'se', x0 = xl(2)-m*dx-lenX;  y0 = yl(1)+m*dy;           ytxtOff = -0.02*dy; xlblVA = 'top';
    case 'nw', x0 = xl(1)+m*dx;       y0 = yl(2)-m*dy;           ytxtOff = +0.02*dy; xlblVA = 'bottom';
    case 'ne', x0 = xl(2)-m*dx-lenX;  y0 = yl(2)-m*dy;           ytxtOff = +0.02*dy; xlblVA = 'bottom';
    otherwise, error('location must be ''sw'',''se'',''nw'',''ne''.');
end

% Horizontal bar + label
h.hx = plot(ax, [x0 x0+lenX], [y0 y0], '-', 'Color', col, 'LineWidth', lw, 'Clipping','off');
h.tx = text(ax, x0 + lenX/2, y0 + ytxtOff, labelX, ...
    'HorizontalAlignment','center','VerticalAlignment',xlblVA,'Color',col,'FontWeight','bold');

% Optional vertical bar (drawn at the right end of horizontal bar)
if ~isempty(lenY) && lenY > 0
    % place upward or downward depending on corner
    if any(strcmpi(location, {'sw','se'}))
        y1 = y0; y2 = y0 + lenY; vtxtVA = 'bottom'; vtxtOff = +0.01*dy;
    else
        y1 = y0; y2 = y0 - lenY; vtxtVA = 'top';    vtxtOff = -0.01*dy;
    end
    h.hy = plot(ax, [x0+lenX x0+lenX], [y1 y2], '-', 'Color', col, 'LineWidth', lw, 'Clipping','off');
    h.ty = text(ax, x0+lenX + 0.01*dx, (y1+y2)/2 + vtxtOff, labelY, ...
        'HorizontalAlignment','left','VerticalAlignment',vtxtVA,'Color',col,'FontWeight','bold');
else
    [h.hy, h.ty] = deal(gobjects(0));
end
end
