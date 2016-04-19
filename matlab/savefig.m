% Nicholas M. Rathmann, 2014

function savefig(h,FILENAME,FONTSIZE, WIDTHHEIGHT, varargin) 

% THIS IS THE REAL PROTOTYPE, BUT WE USE THE ABOVE WRAPPING
%function savefig(h,FILENAME,FONTSIZE, PAPER, INNER, OUTER) % [left bottom width height]

%-------------------
% INPUTS
%-------------------

if length(varargin) == 1, opts = varargin{1}; else opts = struct(); end
if ~isfield(opts,'crop'),           opts.crop           = 1; end
if ~isfield(opts,'cropmargins'),    opts.cropmargins    = [0 0 0 0]; end
if ~isfield(opts,'legshrink'),      opts.legshrink      = 0; end

htxt = [];
FONTSIZEtxt = [];
if length(h) > 1, 
    htxt = h(2:end); h = h(1);
    FONTSIZEtxt = FONTSIZE(2:end); FONTSIZE = FONTSIZE(1);
end

PAPER = [1 1];
INNER = [1 1 WIDTHHEIGHT(1) WIDTHHEIGHT(2)]; % 1, 1, WIDTH, HEIGHT
OUTER = [1 1 1 1];

FONT = 'Times';
BOLD = 0;

%-------------------
% SET PAPER/FRAME SIZES
%-------------------
% http://se.mathworks.com/help/matlab/ref/axes-properties.html#prop_PlotBoxAspectRatio
set(h, 'units','centimeters');
set(h, 'PaperUnits','centimeters');
set(h, 'PaperPositionMode', 'manual');
pos = get(h,'Position');

A4 = [210/10 297/10];
W=PAPER(1)*A4(1);
H=PAPER(2)*A4(2);

set(h, 'PaperSize', [W H]); 
set(h, 'PaperPosition',[0 0 W H]); 

OP = get(gca,'OuterPosition');% left = OP(); [left bottom width height]
set(gca,'OuterPosition',[OP(1)*OUTER(1) OP(2)*OUTER(2) OP(3)*OUTER(3) OP(4)*OUTER(4)]); % Moves frame around within page boudary

pos = get(gca, 'Position');
set(gca, 'Position', [pos(1)*INNER(1) pos(2)*INNER(2) pos(3)*INNER(3) pos(4)*INNER(4)]); % Scales internal image within its frame

%-------------------
% UPDATE HANDLES
%-------------------

set(0,'ShowhiddenHandles','on');
CHILDEREN = allchild(h);
set(CHILDEREN,'handlevisibility', 'on');
%get(findobj(CHILDEREN,'type','text'),'DisplayName')
TEXTSA = findobj(CHILDEREN,'type','text');
set(TEXTSA,'FontSize',FONTSIZE);
set(findobj(CHILDEREN,'type','axes'),'FontSize',FONTSIZE);
if ~isempty(FONT)
    set(TEXTSA,'FontName',FONT);
    set(findobj(CHILDEREN,'type','axes'),'FontName',FONT);
end
if BOLD
    set(TEXTSA,'FontWeight','bold');
    set(findobj(CHILDEREN,'type','axes'),'FontWeight','bold');
end

for ii = 1:length(htxt), set(htxt(ii),'FontSize',FONTSIZEtxt(ii)); end
if opts.legshrink, legendshrink(0.6); end

%-------------------
% SAVE IT
%-------------------

fprintf('saving figure %s...\n',FILENAME);
print(gcf, '-dpdf','-r600',FILENAME) % USE THIS
if opts.crop
    [status,result]=system(sprintf('pdfcrop --margins "%i %i %i %i" %s %s', ... 
                    opts.cropmargins(1), opts.cropmargins(2), opts.cropmargins(3), opts.cropmargins(4), FILENAME, FILENAME)); 
    if status, disp(result); end
end

