function img = renderText(someText, fontName, fontSize, sampsPerPt, antiAlias, fractionalMetrics, bold)
%
% img = renderText(someText, fontName, fontSize, sampsPerPt, [antiAlias], [fractionalMetrics], [bold])
%
% Makes some Java calls to render text into a 2d image. Make sure you launch
% matlab with the JVM enabled!
%
% Example:
%  img = renderText('WHEN we look to the individuals of the same variety...','Sans',18,2);
%  figure; imagesc(img); axis equal tight off; colormap gray;
%
% HISTORY:
% 2003.??.?? Bob Dougherty (RFD) wrote it.
% 2005.10.31 RFD: added some comments and added it to exptTools2.

if(~exist('fontName','var') || isempty(fontName))
    fontName = 'Courier';
end
if(~exist('fontSize','var') || isempty(fontSize))
    fontSize = 24;
end
if(~exist('sampsPerPt','var') || isempty(sampsPerPt))
    sampsPerPt = 2;
end
if(~exist('antiAlias','var') || isempty(antiAlias))
    antiAlias= 0;
end
if(~exist('fractionalMetrics','var') || isempty(fractionalMetrics))
    fractionalMetrics = 0;
end
if(~exist('bold','var') || isempty(bold))
    bold = false;
end

import java.awt.*;

if(bold)
    f = Font(fontName, Font.BOLD, fontSize);
else
    f = Font(fontName, Font.PLAIN, fontSize);   
end

idXform = geom.AffineTransform;
frc = font.FontRenderContext(idXform, antiAlias, fractionalMetrics);

glyphVec = f.createGlyphVector(frc, someText);
r = glyphVec.getLogicalBounds;
glyphShape = glyphVec.getOutline;
%r = glyphShape.getBounds;

% I couln't figure out how to avoid this ugly, painful loop. Matlab doesn't
% seem to allow me to use more efficient Java calls (see dead code below).
img = zeros(ceil(r.height*sampsPerPt), ceil(r.width*sampsPerPt));
for(x=linspace(r.x, r.x+r.width, (r.width+1)*sampsPerPt))
    for(y=linspace(r.y, r.y+r.height, (r.height+1)*sampsPerPt))
        img(floor((y-r.y+1)*sampsPerPt),floor((x-r.x+1)*sampsPerPt)) = glyphShape.contains(x,y);
    end
end

return

% % DOESN'T WORK- no way to convince matlab to allow the currentSegment method
% % to return a double array.
% pi = glyphShape.getPathIterator(idXform);
% 
% %tmpCoords = javaArray('java.lang.Double', 6)
% tmpCoords = zeros(6,1);
% ii = 1;
% clear shape;
% while(~pi.isDone)
%     shape(ii).segType = pi.currentSegment(tmpCoords);
%     shape(ii).coords = tmpCoords;
%     pi.next;
%     ii = ii+1;
% end
