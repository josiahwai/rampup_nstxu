function gs_resize_fonts(h)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  USAGE:   gs_resize_fonts(h)
%           Use in figure's ResizeFcn by:
%           set(gcf,'UserData',h)
%           set(gcf,'ResizeFcn','gs_resize_fonts(get(gcbo,''UserData''))')
%
%  PURPOSE: Resize fonts when width of figure window changes
%
%  INPUTS: h, structure with fields:
%          h.figure = handle to the figure window
%          h.h1 = handles to text of font size 1
%          h.h2 = handles to text of font size 2
%          h.hs = handles to subplots
%          h.f1 = (font size 1)/(window width in pixels)
%          h.f2 = (font size 2)/(window width in pixels)
%          h.fs = (font sizes for subplots)/(window width in pixels)
%
%  OUTPUTS: fonts are resized in the figure
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
%
%  WRITTEN BY:  Anders Welander ON 2015-10-08
%
%  MODIFICATION HISTORY:				
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isfield(h,'figure')
  try
    p = get(h.figure,'Position');
  catch
    return
  end
else
  return
end

if isfield(h,'h1') & isfield(h,'f1')
  for i = 1:numel(h.h1)
    try
      set(h.h1(i),'FontSize',h.f1*p(3))
    end
  end
end

if isfield(h,'h2') & isfield(h,'f2')
  for i = 1:numel(h.h2)
    try
      set(h.h2(i),'FontSize',h.f2*p(3))
    end
  end
end

if isfield(h,'hs') & isfield(h,'fs')
  for i = 1:numel(h.hs)
    j = min(i,numel(h.fs));
    try
      set(h.hs(i),'FontSize',h.fs(j)*p(3))
    end
  end
end

