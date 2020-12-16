  function print_figs(fig_nums,pre,prt_opt);

% prints figures using "print -depsc f#" where # is figure number
%
% SYNTAX:
%         print_figs;                % Prints all open figures
%         print_figs(fig_nums);      % Prints fig_nums vector of figures
%         print_figs(fig_nums,pre);  % Uses pre for name of file prepend (fig1)
%         print_figs(fig_nums,[],prt_opt);  % use different -depsc output
%
% INPUT: [default]
%    fig_nums=  list of figure numbers to print
%               [defaults to all open figures]
%    pre=       name to prepend infront of # in file name ['f']
%    prt_opt=   printer device option ['-dpesc']
%
% OUTPUT: figures are created with f#.eps file names where # are figure numbers

% =======================================
% Jim Leuer 10/24/2005 Leuer@fusion.gat.com
% =======================================
% default input

  if nargin == 0
    fig_nums= sort(get(0,'children'));
    pre= 'f';
    prt_opt= '-depsc';
  elseif nargin == 1
    pre= 'f'
    prt_opt= '-depsc';
  elseif nargin == 2
    prt_opt= '-depsc';
  end
  
  if isempty(pre)
     pre= 'f';
  end

  if isempty(fig_nums)
     fig_nums= sort(get(0,'children'));
  end

  if isempty(prt_opt)
     prt_opt= '-depsc';
  end

% ========================
% start print
  for ii=1:length(fig_nums)
     print(fig_nums(ii),prt_opt,[pre int2str(fig_nums(ii))])
     disp([' % print_figs outputing graphics file: ', pre,...
            int2str(fig_nums(ii)),'.eps'])
  end

  return

% =======================
% testing
%   print_figs(1)
%   print_figs
%   print_figs([],'fig')
