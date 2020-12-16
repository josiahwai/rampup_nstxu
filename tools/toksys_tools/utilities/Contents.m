% Utility functions
% 
% Input/Output functions:
%  hdcopy            - get a hard copy of plot(s) on screen
%  hdcopysys         - get a hard copy of SIMULINK block diagram
%  in_script         - reads variable from call/base area - use for script/fun
%  ot_script         - writes variables to call area      - use for script/fun
%  print_figs        - print all figures (or list of figs) to .eps files
%  read_namelist     - read FORTRAN namelist file
%  read_point        - reads a file (.ptn) which contains list of point names
%  set_save_v6       - makes "save" produce Version 6 compatable .mat files
%  space             - print a blank line to screen
%  wait              - pause with message to terminal
%  write_ascii_table - Write out table to ascii file for easy import into MSWord
%  write_namelist    - Script to write namelist see also read_namelist
%  netcdf_R14        - reads (some) netcdf written files (only Matlab_R14)
%  dlmcell	     - writes a cell array to a file
%
% File/variable management
%  clash             - prints strings that collide with M, Mex, or MAT file name
%  diff_mat          - Compares two .mat files for differences
%  disp_list         - displays a row number in front of input string array
%  dispp             - similar to disp but print to multiple files,and/or screen
%  exist_f_pwd       - does file exist in present working directory?
%  find_conflicts    - find potential conflicts variable name conflicts
%  keep		     - keeps specified variables and clears the rest
%  mk_var            - makes a variable with value only if it does not exist
%  parse_filename    - Returns [dir,name,ext] of input string filename
%  save_mat          - save list of variables into a .mat file 
%  struct_names      - list of all variables & structure names in a structure
%  struct_to_ws      - puts 1st fields of structure into variables in workspace
%  struct_to_ws_all  - Puts ALL Variables of structure into workspace
%  whoss             - extension of whos that shows variable values
%
% String manipulation:
%  fix_undscr        - exchanges "_" with "\_" to eliminate subscripting in text
%  insert_string     - insert a substring into another string at given indices
%  intvec2str        - convert vector of integers to string
%  isblank           - determine if string(s) is (are) blank
%  remove_space      - remove space from different locations in string (deblank)
%  strlen            - find length of string
%  strsfind          - like strfind except operates on multiple rows of patterns
%  strcompare        - single string from 1st identical charcters of string array
%  strsmatch         - find matching strings in two arrays of strings
%  strsmatchs        - like strsmatch but outputs matches in different columns
%  strmatch_anywhr   - finds rows of "strs" which have string "str" anywhere
%  strmatch_end      - like strmatch but looks from end of each row
%  cell2str          - convert a cell array to a string array
%
%
% Other:
%  det_poly
%  dijkstra          - Calculates shortest path between 2 vertices in graph
%  find_near         - find index of vector element with value nearest to input 
%  gfile_shot_time   - return shot and time from gfile or afile name
%  intersect_bdry    - find intersection of line with curve
%  search_path       - string search of functions in matlabpath
%  submatrix_poly
%  tok_from_pwd      - determines tokamak string based on "pwd" and runs startup
%  whatfunction      - find function whose name contains partial string
%  zero_crossing     - estimate value(s) of x where data-defined fn. crosses y=0 

% do these need to be in this directory and/or contents file?
%dijkstra2.m
%dlmwrite.m
%error_data2.m
%find_latest_ver.m
%ineqconstr_lsqs.m
%intvec2str.m
%is_inuse.m
%keep.m
%keep_v5.m
%make_HP.m
%make_psc_file.m
%max_max.m
%min_min.m
%packnlstr.m
%modulo.m
%parse_filename.m
%pause_fig.m
%quadg.m
%read_lines.m
%regstrmatch.m
%removes_space.m
%rmchar.m
%save_whos_data.m
%set_mat_v6.m
%slope.m
%strcombine.m
%stringlen.m
%strwcmp.m
