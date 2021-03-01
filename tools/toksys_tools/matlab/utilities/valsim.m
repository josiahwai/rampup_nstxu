function valsim(ptns, source_array, dia, output_ps)
% valsim(source_array)
%
% INPUT
%   ptns = strlist object containing the names of the points to plot
%   source_array = cell array of data_source objects
%   dia = txttab object
%   output_ps = (OPTIONAL) File name of output ps file.  If specified, plots will be written to a %               *.ps file specified by this argument.


no_alias = strlist;

pldat_group_array = [];
ptns_length = get_length(ptns);
for ptn_index = 1:ptns_length
    pldat_array = [];
    for source_index = 1:length(source_array)
        src = source_array{source_index};
        
        diary_field = get_source_type(src);
        ptnlist_type = get_list_type(ptns);

        ptn_name = get_element_by_index(ptns,ptn_index);
 
        alias = get_alias(dia, ptnlist_type, ptn_name, diary_field); 
        
        if ~strcmp(alias, 'NO_MATCH')
            pld = pldat(src, alias);
            pldat_array = [pldat_array pld];
        else
            no_alias = add_element(no_alias, ptn_name);
        end
    end
    current_group = pldat_group(pldat_array);
    pldat_group_array = [pldat_group_array current_group];
end

%pldat_group_array_len = length(pldat_group_array)
plot_man = plot_manager(pldat_group_array);
plot_man = set_plots_per_page(plot_man, 3);
plot_man = set_stack_plots(plot_man, 1);

if exist('output_ps', 'var')
    create_plot(plot_man, output_ps);
else
    create_plot(plot_man);
end

if (get_length(no_alias) > 0)

    fprintf('\n\n\nWarning: No aliases were found for the following points.\n');

    for j=1:get_length(no_alias)
        fprintf('%s\n', get_element_by_index(no_alias,j));    
    end
end



