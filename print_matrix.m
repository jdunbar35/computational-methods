%% Prints matrices with desired # of decimal points, elements separated by tabs
% Jack Dunbar
% October 31, 2024

function [] = print_matrix(matrix, deci)
    % Create a format string for the desired number of decimal places
    formatSpec = ['%.', num2str(deci), 'f'];
    
    for i = 1:size(matrix, 1)
        fprintf(formatSpec, matrix(i, 1));  % Print first element
        for j = 2:size(matrix, 2)
            fprintf(['\t', formatSpec], matrix(i, j));  % Print remaining elements with tabs
        end
        fprintf('\n');  % New line after each row
    end
end
