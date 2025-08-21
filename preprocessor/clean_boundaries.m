function pore = clean_boundaries(pore)
    boundary_file = "../include/system_files/bnd_types.txt";
    bnd_types = loadTxtFile(boundary_file);
    
    % add periodic ghost layers
    [ny, nx] = size(pore);
    pore_perio = zeros(ny + 2, nx+ 2);
    pore_perio(2:end-1, 2:end-1) = pore;
    % columns
    pore_perio(2:end-1, 1) = pore(:, end);
    pore_perio(2:end-1, end) = pore(:, 1);
    % rows
    pore_perio(1, 2:end-1) = pore(end, :);
    pore_perio(end, 2:end-1) = pore(1, :);
    pore_perio_cp = pore_perio;
    diff = 1;
    cleaning_attempt = 1;
    % loop thru the domain
    while diff > 0
        disp(strcat('Cleaning attemp- ', num2str(cleaning_attempt)));
        cleaning_attempt = cleaning_attempt + 1;
        for i = 2:ny + 1
            for j = 2:nx + 1
                kernel = pore_perio(i-1:i+1, j-1:j+1); % the kernel to determine the boundary types
                % check if it is a part of bnd_types
                if kernel(2, 2) == 1 && sum(sum(kernel)) < 9
                    kernel_lin = num2str(reshape(kernel, [9, 1]));
                    kernel_lin = sprintf('%s', kernel_lin);
                    decimalNumber = bin2dec(kernel_lin); % Convert binary string to decimal
                    if isempty(find(bnd_types == decimalNumber))
                        pore_perio(i, j) = 0;
                    end
                end
            end
        end
        diff = sum(sum(pore_perio_cp - pore_perio));
        pore_perio_cp = pore_perio;
    end
    
    % refresh pore
    pore = pore_perio(2:end-1, 2:end-1);    
end