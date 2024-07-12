function pore = simplePore(ny, nx)
    pore = ones(ny, nx); 
    % circle
    ic = 256; % floor(ny/2); 
    jc = 256; %floor(nx/2);
    r = 50;
    for i = 1:ny
        for j = 1:nx
            if sqrt((i-ic)^2 +(j - jc)^2) <=r
                pore(i, j) = 0;
            end
        end
    end
    
%     ic = 700; % floor(ny/2); 
%     jc = 700; %floor(nx/2);
%     for i = 1:ny
%         for j = 1:nx
%             if sqrt((i-ic)^2 +(j - jc)^2) <=r
%                 pore(i, j) = 0;
%             end
%         end
%     end

    pore(1, :) = 0;
    pore(end, :) = 0;
end