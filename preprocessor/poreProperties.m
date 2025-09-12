function [Kn, locporesize, pore, BW3] = poreProperties(pore, lamb, Cl)
    %% Initialization   
    [ny, nx] = size(pore);
    diff = 1.0;
    cleaning_attempt = 1;
    while diff > 0
        disp(strcat('====================== Smoothing pores-', num2str(cleaning_attempt)));
        cleaning_attempt = cleaning_attempt + 1;
        pore_cp = pore;
    
        % local poresize
        % Extend the length of the channel on the boundary to the outside of the boundary
        w_bnd = floor(max(size(pore)));
        pore_expand = ones(ny+2*w_bnd,nx+2*w_bnd);
        pore_expand(w_bnd+1:w_bnd+ny,w_bnd+1:w_bnd+nx) = pore;
        % for i = 1:ny
        %     %left boundary
        %     if pore(i,1)
        %         pore_expand(w_bnd+i,1:w_bnd) = 1;
        %     end
        %     %right boundary
        %     if pore(i,nx)
        %         pore_expand(w_bnd+i,w_bnd+nx+1:end) = 1;
        %     end
        % end
        pore_expand(1, :) = 0;
        pore_expand(end, :) = 0;
        pore_expand(:, 1) = 0;
        pore_expand(:, 1) = 0;
        BW = pore_expand;
        BW1 = bwmorph(BW,'thin',Inf);
        BW2 = BW1(w_bnd+1:w_bnd+ny,w_bnd+1:w_bnd+nx);   % midaxis
        [D,~] = bwdist(~pore);  % distanct to wall
        BW3 = D.*BW2;           % midaxis to wall (local poresize of midaxis)
        [~,L] = bwdist(BW2);    % index of closet midaxis
        locporesize = double(BW3);
        locpore_midaxis = int64(pore).*int64(L);
        [row,col] = find(locpore_midaxis);
        for i = 1:length(row)
            temp_row = row(i);
            temp_col = col(i);
            temp_midaxis = locpore_midaxis(temp_row,temp_col);
            locporesize(temp_row,temp_col) = BW3(temp_midaxis)*2;
        end
    
        Kn = lamb./locporesize(:, :)./Cl;
        Kn(Kn==inf)  = 0;
        
        % clean
        for i =1:ny
            for j =1:nx
                % the case where pore size is too small
                if pore(i, j) == 1 && locporesize(i, j) < 2
                    pore(i, j) = 0;
                end
            end
        end
        pore = clean_boundaries(pore);
        diff = sum(sum(abs(pore_cp - pore)));
    end
end

