function img = porous_media_circles(N, num_circles, radius_range)

    % Parameters
    buffer_const = 10; % extra spacing margin

    % Initialize image and masks
    img = ones(N, N);
    img(1, :) = 0;
    img(end, :) = 0;

    obstacle_mask = false(N, N);
    [xx, yy] = meshgrid(1:N, 1:N);

    % Candidate centers: initially all valid pixels
    max_r = radius_range(2);
    option_mask = true(N, N);
    option_mask(1:max_r, :) = false;
    option_mask(end-max_r+1:end, :) = false;
    option_mask(:, 1:max_r) = false;
    option_mask(:, end-max_r+1:end) = false;

    % Radii sorted largest → smallest
    radii_list = radius_range(2):-1:radius_range(1);

    attempts = 0;

    while num_circles > 0 && any(option_mask(:)) && ~isempty(radii_list)

        placed = false;

        % Try radii from largest to smallest
        for r_idx = 1:length(radii_list)
            radius = radii_list(r_idx);
            buffer = radius + buffer_const;

            % Compute feasible centers: option_mask minus buffered obstacles
            exclusion_from_obstacles = imdilate(obstacle_mask, strel('disk', buffer));
            feasible_mask = option_mask & ~exclusion_from_obstacles;

            if ~any(feasible_mask(:))
                % Radius not possible anymore → log and remove from list
                fprintf('[INFO] Radius %d removed (no feasible space left).\n', radius);
                radii_list(r_idx) = NaN; % mark for removal
                continue
            end

            % --- Pick random feasible center ---
            idx = find(feasible_mask);
            k = idx(randi(numel(idx)));
            [y_center, x_center] = ind2sub(size(feasible_mask), k);

            % Place circle
            circle_mask = sqrt((xx - x_center).^2 + (yy - y_center).^2) <= radius;
            obstacle_mask = obstacle_mask | circle_mask;
            num_circles = num_circles - 1;

            % Remove option space with buffer
            exclusion_mask = sqrt((xx - x_center).^2 + (yy - y_center).^2) <= buffer;
            option_mask(exclusion_mask) = false;

            fprintf('[PLACED] Circle radius %d at (%d, %d). Circles left: %d\n', ...
                    radius, x_center, y_center, num_circles);

            placed = true;
            break % stop after placing one circle
        end

        % Clean up removed radii
        radii_list = radii_list(~isnan(radii_list));

        attempts = attempts + 1;
        if mod(attempts, 50) == 0
            fprintf('Attempts: %d, Remaining radii: %s\n', attempts, mat2str(radii_list));
        end

        % If no placement was possible in this round, break
        if ~placed && isempty(radii_list)
            fprintf('[STOP] No more feasible radii left.\n');
            break
        end
    end

    img(obstacle_mask) = 0;

    % Keep vertical inlet/outlet channels
    img(2:end-1, 1:floor(0.1*N)) = 1;
    img(2:end-1, N-floor(0.1*N):end) = 1;

    % figure; imshow(img);
    % title('Porous Media Image');

end
