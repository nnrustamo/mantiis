function img = porous_media_circles(N, num_circles, radius_range)
    % Create a porous media image with circles representing obstacles.
    % N: size of the image (NxN)
    % num_circles: number of circles to place in the image
    % radius_range: a two-element vector [rmin, rmax] for radius range
    
    img = ones(N, N);
    img(1, :) = 0;
    img(end, :) = 0;
    
    % Randomly place non-overlapping circles
    obstacle_mask = false(N, N);
    
    max_attempts = 10000;
    attempts = 0;
    
    while num_circles > 0 || attempts < max_attempts
        % Randomly choose a radius within the specified range.
        radius = randi(radius_range);

        % Randomly choose a center for the circle.
        x_center = randi([1, N]);
        y_center = randi([2, N-1]);
        
        % Create the circle using a mask.
        [xx, yy] = meshgrid(1:N, 1:N);
        circle_mask = sqrt((xx - x_center).^2 + (yy - y_center).^2) <= radius;
        
        % Check if this circle overlaps with existing obstacles.
        if ~any(obstacle_mask(circle_mask))
            obstacle_mask = obstacle_mask | circle_mask;
            num_circles = num_circles - 1;
        end
        
        attempts = attempts + 1;
        fprintf('circles left %d\n', num_circles);
        fprintf('number of attempts %d\n', attempts);
    end
    
    img(obstacle_mask) = 0;

    % Add buffer
    img(2:end-1, 1:floor(0.1*N)) = 1;
    img(2:end-1, N-floor(0.1*N):end) = 1;
    
    % Display the generated image.
    figure;
    imshow(img);
    title('Porous Media Image');

end
