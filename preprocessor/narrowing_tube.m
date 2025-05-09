function pore = narrowing_tube(N, maxWidth, minWidth)

    % Parameters
    imageSize = N;
    
    waveAmplitude = 10; % Amplitude of the wave deformation
    waveFrequency = 0.006; % Decreased frequency of the wave deformation
    binaryImage = zeros(imageSize, imageSize);

    for y = 1:imageSize
        % Compute the width of the tube at position y with a wavy deformation
        baseWidth = maxWidth - ((maxWidth - minWidth) * (y / imageSize));
        wave = waveAmplitude * sin(2 * pi * waveFrequency * y);
        currentWidth = baseWidth + wave;
        
        % Calculate the range of x-coordinates that are inside the tube
        xMin = round((imageSize - currentWidth) / 2);
        xMax = round((imageSize + currentWidth) / 2);
        xMin = max(1, xMin);
        xMax = min(imageSize, xMax);
        
        % Set pixels inside the tube to 1
        binaryImage(y, xMin:xMax) = 1;
    end

    pore = binaryImage';
    pore(:, 1:50) = 1;
    pore(:, end-51:end) = 1;
    pore(1, :) = 0;
    pore(end, :) = 0;

end
