function data = loadTxtFile(filename)
    if exist(filename, 'file') ~= 2
        error('File does not exist.');
    end
    fileID = fopen(filename, 'r');
    data = fscanf(fileID, '%f');
    fclose(fileID);
end
