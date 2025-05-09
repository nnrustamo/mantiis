function savetoFile(P, fname)
    [ny, nx] = size(P);
    P = reshape(P, [ny*nx,1]);
    writematrix(P, fname);
    fprintf('file has been printed \n')
end