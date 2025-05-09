function binaryMatrix = createRectangleObstacleFlowGeom(H)
    binaryMatrix = ones(H, H);
    center = ceil(H / 2);
    leftedge = ceil(center - 0.1*H);
    rightedge = ceil(center + 0.1*H);
    binaryMatrix(leftedge:rightedge, leftedge:rightedge) = 0;
    binaryMatrix(1, :) = 0;
    binaryMatrix(end, :) = 0;
end
