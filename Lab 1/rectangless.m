function  rectangles = rectangless(xVec,f)
h = xVec(2) - xVec(1);
rect = h*f((xVec(1:end-1)+xVec(2:end))/2);
rect0 = [0.0];
rectangles = cumsum([rect0 rect]);
end

%(xVec(1:end-1)+xVec(2:end))/2