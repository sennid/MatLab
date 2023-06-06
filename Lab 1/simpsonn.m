function simpson = simpsonn(xVec, f)
h = xVec(2) - xVec(1);
simp = h/6*(f(xVec(1:end-1))+4*f((xVec(1:end-1)+h/2))+f(xVec(1:end-1)+h));
simp0 = [0.0];
simpson = cumsum([simp0 simp]);
end