function [value, isterminal, direction] = funcEvent6(t,x,ax,bx,ay,by)
    value = (x(1)-ax).*(x(1)-bx).*(x(3)-ay).*(x(3)-by);
    isterminal = 1;
    direction = 0;