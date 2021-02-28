function [n] = recCoordinates(i,j, nx, ny)
%RECCOORDINATES maps between rectangular and linear cordinates
    global C;
    n = j + (i - 1)*ny;
end