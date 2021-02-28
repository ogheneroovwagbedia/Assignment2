function V = numericSolution(nx, ny, crecCoordinates, bc_left, bc_right, bc_top, bc_bottom)
%NUMERICSOLUTION Calculates the numeric solution using the finite difference method
%       Infinite BC means insolation
    global C;
    G = sparse(nx*ny, ny*nx);
    B = zeros(1, nx*ny);
    for i=1:nx
        for j=1:ny
            n = recCoordinates(i,j, nx, ny);
            nxm = recCoordinates(i-1,j, nx, ny);
            nxp = recCoordinates(i+1,j, nx, ny);
            nym = recCoordinates(i,j-1, nx, ny);
            nyp = recCoordinates(i,j+1, nx, ny);

            if (i == 1 && j == 1)
                if (bc_left == Inf)
                    rxp = (crecCoordinates(i,j) + crecCoordinates(i+1,j))/2;
                    ryp = (crecCoordinates(i,j) + crecCoordinates(i,j+1))/2;

                    G(n,n)   = -(rxp+ryp);
                    G(n,nxp) =  rxp;
                    G(n,nyp) =  ryp;
                else
                    G(n,n) = 1;
                    B(n) = bc_left;
                end
            elseif (i == 1 && j == ny)
                if (bc_left == Inf)
                    rxp = (crecCoordinates(i,j) + crecCoordinates(i+1,j))/2;
                    rym = (crecCoordinates(i,j) + crecCoordinates(i,j-1))/2;

                    G(n,n)   = -(rxp+rym);
                    G(n,nxp) =  rxp;
                    G(n,nym) =  rym;
                else
                    G(n,n) = 1;
                    B(n) = bc_left;
                end
            elseif i == nx && j == 1 % Right side
                if (bc_right == Inf)
                    rxm = (crecCoordinates(i,j) + crecCoordinates(i-1,j))/2;
                    ryp = (crecCoordinates(i,j) + crecCoordinates(i,j+1))/2;
                    G(n,n)   = -(rxm+ryp);
                    G(n,nxm) =  rxm;
                    G(n,nyp) =  ryp;
                else
                    G(n,n) = 1;
                    B(n) = bc_right;
                end
            elseif i == nx && j == ny % Right side
                if (bc_right == Inf)
                    rxm = (crecCoordinates(i,j) + crecCoordinates(i-1,j))/2;
                    rym = (crecCoordinates(i,j) + crecCoordinates(i,j-1))/2;
                    G(n,n)   = -(rxm+rym);
                    G(n,nxm) =  rxm;
                    G(n,nym) =  rym;
                else
                    G(n,n) = 1;
                    B(n) = bc_right;
                end
            elseif (i == 1) % Left Side
                if (bc_left == Inf)
                    rxp = (crecCoordinates(i,j) + crecCoordinates(i+1,j))/2;
                    rym = (crecCoordinates(i,j) + crecCoordinates(i,j-1))/2;
                    ryp = (crecCoordinates(i,j) + crecCoordinates(i,j+1))/2;

                    G(n,n)   = -(rxp+rym+ryp);
                    G(n,nxp) =  rxp;
                    G(n,nym) =  rym;
                    G(n,nyp) =  ryp;
                else
                    G(n,n) = 1;
                    B(n) = bc_left;
                end
            elseif i == nx % Right side
                if (bc_right == Inf)
                    rxm = (crecCoordinates(i,j) + crecCoordinates(i-1,j))/2;
                    rym = (crecCoordinates(i,j) + crecCoordinates(i,j-1))/2;
                    ryp = (crecCoordinates(i,j) + crecCoordinates(i,j+1))/2;
                    G(n,n)   = -(rxm+rym+ryp);
                    G(n,nxm) =  rxm;
                    G(n,nym) =  rym;
                    G(n,nyp) =  ryp;
                else
                    G(n,n) = 1;
                    B(n) = bc_right;
                end
            elseif j == 1 % Top Side
                if (bc_top == Inf)
                    rxm = (crecCoordinates(i,j) + crecCoordinates(i-1,j))/2;
                    rxp = (crecCoordinates(i,j) + crecCoordinates(i+1,j))/2;
                    ryp = (crecCoordinates(i,j) + crecCoordinates(i,j+1))/2;
                    G(n,n) = -(rxm+rxp+ryp);
                    G(n,nxm) =  rxm;
                    G(n,nxp) =  rxp;
                    G(n,nyp) =  ryp;
                else
                    G(n,n) = 1;
                    B(n) = bc_top;
                end
            elseif j == ny % Bottom Side
                if (bc_bottom == Inf)
                    rxm = (crecCoordinates(i,j) + crecCoordinates(i-1,j))/2;
                    rxp = (crecCoordinates(i,j) + crecCoordinates(i+1,j))/2;
                    rym = (crecCoordinates(i,j) + crecCoordinates(i,j-1))/2;
                    G(n,n) = -(rxm+rxp+rym);
                    G(n,nxm) =  rxm;
                    G(n,nxp) =  rxp;
                    G(n,nym) =  rym;
                else
                    G(n,n) = 1;
                    B(n) = bc_bottom;
                end
            else % Bulk Area
                rxm = (crecCoordinates(i,j) + crecCoordinates(i-1,j))/2;
                rxp = (crecCoordinates(i,j) + crecCoordinates(i+1,j))/2;
                rym = (crecCoordinates(i,j) + crecCoordinates(i,j-1))/2;
                ryp = (crecCoordinates(i,j) + crecCoordinates(i,j+1))/2;
                
                G(n,n)   = -(rxm+rxp+rym+ryp);
                G(n,nxm) =  rxm;
                G(n,nxp) =  rxp;
                G(n,nym) =  rym;
                G(n,nyp) =  ryp;
            end
        end
    end
    
    V_temp = G\B';
    
    V = zeros(nx,ny,1);
    for i=1:nx
        for j=1:ny
            V(i,j) = V_temp(recCoordinates(i,j, nx, ny));
        end
    end
end
