clearvars
clearvars -GLOBAL
close all
set(0,'DefaultFigureWindowStyle', 'docked')

% 
nx = 40;
ny = 40;
max=3;
Lb = round(nx/max); % rounds the elements of X to the nearest integers
Wb = round(ny/max);


inside = 1;
outside = 10e-2;

%Conductivity map


conductivityMap = zeros(nx,ny);
% mapping of inside
for i = 1:nx
    for j = 1: ny
        conductivityMap(i,j) = inside;
    end
end
%mapping of outside
% for i = 1:nx
%     for j = 1: ny
%         conductivityMap(i,j) = outside;
%     end
% end
max1=2;
for i = 1:nx
    for j = 1:ny
     
        if (i>=1 && i<=Wb && j>Lb && j<=(max1*Lb) )||(i<=ny && i>=(ny-Wb) && j>Lb && j<=(max1*Lb));
            %max1*Lb
            conductivityMap(i,j) = outside;
        end
%         
        % if (i<=ny && i>=(ny-Wb) && j>Lb && j<=(max1*Lb))
           %  conductivityMap(i,j) = outside;
        %end
    end
end 
conductivityMap=conductivityMap';
G = sparse(nx*ny);
b = zeros(1, nx*ny);

% 

% G - Matrix Formulation from slides
for i = 1:nx
    for j = 1:ny
        n = j + (i - 1) * ny;
        
        
        if i == 1
            G(n,:) = 0;
            G(n,n) = 1;
             b(n) = 1;
            
        elseif i == nx %(i == 1 && i > 1 && i < nx)
             G(n,:) = 0;
             G(n,n) = 1;
        
        elseif j == 1 
            nxm = j + (i-2) * ny;
            nxp = j + (i) * ny;
            nyp = (j+1) + (i-1) * ny;
            
            rxm = (conductivityMap(i,j) + conductivityMap(i-1,j))/2.0;
            rxp = (conductivityMap(i,j) + conductivityMap(i+1,j))/2.0;
            ryp = (conductivityMap(i,j) + conductivityMap(i,j+1))/2.0;
            
            G(n,n) = -(rxm + rxp + ryp); 
            G(n,nxm) = rxm;
            G(n,nxp) = rxp;
            G(n, nyp) = ryp;
              
        elseif j == ny %(j == ny && i > 1 && i < nx)
            nxm = j + (i-2) * ny;
            nxp = j + (i) * ny;
            nym = (j-1) + (i-1) * ny;
            
            rxm = (conductivityMap(i,j) + conductivityMap(i-1,j))/2.0;
            rxp = (conductivityMap(i,j) + conductivityMap(i+1,j))/2.0;
            rym = (conductivityMap(i,j) + conductivityMap(i,j-1))/2.0;
            
            G(n,n) = -(rxm + rxp + rym);
            G(n,nxm) = rxm;
            G(n,nxp) = rxp;
            G(n, nym) = rym;
      
        else
            nxm = j + (i-2) * ny;
            nxp = j + (i) * ny;
            nym = (j-1) + (i-1) * ny;
            nyp = (j+1) + (i-1) * ny;
            
            rxm = (conductivityMap(i,j) + conductivityMap(i-1,j))/2.0;
            rxp = (conductivityMap(i,j) + conductivityMap(i+1,j))/2.0;
            rym = (conductivityMap(i,j) + conductivityMap(i,j-1))/2.0;
            ryp = (conductivityMap(i,j) + conductivityMap(i,j+1))/2.0;
            
            G(n,n) = -(rxm + rxp + rym + ryp);
            G(n,nxm) = rxm;
            G(n,nxp) = rxp;
            G(n, nym) = rym;
            G(n, nyp) = ryp;
        end
    end
end
V = G\b';
%vmap= zeros (nx,ny);
for i = 1:nx %Converting V to matrix to plot
    for j = 1:ny
        n = j + (i-1) * ny;
        VG(i,j) = V(n);
    end   
end
for i = 1:nx
    for j = 1:ny
        if i == 1
            Ex(i, j) = (VG(i + 1, j) - VG(i, j));
        elseif i == nx
            Ex(i, j) = (VG(i, j) - VG(i - 1, j));
        else
            Ex(i, j) = (VG(i + 1, j) - VG(i - 1, j)) * 0.5;
        end
        if j == 1
            Ey(i, j) = (VG(i, j + 1) - VG(i, j));
        elseif j == ny
            Ey(i, j) = (VG(i, j) - VG(i, j - 1));
        else
            Ey(i, j) = (VG(i, j + 1) - VG(i, j - 1)) * 0.5;
        end
    end
end

Ex = -Ex;
Ey = -Ey;

figure(2);
surf(VG);
title('Voltage Map with Bottleneck');
figure(3); 
quiver(Ex', Ey');
title('Electric field Plot')