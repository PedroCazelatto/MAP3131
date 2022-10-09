%%  Matrizes
clc; clear;

arquivo = importdata("input-c.txt", " ");

Total = arquivo(1,1);
livres = arquivo(1,2);
barras = arquivo(1,3);

rho = arquivo(2,1);
A = arquivo(2,2); % Area das barras
E = arquivo(2,3); % 

x = zeros(2*livres,1); % Deslocamento dos nós, horizontal nos impares e vertical nos pares
K = zeros(2*livres);

M = zeros(2*livres,1);

for iter = 1:1:barras
    i   = arquivo(iter+2,1);
    j   = arquivo(iter+2,2);
    ang = arquivo(iter+2,3);
    L   = arquivo(iter+2,4);
    
    pos(1) = 2*i-1;
    pos(2) = 2*i;
    pos(3) = 2*j-1;
    pos(4) = 2*j;
    
    coss = cosd(ang);
    seno = sind(ang);
    
    for x = 1:1:4
        if pos(x) <= 2*livres
            M(pos(x),1) = M(pos(x),1) + 0.5*rho*A*L;
        endif
        for y = 1:1:4
            somar = coss^rem(x,2)*coss^rem(y,2)*seno^rem(x+1,2)*seno^rem(y+1,2)*(-1)^floor(x/2.1)*(-1)^floor(y/2.1);
            if and((pos(x) <= 2*livres),(pos(y) <= 2*livres))
                K(pos(x), pos(y)) = K(pos(x), pos(y)) + somar;
            endif
        endfor
    endfor
    
endfor

for x = 1:1:2*livres
    M(x,1) = 1/sqrt(M(x,1));
endfor

for x = 1:1:2*livres
    for y = 1:1:2*livres
        Ktil(x,y) = M(x,1)*K(x,y)*M(y,1);
    endfor
endfor

Ktil
