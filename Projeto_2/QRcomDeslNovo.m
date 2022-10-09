%% Algoritmo QR com deslocamento otimizado
clc; clear;

function saida = sinal(parametro)
    if parametro >= 0
        saida = 1;
        return
    endif
    saida = -1;
endfunction;

disp("Apenas para matriz tridiagonal simétrica");
tamanho = input("Qual o tamanho da matriz: ");
precisao = 1E-6;

disp("\nDiagonal inferior:");
for i = 1:1:tamanho-1
%%    printf("Elemento a%d%d:",i+1,i)
%%    diagInf(i) = input(" ");
    diagInf(i) = -1;
endfor

disp("\nDiagonal principal:\n");
for i = 1:1:tamanho
%%    printf("Elemento a%d%d:",i,i)
%%    diagPri(i) = input(" ");
    diagPri(i) = 2;
endfor

diagSup = diagInf;
A = full(gallery("tridiag",diagInf,diagPri,diagSup));

V = eye(tamanho);
iter = 0;
mi = 0;

for m = tamanho:-1:2
    while abs(diagInf(m-1)) > precisao
        if iter != 0
            d  = (diagPri(m-1) - diagPri(m))/2;
            mi = diagPri(m) + d - sinal(d)*sqrt(d^2 + diagInf(m-1)^2);
        endif
        
        diagPri = diagPri - mi;
        
        for i = 1:1:tamanho-1
            if abs(diagPri(i)) > abs(diagInf(i))
                tau  = -diagInf(i)/diagPri(i);
                coss(i) = 1/sqrt(1+tau^2);
                seno(i) = coss(i)*tau;
            else
                tau  = -diagPri(i)/diagInf(i);
                seno(i) = 1/sqrt(1+tau^2);
                coss(i) = seno(i)*tau;
            endif
            
            diagPri(i  ) = coss(i)*diagPri(i  ) - seno(i)*diagInf(i  );
            diagSupExtra = diagSup(i);
            diagSup(i  ) = coss(i)*diagSupExtra - seno(i)*diagPri(i+1);
            diagPri(i+1) = seno(i)*diagSupExtra + coss(i)*diagPri(i+1);
            if i != tamanho-1
                diagSup(i+1) = coss(i)*diagSup(i+1);
            endif
        endfor
        
        for i = 1:1:tamanho-1
            diagInf(i  ) = -seno(i)*diagPri(i+1);
            diagPri(i  ) =  coss(i)*diagPri(i  ) - seno(i)*diagSup(i  );
            diagPri(i+1) =  coss(i)*diagPri(i+1);
            diagSup(i  ) =  diagInf(i);
        endfor
        
        diagPri = diagPri + mi;
        
        for j = 1:1:tamanho-1
            for i = 1:1:tamanho
                Vextra   = V(i,j);
                V(i,j  ) = coss(j)*Vextra - seno(j)*V(i,j+1);
                V(i,j+1) = seno(j)*Vextra + coss(j)*V(i,j+1);
            endfor
        endfor
        
        iter = iter +1;
    endwhile
    diagInf(m-1) = 0;
    diagSup(m-1) = 0;
endfor
autovalores = diagPri;
iter