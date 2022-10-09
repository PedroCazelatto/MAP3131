%% Algoritmo QR com deslocamento
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
%%precisao = input("Qual a precisao: ");
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

%%disp("\nDiagonal superior:\n");
%%for i = 1:1:tamanho-1
%%    printf("Elemento a%d%d:",i,i+1)
%%    diagSup(i) = input(" ");
%%endfor
diagSup = diagInf;

inicial = full(gallery("tridiag",diagInf,diagPri,diagSup));

A = inicial;
V = eye(tamanho);
iter = 0;
mi(1) = 0;

for m = tamanho:-1:2
    while abs(A(m,m-1)) > precisao
        if iter != 0
            d = (A(m-1,m-1) - A(m,m))/2;
            mi(iter+1) = A(m,m) + d - sinal(d)*sqrt(d^2+A(m,m-1)^2);
        endif
        Aprox = A - mi(iter+1)*eye(tamanho);
        Qinv = eye(tamanho);
        for i = 1:1:tamanho-1
            if abs(Aprox(i,i)) > abs(Aprox(i+1,i))
                tau  = -Aprox(i+1,i)/Aprox(i,i);
                coss(i) = 1/sqrt(1+tau^2);
                seno(i) = coss(i)*tau;
            else
                tau  = -Aprox(i,i)/Aprox(i+1,i);
                seno(i) = 1/sqrt(1+tau^2);
                coss(i) = seno(i)*tau;
            endif
            Aprox(i+1,i+1) = seno(i)*Aprox(i,i+1) + coss(i)*Aprox(i+1,i+1);
            if i != tamanho-1
                Aprox(i+1,i+2) = Aprox(i+1,i+2)*coss(i);
            endif
            aux = eye(tamanho);
            aux( i , i ) =  coss(i);
            aux(i+1,i+1) =  coss(i);
            aux( i ,i+1) = -seno(i);
            aux(i+1, i ) =  seno(i);
            Qinv = aux*Qinv;
        endfor
        Q = transpose(Qinv);
        R = Qinv*(A - mi(iter+1)*eye(tamanho));
        A = R*Q + mi(iter+1)*eye(tamanho);
        V = V*Q;
        for x = 1:1:tamanho
            for y = 1:1:tamanho
                if abs(A(x,y)) < precisao
                    A(x,y) = 0;
                endif
            endfor
        endfor
        iter = iter +1
    endwhile
endfor
autovalores = diag(A);

autovalores;
V;
iter
