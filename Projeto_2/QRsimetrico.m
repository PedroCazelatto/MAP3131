%% Algoritmo QR com deslocamento para qualquer matriz simetrica
clc; clear;

function saida = sinal(parametro)
    if parametro >= 0
        saida = 1;
        return
    endif
    saida = -1;
endfunction;

precisao = 1E-6;
disp("Apenas para matriz simétrica");

tamanho = 4;
entrada = [2,4,1,1;4,2,1,1;1,1,1,2;1,1,2,1]; %% Teste A
%%entrada = [2,-1,1,3;-1,1,4,2;1,4,2,-1;3,2,-1,1];

%%for i = 1:1:tamanho %% Teste B
%%    for j = i:1:tamanho
%%        entrada(j,i) = tamanho+1-max(i,j);
%%        entrada(i,j) = entrada(j,i);
%%    endfor
%%endfor

%%tamanho = input("Qual o tamanho da matriz: ");
%%for i = 1:1:tamanho
%%    printf("\n%da coluna:\n", i);
%%    for j = i:1:tamanho
%%        printf("Elemento a%d%d:",j,i)
%%        entrada(j,i) = input(" ");
%%        entrada(i,j) = entrada(j,i);
%%    endfor
%%endfor

A = entrada;
Ht = eye(tamanho);

for i = 1:1:tamanho-2
    for j = 1:1:tamanho-i
        w(j) = entrada(j+i,i);
    endfor
    
    w(1) = w(1) + norm(w)*sinal(w(1));
    ww = dot(w,w);
    
    for j = i:1:tamanho
        for k = 1:1:tamanho-i
            x(k) = entrada(k+i,j);
        endfor
        
        wx = dot(w,x);
        alfax = -2*wx/ww;
        
        for k = 1:1:tamanho-i
            entrada(k+i,j) = entrada(k+i,j) + alfax*w(k);
        endfor
    endfor
    
    for j = i:1:tamanho
        for k = 1:1:tamanho-i
            x(k) = entrada(j,k+i);
        endfor
        
        wx = dot(w,x);
        alfax = -2*wx/ww;
        
        for k = 1:1:tamanho-i
            entrada(j,k+i) = entrada(j,k+i) + alfax*w(k);
        endfor
    endfor
    
    for j = 2:1:tamanho
        for k = 1:1:tamanho-i
            h(k) = Ht(j,k+i);
        endfor
        
        wh = dot(w,h);
        alfah = -2*wh/ww;
        
        for k = 1:1:tamanho-i
            Ht(j,k+i) = Ht(j,k+i) + alfah*w(k);
        endfor
    endfor
    
    w(tamanho-i) = 0;
    x(tamanho-i) = 0;
    h(tamanho-i) = 0;
endfor

entrada

for i = 1:1:tamanho-1
    diagSup(i) = entrada(i+1,i  );
    diagPri(i) = entrada(i  ,i  );
    diagInf(i) = entrada(i  ,i+1);
endfor
diagPri(tamanho) = entrada(tamanho,tamanho);

V = Ht;

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
autovalores = diagPri
V;
iter;