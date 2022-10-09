%% Tranformacao de Householder
clc; clear;

function saida = sinal(parametro)
    if parametro >= 0
        saida = 1;
        return
    endif
    saida = -1;
endfunction;

tamanho = 4;
entrada = [2,4,1,1;4,2,1,1;1,1,1,2;1,1,2,1];

disp("Apenas para matriz simétrica");
%%tamanho = input("Qual o tamanho da matriz: ");
%%for i = 1:1:tamanho
%%    printf("\n%da coluna:\n", i);
%%    for j = i:1:tamanho
%%        printf("Elemento a%d%d:",j,i)
%%        entrada(j,i) = input(" ");
%%        entrada(i,j) = entrada(j,i);
%%    endfor
%%endfor

Ht = eye(tamanho);
entrada

for i = 1:1:tamanho-2
    % Calculo do vetor omega
    for j = 1:1:tamanho-i
        w(j) = entrada(j+i,i);
    endfor
    w(1) = w(1) + norm(w)*sinal(w(1));
    ww = dot(w,w);
    % Transforma as colunas da matriz de entrada
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
    % Transforma as linhas da matriz de entrada
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
    % Transforma as linhas da matriz Ht
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
    w(tamanho-i) = 0; % Zera o ultimo elemento de
    x(tamanho-i) = 0; % cada vetor utilizado, pois
    h(tamanho-i) = 0; % sera reusado com n-1 elementos
endfor

entrada
Ht
%%
%%for i = 1:1:tamanho-1
%%    diagInf(i) = entrada(i+1,i  );
%%    diagPri(i) = entrada(i  ,i  );
%%    diagSup(i) = entrada(i  ,i+1);
%%endfor
%%diagPri(tamanho) = entrada(tamanho,tamanho);