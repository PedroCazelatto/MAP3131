%% Metodo da potencia
clc; clear;

precisao = input("Qual a precisao: ");
tamanho = input("Qual o tamanho da matriz: ");

for i = 1:1:tamanho
    for j = i:1:tamanho
        printf("Elemento a%d%d:",i,j)
        matriz(i,j) = input(" ");
        matriz(j,i) = matriz(i,j);
    endfor
endfor

for i = 1:1:tamanho
    printf("Elemento v%d:",i)
    vAnt(i,1) = input(" ");
endfor

for k = 1:1:tamanho
    vAntNorm = vAnt ./ norm(vAnt);
    autovalorAnt = 0;
    diferenca = inf;
    
    while diferenca > precisao
        vNovo = matriz * vAntNorm;;
        vNovoNorm = vNovo ./ norm(vNovo);
        autovalorNovo = transpose(vNovoNorm)*matriz*vNovoNorm;
        diferenca = abs(autovalorNovo - autovalorAnt);
        autovalorAnt = autovalorNovo;
        vAntNorm = vNovoNorm;
    endwhile
    
    printf("O %do maior autovalor é: %f\n",k,autovalorNovo);
    matriz = matriz - autovalorNovo*vNovoNorm*transpose(vNovoNorm);
endfor