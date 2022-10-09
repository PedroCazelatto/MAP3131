%% Fatoracao QR
clc; clear;

disp("Apenas para matriz tridiagonal");
tamanho = input("Qual o tamanho da matriz: ");

disp("\nDiagonal inferior:");
for i = 1:1:tamanho-1
    printf("Elemento a%d%d:",i+1,i)
    diagInf(i) = input(" ");
endfor

disp("\nDiagonal principal:\n");
for i = 1:1:tamanho
    printf("Elemento a%d%d:",i,i)
    diagPri(i) = input(" ");
endfor

disp("\nDiagonal superior:\n");
for i = 1:1:tamanho-1
    printf("Elemento a%d%d:",i,i+1)
    diagSup(i) = input(" ");
endfor

Qinv = eye(tamanho);
A = full(gallery("tridiag",diagInf,diagPri,diagSup));

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
    
    diagPri(i+1) = seno(i)*diagSup(i) + coss(i)*diagPri(i+1);
    if i != tamanho-1
        diagSup(i+1) = diagSup(i+1)*coss(i);
    endif
    
    aux = eye(tamanho);
    aux( i , i ) =  coss(i);
    aux(i+1,i+1) =  coss(i);
    aux( i ,i+1) = -seno(i);
    aux(i+1, i ) =  seno(i);
    
    Qinv = aux*Qinv;
    
endfor

Q = transpose(Qinv)
R = Qinv*A

