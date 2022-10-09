%% Rotacao
clc; clear;
format long;
disp("Apenas para matriz tridiagonal");
tamanho = input("Qual o tamanho da matriz: ");

for i = 1:1:tamanho-1
    matQ(:,:,i) = eye(tamanho);
endfor

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

Q = eye(tamanho);
A = full(gallery("tridiag",diagInf,diagPri,diagSup));
R = full(gallery("tridiag",diagInf,diagPri,diagSup));

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
    
    matQ(  i,   i, i) =  coss(i);
    matQ(i+1, i+1, i) =  coss(i);
    matQ(  i, i+1, i) = -seno(i);
    matQ(i+1,   i, i) =  seno(i);
    
    R = matQ(:,:,i)*R;
    
    diagInf = diag(R,-1);
    diagPri = diag(R, 0);
    diagSup = diag(R, 1);
    Q = Q*transpose(matQ(:,:,i));
endfor
format short;
