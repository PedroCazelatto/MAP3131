%% Teste da trelica sem otimizacao
clc; clear;

function saida = sinal(parametro)
    if parametro >= 0
        saida = 1;
        return
    endif
    saida = -1;
endfunction;

arquivo = importdata("input-c.txt", " ");

Total    = arquivo(1,1);
livres   = arquivo(1,2);
barras   = arquivo(1,3);
tamanho  = 2*livres;
precisao = 1E-6;

rho = arquivo(2,1);
Are = arquivo(2,2); % Area das barras
E   = arquivo(2,3)*1E9; % Elasticidade
x   = zeros(tamanho,1); % Deslocamento dos nós, horizontal nos impares e vertical nos pares
M   = zeros(tamanho,1);
K   = zeros(tamanho);
top = Are*E;

%% Gera o K
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
    mult = top/L;
    
    for x = 1:1:4
        if pos(x) <= tamanho
            M(pos(x),1) = M(pos(x),1) + 0.5*rho*Are*L;
        endif
        for y = 1:1:4
            somar = 1;
            somar = somar*coss^rem(x,2)*coss^rem(y,2);
            somar = somar*seno^rem(x+1,2)*seno^rem(y+1,2);
            somar = somar*sinal(x-2.5)*sinal(y-2.5);
            somar = somar*mult;
            
            if and((pos(x) <= tamanho),(pos(y) <= tamanho))
                K(pos(x), pos(y)) = K(pos(x), pos(y)) + somar;
            endif
        endfor
    endfor
    
endfor


%% Calcula M^-1/2
for x = 1:1:tamanho
    Mt(x,1) = M(x,1);
    M(x,1) = 1/sqrt(M(x,1));
endfor

%% Gera o Ktil
for x = 1:1:tamanho
    for y = 1:1:tamanho
        Ktil(x,y) = M(x,1)*K(x,y)*M(y,1);
    endfor
endfor

% ok


%% Tridiagonaliza a matriz Ktil
entrada = Ktil;
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
            h(k) = Ht(j,k+i);
        endfor
        
        wx = dot(w,x);
        alfax = -2*wx/ww;
        wh = dot(w,h);
        alfah = -2*wh/ww;
        
        for k = 1:1:tamanho-i
            entrada(k+i,j) = entrada(k+i,j) + alfax*w(k);
            Ht(j,k+i) = Ht(j,k+i) + alfah*w(k);
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
    w(tamanho-i) = 0;
endfor

%% Aqui, entrada eh tridiagonal de Ktil

%% Separa as tres diagonais
for i = 1:1:tamanho-1
    diagSup(i) = entrada(i+1,i  );
    diagPri(i) = entrada(i  ,i  );
    diagInf(i) = entrada(i  ,i+1);
endfor
diagPri(tamanho) = entrada(tamanho,tamanho);


%% Calcula os autovalores e autovetores;
V = eye(tamanho);
iter = 0;
mi = 0;

A = entrada;
for m = tamanho:-1:2
    while abs(A(m,m-1)) > precisao
%%      Calcula o mi
        if iter != 0
            d = (A(m-1,m-1) - A(m,m))/2;
            mi = A(m,m) + d - sinal(d)*sqrt(d^2+A(m,m-1)^2);
        endif
        
        Aprox = A - mi*eye(tamanho);
        
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
        R = Qinv*(A - mi*eye(tamanho));
        A = R*Q + mi*eye(tamanho);
        V = V*Q;
        iter = iter +1;
    endwhile
    A(m  ,m-1) = 0;
    A(m-1,m  ) = 0;
endfor
autovalores = diag(A);


%% Separa as menores frequencias
%%frequencias = sqrt(autovalores);
%%for i = 1:1:5
%%    aux = frequencias(i);
%%    posicao(i) = i;
%%    for j = (i+1):1:tamanho
%%        if frequencias(j) < aux
%%            posicao(i) = j;
%%            aux = frequencias(j);
%%        endif
%%    endfor
%%    frequencias(posicao(i)) = frequencias(i);
%%    frequencias(i) = aux;
%%endfor
%%frequencias;
%%posicao;