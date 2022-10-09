%%  Teste da trelica com otimizacao
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

rho = arquivo(2,1);     % Densidade volumetrica das barras
A   = arquivo(2,2);     % Area das barras
E   = arquivo(2,3)*1E9; % Elasticidade
M   = zeros(tamanho,1); % Vetor de massas
K   = zeros(tamanho);   % Matriz de rigidez total
top = A*E;

%% Gera o K e monta a M
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
            M(pos(x),1) = M(pos(x),1) + 0.5*rho*A*L;
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
    M(x,1) = 1/sqrt(M(x,1));
endfor

%% Gera o Ktil no K
for x = 1:1:tamanho
    for y = 1:1:tamanho
        K(x,y) = M(x,1)*K(x,y)*M(y,1);
    endfor
endfor

%% Tridiagonaliza a matriz K
entrada = K;
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

%% Separa as tres diagonais
for i = 1:1:tamanho-1
    diagSup(i) = entrada(i+1,i  );
    diagPri(i) = entrada(i  ,i  );
    diagInf(i) = entrada(i  ,i+1);
endfor
diagPri(tamanho) = entrada(tamanho,tamanho);

%% Calcula os autovalores e autovetores;
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
autovalores = transpose(diagPri);

%% Separa as menores frequencias
anterior = 0;
for i = 1:1:5
    minFreq(i) = autovalores(1);
    posicao(i) = 1;
    for j = 2:1:tamanho
        if and(autovalores(j) < minFreq(i), anterior < autovalores(j))
            minFreq(i) = autovalores(j);
            posicao(i) = j;
        endif
    endfor
    anterior = minFreq(i);
endfor

minFreq = sqrt(minFreq);

%% Junta os modos das menores frequencias
for i = 1:1:5
    for j = 1:1:tamanho
        minModo(j,i) = M(j,1)*V(j,posicao(i));
    endfor
endfor
