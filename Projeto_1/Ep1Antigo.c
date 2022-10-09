/***************************************************************/
/**                                                           **/
/**    Pedro Henrique Galhardi Cazelatto - n° USP 11261090    **/
/**    Exercício-Programa 01 - MAP3121                        **/
/**    Autovalores e Autovetores de Matrizes                  **/
/**    Tridiagonais Simetricas - O Algoritmo QR               **/
/**                                                           **/
/***************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

/* Define o valor máximo das matrizes do codigo */
#define max 128

/* Calcula o sinal de um número, 1 para positivos e -1 para negativos */
int Sinal(double valor) {
    return 1-2*(valor < 0);
}

/* Imprime uma matriz quadrada de tamanho n na tela */
void ImprimeMatriz(double matriz[max][max], int n) {
    int i, j;
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            printf(" % lf ", matriz[i][j]);
        }
        printf("\n");
    }
}

/* Transforma a matriz em uma identidade de tamanho n multiplicada por valor (Identidade Especial) */
void IdentEsp(double matriz[max][max], int n, double valor) {
    int i, j;
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            matriz[i][j] = (i == j)*valor;
        }
    }
}

/* resultado = A + B */
void SomaMatriz(double resultado[max][max], double A[max][max], double B[max][max], int n) {
    int i, j;
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            resultado[i][j] = A[i][j] + B[i][j];
        }
    }
}

/* resultado = A - B */
void SubMatriz(double resultado[max][max], double A[max][max], double B[max][max], int n) {
    int i, j;
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            resultado[i][j] = A[i][j] - B[i][j];
        }
    }
}

/* resultado = aux = A x B */
void MultMatriz(double resultado[max][max], double A[max][max], double B[max][max], int n1, int n2, int n3) {
    int i, j, k;
    double aux[max][max];
    for (i = 0; i < n1; i++) {
        for (j = 0; j < n3; j++) {
            aux[i][j] = 0;
            for (k = 0; k < n2; k++) {
                aux[i][j] = aux[i][j] + A[i][k]*B[k][j];
            }
        }
    }
    /* Criei a distincao entre resultado e aux para fazer contas do tipo C = D x C */
    for (i = 0; i < n1; i++) {
        for (j = 0; j < n3; j++) {
            resultado[i][j] = aux[i][j];
        }
    }
}

/* Transpoe a matriz resultado */
void Transpor(double resultado[max][max], int n) {
    int i, j;
    double aux;
    for (i = 1; i < n; i++) {
        for (j = 0; j < i; j++) {
            aux = resultado[i][j];
            resultado[i][j] = resultado[j][i];
            resultado[j][i] = aux;
        }
    }
}

/* Monta a matriz A com n+1 molas */
void MassaMola(double A[max][max], double ks[max+1], double massa, int n) {
    int i, j;
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            A[i][j] = ((ks[i]+ks[i+1])*(i == j) - ks[j]*(i == j-1) - ks[i]*(i-1 == j))/massa;
        }
    }
}

/* Computa o algoritmo QR sem deslocamento espectral e retorna o numero de iteracoes */
int AlgQRsemDesl(double autovalores[max], double autovetores[max][max], double A[max][max], double precisao, int tamanho) {
    int iteracoes, m, i;
    double tau, coss, seno;
    double Aprox[max][max], aux[max][max], Q[max][max], R[max][max];

    iteracoes = 0;
    for (m = tamanho-1; m > 0; m--) {
        while (fabs(A[m][m-1]) > precisao) {

            IdentEsp(aux, tamanho, 0); // Zera a matriz aux
            SubMatriz(Aprox, A, aux, tamanho); // Copia a matriz A (entrada) para Aprox, pois A ainda sera necessaria
            IdentEsp(Q, tamanho, 1); // Faz Q = identidade

            for (i = 0; i < tamanho-1; i++) {
                /* Calcula os valores de seno e cosseno para a rotacao de Givens */
                if (fabs(Aprox[i][i]) > fabs(Aprox[i+1][i])) {
                    tau = -Aprox[i+1][i]/Aprox[i][i];
                    coss = 1/sqrt(1+tau*tau);
                    seno = coss*tau;
                }
                else {
                    tau = -Aprox[i][i]/Aprox[i+1][i];
                    seno = 1/sqrt(1+tau*tau);
                    coss = seno*tau;
                }
                /* Calcula os termos uteis ao proximo calculo de seno e cosseno  */
                Aprox[i+1][i+1] = seno*Aprox[i][i+1] + coss*Aprox[i+1][i+1];
                if (i != tamanho-1) {
                    Aprox[i+1][i+2] = coss*Aprox[i+1][i+2];
                }
                /* Monta a matriz de rotacao de Givens */
                IdentEsp(aux, tamanho, 1);
                aux[ i ][ i ] =  coss;
                aux[i+1][i+1] =  coss;
                aux[ i ][i+1] = -seno;
                aux[i+1][ i ] =  seno;
                /* Executa a rotacao de Givens */
                MultMatriz(Q, aux, Q, tamanho, tamanho, tamanho);
                /* Optei pelo caminho mais longo pois acredito que para matrizes pequenas, ate 32x32, o tempo de execucao em C nao melhoraria tanto assim */
            }

            MultMatriz(R, Q, A, tamanho, tamanho, tamanho); // Calcula a matriz R
            Transpor(Q, tamanho); // Como Q eh ortogonal, transpor eh o mesmo que inverter, e as rotacoes de Givens dao a matriz inversa
            MultMatriz(A, R, Q, tamanho, tamanho, tamanho); // Calcula a nova matriz A
            MultMatriz(autovetores, autovetores, Q, tamanho, tamanho, tamanho); // Atualiza a matriz de autovetores
            iteracoes++;
        }
        /* Zera os elementos que estao menores que precisao */
        A[m][m-1] = 0;
        A[m-1][m] = 0;
    }
    /* Separa os autovalores da matriz A */
    for (i = 0; i < tamanho; i++) {
        autovalores[i] = A[i][i];
    }
    return iteracoes;
}

/* Computa o algoritmo QR com deslocamento espectral e retorna o numero de iteracoes */
int AlgQRcomDesl(double autovalores[max], double autovetores[max][max], double A[max][max], double precisao, int tamanho) {
    int iteracoes, m, i;
    double d, mi, tau, coss, seno;
    double Aprox[max][max], aux[max][max], Q[max][max], R[max][max];

    iteracoes = 0;
    mi = 0;
    for (m = tamanho-1; m > 0; m--) {
        while (fabs(A[m][m-1]) > precisao) {
            /* Calcula o valor de mi pela heuristica de Wilkinson, e separei em funcoes diferentes para nao precisar de um if a mais a cada iteracao */
            if (iteracoes != 0) {
                d = (A[m-1][m-1] - A[m][m])/2;
                mi = A[m][m] + d - Sinal(d)*sqrt(d*d+A[m][m-1]*A[m][m-1]);
            }

            IdentEsp(aux, tamanho, mi); // Cria a matriz que faz o deslocamento
            SubMatriz(Aprox, A, aux, tamanho); // Subtrai o deslocamento da matriz de entrada
            IdentEsp(Q, tamanho, 1); // Prepara a matriz Q

            for (i = 0; i < tamanho-1; i++) {
                /* Calcula os valores de seno e cosseno para a rotacao de Givens */
                if (fabs(Aprox[i][i]) > fabs(Aprox[i+1][i])) {
                    tau = -Aprox[i+1][i]/Aprox[i][i];
                    coss = 1/sqrt(1+tau*tau);
                    seno = coss*tau;
                }
                else {
                    tau = -Aprox[i][i]/Aprox[i+1][i];
                    seno = 1/sqrt(1+tau*tau);
                    coss = seno*tau;
                }
                /* Calcula os termos uteis ao proximo calculo de seno e cosseno  */
                Aprox[i+1][i+1] = seno*Aprox[i][i+1] + coss*Aprox[i+1][i+1];
                if (i != tamanho-1) {
                    Aprox[i+1][i+2] = coss*Aprox[i+1][i+2];
                }
                /* Monta a matriz de rotacao de Givens */
                IdentEsp(aux, tamanho, 1);
                aux[ i ][ i ] =  coss;
                aux[i+1][i+1] =  coss;
                aux[ i ][i+1] = -seno;
                aux[i+1][ i ] =  seno;
                /* Executa a rotacao de Givens */
                MultMatriz(Q, aux, Q, tamanho, tamanho, tamanho);
                /* Optei pelo caminho mais longo pois acredito que para matrizes pequenas, ate 32x32, o tempo de execucao em C nao melhoraria tanto assim */
            }
            IdentEsp(aux, tamanho, mi); // Novamente calcula a matriz responsavel pelo deslocamento
            SubMatriz(Aprox, A, aux, tamanho); // Calcula a matriz que entrou na fatoracao QR
            MultMatriz(R, Q, Aprox, tamanho, tamanho, tamanho); // Calcula a matriz R
            Transpor(Q, tamanho); // Como Q eh ortogonal, transpor eh o mesmo que inverter, e as rotacoes de Givens dao a matriz inversa
            MultMatriz(A, R, Q, tamanho, tamanho, tamanho); // Calcula a nova matriz A
            SomaMatriz(A, A, aux, tamanho); // Soma o deslocamento
            MultMatriz(autovetores, autovetores, Q, tamanho, tamanho, tamanho); // Atualiza os autovetores
            iteracoes++;
        }
        /* Zera os elementos que estao menores que precisao */
        A[m][m-1] = 0;
        A[m-1][m] = 0;
    }
    /* Separa os autovalores da matriz A */
    for (i = 0; i < tamanho; i++) {
        autovalores[i] = A[i][i];
    }
    return iteracoes;
}

/***** Funcao principal ****************************************************************/
int main() {
    /* Definicao das variaveis */
    FILE *fptr;
    time_t start, end;
    int tamanho, i, j, t, iter, quant, teste, desloc, opcaoB, testeB, opcaoC, testeC;
    double precisao, massa, passo, omega, moment;
    double autovalores[max], ks[max];
    double entrada[max][max], autovetores[max][max], aux[max][max], X0[max][max], Y0[max][max];
    /* Usei variaveis int para guardar valores ASCII pois tive problemas com o scanf lendo dois caracteres ao mesmo tempo */

    /* Inicializacao de algumas variaveis */
    iter = 0;
    precisao = 0;
    tamanho = 0;
    teste = 'z';
    desloc = 'z';
    opcaoB = 'z';
    testeB = 'z';
    opcaoC = 'z';
    testeC = 'z';
    passo = 0.025;
    quant = 401;
    massa = 0;

    /* A cada while, faz uma pergunta ate receber uma resposta valida */
    while (!(teste == 'a' || teste == 'b' || teste == 'c')) {
        printf("Qual teste deseja rodar? [a;b;c] ");
        scanf(" %c", &teste);
        fflush(stdin);
    }

    while (!(desloc == 'c' || desloc == 's')) {
        printf("Com ou sem deslocamento espectral? [c;s] ");
        scanf(" %c", &desloc);
        fflush(stdin);
    }

    while (!(0 < precisao && precisao < 1)) {
        printf("Qual a precisao do teste? [0 < decimal com ponto < 1] ");
        scanf(" %lf", &precisao);
        fflush(stdin);
    }
    /* Se o teste escolhido for o A */
    if (teste == 'a') {
        printf("    Realizando o teste A\n");
        while (!(2 <= tamanho && tamanho <= max)) {
            printf("Qual o tamanho da matriz? [2 <= inteiro <= %d] ", max);
            scanf(" %d", &tamanho);
            fflush(stdin);
        }
        /* Cria a matriz de entrada com 2 na diagonal principal e -1 nas diagonais adjacentes */
        for (i = 0; i < tamanho; i++) {
            for (j = 0; j < tamanho; j++) {
                entrada[i][j] = 2*(i == j) - (i == j-1) - (i-1 == j);
            }
        }
        /* Imprime a matriz de entrada e ajeita a matriz de autovetores */
        printf("Matriz de entrada:\n");
        ImprimeMatriz(entrada, tamanho);
        IdentEsp(autovetores, tamanho, 1);
        /* Chama o algoritmo certo, com ou sem deslocamento */
        time(&start);
        if (desloc == 's') {
            iter = AlgQRsemDesl(autovalores, autovetores, entrada, precisao, tamanho);
        }
        else {
            iter = AlgQRcomDesl(autovalores, autovetores, entrada, precisao, tamanho);
        }
        time(&end);
        /* Imprime os autovalores, os autovetores e o numero de iteracoes */
        printf("\nAutovalores:\n");
        for (i = 0; i < tamanho; i++) {
            printf("    Autovalor[%03d] = % lf\n", i+1, autovalores[i]);
        }
        printf("Matriz de autovetores:\n");
        ImprimeMatriz(autovetores, tamanho);
        printf("Numero de iteracoes: %d\n", iter);
        printf("Tempo de execucao: %ld\n", start);
        printf("Tempo de execucao: %ld\n", end);
    }
    /* Se nao era o teste A, pode ser o B */
    else if (teste == 'b') {
        printf("    Realizando o teste B\n");
        while (!(opcaoB == 'v' || opcaoB == 'p')) {
            printf("Deseja entrar com valores(v) ou executar um teste predefinido(p)? [v;p] ");
            scanf(" %c", &opcaoB);
            fflush(stdin);
        }
        /* Se optou por entrar valores predefinidos */
        if (opcaoB == 'p') {
            /* Configura o tamanho, a massa e as constantes elasticas */
            tamanho = 5;
            massa = 2;
            for (i = 0; i < tamanho+1; i++) {
                ks[i] = 40+2*(i+1);
            }
            /* Cria a matriz de entrada e calcula os autovalores e autovetores */
            MassaMola(entrada, ks, massa, tamanho);
            IdentEsp(autovetores, tamanho, 1);
            if (desloc == 's') {
                iter = AlgQRsemDesl(autovalores, autovetores, entrada, precisao, tamanho);
            }
            else {
                iter = AlgQRcomDesl(autovalores, autovetores, entrada, precisao, tamanho);
            }
            while (!(testeB == '1' || testeB == '2' || testeB == '3')) {
                printf("Qual teste predefinido? [1;2;3] ");
                scanf(" %c", &testeB);
                fflush(stdin);
            }
            if (testeB == '1') { // Primeira opcao para X0
                X0[0][0] = -2;
                X0[1][0] = -3;
                X0[2][0] = -1;
                X0[3][0] = -3;
                X0[4][0] = -1;
                fptr = fopen(".\\saidaB1.txt", "w");
            }
            else if (testeB == '2') { // Segunda opcao para X0
                X0[0][0] =  1;
                X0[1][0] = 10;
                X0[2][0] = -4;
                X0[3][0] =  3;
                X0[4][0] = -2;
                fptr = fopen(".\\saidaB2.txt", "w");
            }
            else { // Terceira opcao para X0, os modos de vibracao para maior frequencia
                for (i = 0; i < tamanho; i++) {
                    X0[i][0] = autovetores[i][0];
                }
                fptr = fopen(".\\saidaB3.txt", "w");
            }
            /* Imprime os parametros de entrada e a saida */
            printf("\nTeste com 5 massas de 2kg e 6 molas de constantes elasticas:\n");
            for (i = 0; i < tamanho+1; i++) {
                printf("    k[%03d] = %lf\n", i+1, ks[i]);
            }
            printf("Posicao inicial:\n");
            for (i = 0; i < tamanho; i++) {
                printf("    X0[%03d] = % lf\n", i+1, X0[i][0]);
            }
            printf("\nAutovalores:\n");
            for (i = 0; i < tamanho; i++) {
                printf("    Autovalor[%03d] = % lf\n", i+1, autovalores[i]);
            }
            printf("Matriz de autovetores:\n");
            ImprimeMatriz(autovetores, tamanho);
            printf("Numero de iteracoes: %d\n", iter);
            printf("Os dados de saida estao no arquivo 'saidaB%d.txt'\n", testeB-48);
        }
        else {/* Se opcaoB == 'v' */
            while (!(2 <= tamanho && tamanho <= max)) {
                printf("Quantas massas? [2 <= inteiro <= %d] ", max);
                scanf(" %d", &tamanho);
                fflush(stdin);
            }
            while (!(massa > 0)) {
                printf("Qual a massa? [real > 0] ");
                scanf(" %lf", &massa);
                fflush(stdin);
            }
            for (i = 0; i < tamanho; i++) {
                printf("Qual a posicao inicial da massa[%03d]? [real] ", i+1);
                scanf(" %lf", &X0[i][0]);
                fflush(stdin);
            }
            /* As constantes elasticas foram definidas pelo enunciado do EP1 */
            for (i = 0; i < tamanho+1; i++) {
                ks[i] = 40+2*(i+1);
            }

            MassaMola(entrada, ks, massa, tamanho);
            IdentEsp(autovetores, tamanho, 1);
            if (desloc == 's') {
                iter = AlgQRsemDesl(autovalores, autovetores, entrada, precisao, tamanho);
            }
            else {
                iter = AlgQRcomDesl(autovalores, autovetores, entrada, precisao, tamanho);
            }
            printf("\nTeste com %d massas de %.3lfkg e %d molas de constantes elasticas:\n", tamanho, massa, tamanho+1);
            for (i = 0; i < tamanho+1; i++) {
                printf("    k[%03d] = %lf\n", i+1, ks[i]);
            }
            printf("Posicao inicial:\n");
            for (i = 0; i < tamanho; i++) {
                printf("    X0[%03d] = % lf\n", i+1, X0[i][0]);
            }
            printf("\nAutovalores:\n");
            for (i = 0; i < tamanho; i++) {
                printf("    Autovalor[%03d] = % lf\n", i+1, autovalores[i]);
            }
            printf("Matriz de autovetores:\n");
            ImprimeMatriz(autovetores, tamanho);
            printf("Numero de iteracoes: %d\n", iter);
            printf("Os dados de saida estao no arquivo 'saidaBv.txt'\n");
            fptr = fopen(".\\saidaBv.txt", "w");
        }
        /* Independente de qual for a entrada */
        IdentEsp(aux, tamanho, 1); // Cria aux = identidade
        MultMatriz(aux, autovetores, aux, tamanho, tamanho, tamanho); // copia os autovetores para aux
        Transpor(aux, tamanho); // transpoe os autovetores
        MultMatriz(Y0, aux, X0, tamanho, tamanho, 1); // Calcula os Y0 multiplicando aux por X0
    }
    /* Se nao era A nem B, so pode ser C */
    else if (teste == 'c') {
        /* Realiza as mesmas perguntas que o teste B, mas com 10 massas */
        printf("    Realizando o teste C\n");
            while (!(opcaoC == 'v' || opcaoC == 'p')) {
            printf("Deseja entrar com valores(v) ou executar um teste predefinido(p)? [v;p] ");
            scanf(" %c", &opcaoC);
            fflush(stdin);
        }
        
        if (opcaoC == 'p') {
            /* predefine valores */
            tamanho = 10;
            massa = 2;
            for (i = 0; i < tamanho+1; i++) {
                ks[i] = 40+2*pow(-1,i+1);
            }
            MassaMola(entrada, ks, massa, tamanho);
            IdentEsp(autovetores, tamanho, 1);
            if (desloc == 's') {
                iter = AlgQRsemDesl(autovalores, autovetores, entrada, precisao, tamanho);
            }
            else {
                iter = AlgQRcomDesl(autovalores, autovetores, entrada, precisao, tamanho);
            }

            while (!(testeC == '1' || testeC == '2' || testeC == '3')) {
                printf("Qual teste predefinido? [1;2;3] ");
                scanf(" %c", &testeC);
                fflush(stdin);
            }
            if (testeC == '1') {
                X0[0][0] = -2;
                X0[1][0] = -3;
                X0[2][0] = -1;
                X0[3][0] = -3;
                X0[4][0] = -1;
                X0[5][0] = -2;
                X0[6][0] = -3;
                X0[7][0] = -1;
                X0[8][0] = -3;
                X0[9][0] = -1;
                fptr = fopen(".\\saidaC1.txt", "w");
            }
            else if (testeC == '2') {
                X0[0][0] =  1;
                X0[1][0] = 10;
                X0[2][0] = -4;
                X0[3][0] =  3;
                X0[4][0] = -2;
                X0[5][0] =  1;
                X0[6][0] = 10;
                X0[7][0] = -4;
                X0[8][0] =  3;
                X0[9][0] = -2;
                fptr = fopen(".\\saidaC2.txt", "w");
            }
            else {
                for (i = 0; i < tamanho; i++) {
                    X0[i][0] = autovetores[i][0];
                }
                fptr = fopen(".\\saidaC3.txt", "w");
            }
            
            printf("\nTeste com 10 massas de 2kg e 11 molas de constantes elasticas:\n");
            for (i = 0; i < tamanho+1; i++) {
                printf("    k[%03d] = %lf\n", i+1, ks[i]);
            }
            printf("Posicao inicial:\n");
            for (i = 0; i < tamanho; i++) {
                printf("    X0[%03d] = % lf\n", i+1, X0[i][0]);
            }
            printf("\nAutovalores:\n");
            for (i = 0; i < tamanho; i++) {
                printf("    Autovalor[%03d] = % lf\n", i+1, autovalores[i]);
            }
            printf("Matriz de autovetores:\n");
            ImprimeMatriz(autovetores, tamanho);
            printf("Numero de iteracoes: %d\n", iter);
            printf("Os dados de saida estao no arquivo 'saidaC%d.txt'\n", testeC-48);
        }
        else {/* Se opcaoC == 'v' */
            while (!(2 <= tamanho && tamanho <= max)) {
                printf("Quantas massas? [2 <= inteiro <= %d] ", max);
                scanf(" %d", &tamanho);
                fflush(stdin);
            }
            while (!(massa > 0)) {
                printf("Qual a massa(kg)? [real > 0] ");
                scanf(" %lf", &massa);
                fflush(stdin);
            }
            for (i = 0; i < tamanho; i++) {
                printf("Qual a posicao inicial da massa[%03d]? [real] ", i+1);
                scanf(" %lf", &X0[i][0]);
                fflush(stdin);
            }
            /* Novamente, as constantes elasticas foram definidas pelo enunciado do EP1 */
            for (i = 0; i < tamanho+1; i++) {
                ks[i] = 40+2*pow(-1,i+1);
            }
            MassaMola(entrada, ks, massa, tamanho);
            IdentEsp(autovetores, tamanho, 1);
            if (desloc == 's') {
                iter = AlgQRsemDesl(autovalores, autovetores, entrada, precisao, tamanho);
            }
            else {
                iter = AlgQRcomDesl(autovalores, autovetores, entrada, precisao, tamanho);
            }

            printf("\nTeste com %d massas de %.3lfkg e %d molas de constantes elasticas:\n", tamanho, massa, tamanho+1);
            for (i = 0; i < tamanho+1; i++) {
                printf("    k[%03d] = %lf\n", i+1, ks[i]);
            }
            printf("Posicao inicial:\n");
            for (i = 0; i < tamanho; i++) {
                printf("    X0[%03d] = % lf\n", i+1, X0[i][0]);
            }
            printf("\nAutovalores:\n");
            for (i = 0; i < tamanho; i++) {
                printf("    Autovalor[%03d] = % lf\n", i+1, autovalores[i]);
            }
            printf("Matriz de autovetores:\n");
            ImprimeMatriz(autovetores, tamanho);
            printf("Numero de iteracoes: %d\n", iter);
            printf("Os dados de saida estao no arquivo 'saidaCv.txt'\n");
            fptr = fopen(".\\saidaCv.txt", "w");
        }
        /* Faz os mesmos calculos como se fosse o teste B */
        IdentEsp(aux, tamanho, 1);
        MultMatriz(aux, autovetores, aux, tamanho, tamanho, tamanho);
        Transpor(aux, tamanho);
        MultMatriz(Y0, aux, X0, tamanho, tamanho, 1);
    }
    /* Guarda os valores em funcao do tempo no arquivo txt */
    if (teste == 'b' || teste == 'c') {
        /* Imprime um cabecalho na primeira linha */
        fprintf(fptr, "Tempo");
        for (i = 0; i < tamanho; i++) {
            fprintf(fptr, "  |   x%03d(t) ", i+1);
        }
        for (t = 0; t < quant; t++) { // Para cada valor de tempo ...
            fprintf(fptr, "\n%02.3lf", passo*t);
            for (i = 0; i < tamanho; i++) { // ... calcula a posicao de cada massa ...
                moment = 0;
                fprintf(fptr, "  |  ");
                for (j = 0; j < tamanho; j++) { // ... dependendo da posicao inicial das outras massas
                    omega = sqrt(autovalores[j]);
                    moment = moment + autovetores[i][j]*Y0[j][0]*cos(omega*passo*t);
                }
                fprintf(fptr, "% .6lf", moment);
            }
        }
        /* Fecha o arquivo */
        fclose(fptr);
    }
    /* Espera uma tecla qualquer para nao fechar o programa assim que terminar a execucao e faz um aviso sonoro */
    printf("\n\n\aPressione qualquer tecla para fechar...");
    getchar();
    return 0;
}