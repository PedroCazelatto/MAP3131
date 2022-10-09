/***************************************************************/
/**                                                           **/
/**    Pedro Henrique Galhardi Cazelatto - n° USP 11261090    **/
/**    Exercício-Programa 02 - MAP3121                        **/
/**    Autovalores e Autovetores de Matrizes                  **/
/**    Reais Simetricas - O Algoritmo QR                      **/
/**                                                           **/
/***************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* Define o valor máximo das matrizes do codigo */
#define max 256
#define DtoR 3.141592654/180

/* Calcula o sinal de um número, 1 para positivos e -1 para negativos */
int Sinal(double valor) {
    return 1-2*(valor < 0);
}

/* Imprime uma matriz quadrada de tamanho n na tela */
void ImprimeMatriz(double matriz[max][max], int n) {
    int i, j;
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            printf(" % .4f ", matriz[i][j]);
        }
        printf("\n");
    }
}

/* Transforma a matriz em uma identidade de tamanho n */
void Ident(double matriz[max][max], int n) {
    int i, j;
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            matriz[i][j] = (i == j);
        }
    }
}

/* Soma um valor a todos os n elementos de um vetor */
void SomaVetor(double A[max], double valor, int n) {
    int i;
    for (i = 0; i < n; i++) {
        A[i] = A[i] + valor;
    }
}

/* Retorna a norma de um vetor de tamanho n */
double Norma(double vetor[max], int n) {
    int i;
    double norma2;
    norma2 = 0;
    for (i = 0; i < n; i++) {
        norma2 = norma2 + vetor[i]*vetor[i];
    }
    return sqrt(norma2);
}

/* Retorna o produto escalar entre dois vetores de tamanho n */
double ProdEsc(double vetorA[max], double vetorB[max], int n) {
    int i;
    double prod;
    prod = 0;
    for (i = 0; i < n; i++) {
        prod = prod + vetorA[i]*vetorB[i];
    }
    return prod;
}

/* Retorna o maior valor (equivalente a max()) */
int Maior(int A, int B) {
    if (A > B) {
        return A;
    }
    return B;
}

/* Computa o algoritmo QR com deslocamento espectral e retorna o numero de iteracoes */
int AlgQRcomDesl(double autovetores[max][max], double diagInf[max-1], double diagPri[max], double diagSup[max-1], double precisao, int tamanho) {
    int iteracoes, m, i, j;
    double d, mi, tau, diagSupExtra, Vextra;
    double coss[max], seno[max];

    iteracoes = 0;
    mi = 0;

    for (m = tamanho-1; m > 0; m--) {
        while (fabs(diagInf[m-1]) > precisao) {
            /* Calcula o valor de mi pela heuristica de Wilkinson */
            if (iteracoes != 0) {
                d = (diagPri[m-1] - diagPri[m])/2;
                mi = diagPri[m] + d - Sinal(d)*sqrt(d*d + diagInf[m-1]*diagInf[m-1]);
            }

            SomaVetor(diagPri, -mi, tamanho);

            for (i = 0; i < tamanho-1; i++) {
                /* Calcula os valores de seno e cosseno para a rotacao de Givens */
                if (fabs(diagPri[i]) > fabs(diagInf[i])) {
                    tau = -diagInf[i]/diagPri[i];
                    coss[i] = 1/sqrt(1+tau*tau);
                    seno[i] = coss[i]*tau;
                }
                else {
                    tau = -diagPri[i]/diagInf[i];
                    seno[i] = 1/sqrt(1+tau*tau);
                    coss[i] = seno[i]*tau;
                }

                diagPri[i  ] = coss[i]*diagPri[i  ] - seno[i]*diagInf[i  ];
                diagSupExtra = diagSup[i];
                diagSup[i  ] = coss[i]*diagSupExtra - seno[i]*diagPri[i+1];
                diagPri[i+1] = coss[i]*diagPri[i+1] + seno[i]*diagSupExtra;
                if (i != tamanho-2) {
                    diagSup[i+1] = coss[i]*diagSup[i+1];
                }
            }

            for (i = 0; i < tamanho-1; i++) {
                diagInf[i  ] = -seno[i]*diagPri[i+1];
                diagPri[i  ] =  coss[i]*diagPri[i  ] - seno[i]*diagSup[i  ];
                diagPri[i+1] =  coss[i]*diagPri[i+1];
                diagSup[i  ] = diagInf[i];
            }

            SomaVetor(diagPri, mi, tamanho);

            for (j = 0; j < tamanho-1; j++) {
                for (i = 0; i < tamanho; i++) {
                    Vextra  = autovetores[i][j];
                    autovetores[i][j  ] = coss[j]*Vextra - seno[j]*autovetores[i][j+1];
                    autovetores[i][j+1] = seno[j]*Vextra + coss[j]*autovetores[i][j+1];
                }
            }
            iteracoes++;
        }
        /* Zera os elementos que estao menores que precisao */
        diagInf[m-1] = 0;
        diagSup[m-1] = 0;
    }
    /* Os autovalores retorna em diagPri */
    return iteracoes;
}

void Tridiagonalizacao(double entrada[max][max], double Ht[max][max], int tamanho) {
    int i, j, k;
    double ww, wh, wx, alfah, alfax;
    double w[max-1], h[max-1], x[max-1];

    for (i = 0; i < tamanho-2; i++) {
        for (j = 1; j < tamanho-i; j++) {
            w[j-1] = entrada[j+i][i];
        }
        w[0] = w[0] + Norma(w, tamanho-1-i)*Sinal(w[0]);
        ww = ProdEsc(w, w, tamanho-1-i);

        for (j = i; j < tamanho; j++) {
            for (k = 1; k < tamanho-i; k++) {
                x[k-1] = entrada[k+i][j];
            }
            wx = ProdEsc(w, x, tamanho-1-i);
            alfax = -2*wx/ww;
            for (k = 1; k < tamanho-i; k++) {
                entrada[k+i][j] = entrada[k+i][j] + alfax*w[k-1];
            }
        }

        for (j = i; j < tamanho; j++) {
            for (k = 1; k < tamanho-i; k++) {
                x[k-1] = entrada[j][k+i];
            }
            wx = ProdEsc(w, x, tamanho-1-i);
            alfax = -2*wx/ww;
            for (k = 1; k < tamanho-i; k++) {
                entrada[j][k+i] = entrada[j][k+i] + alfax*w[k-1];
            }
        }

        for (j = 1; j < tamanho; j++) {
            for (k = 1; k < tamanho-i; k++) {
                h[k-1] = Ht[j][k+i];
            }
            wh = ProdEsc(w, h, tamanho-1-i);
            alfah = -2*wh/ww;
            for (k = 1; k < tamanho-i; k++) {
                Ht[j][k+i] = Ht[j][k+i] + alfah*w[k-1];
            }
        }
    }
}

/***** Funcao principal ****************************************************************/
int main() {
    /* Definicao das variaveis */
    FILE *file;
    int tamanho, i, j, k, iter, totalNos, nosLivres, barras;
    int posicao[5];
    char teste, opcao;
    char arquivo[max];
    double precisao, rho, area, elast, numer, compr, angle, coss, seno, mult, rhoAL, somar, anterior;
    double diagInf[max], diagPri[max], diagSup[max], massas[max], minFreq[5];
    double entrada[max][max], autovetores[max][max], minModo[max][5];

    /* Inicializacao de algumas variaveis */
    iter      = 0;
    tamanho   = 0;
    totalNos  = 0;
    nosLivres = 0;
    barras    = 0;
    numer     = 0;
    rho       = 0;
    area      = 0;
    elast     = 0;
    anterior  = 0;
    precisao  = 1e-6;
    teste     = 'z';
    opcao     = 'z';

    /* A cada while, faz uma pergunta ate receber uma resposta valida */
    while (!(teste == 'm' || teste == 't')) {
        printf("Testar com matriz(m) ou com a trelica(t)? [m;t] ");
        scanf(" %1s", &teste);
    }

    /* Se for testar com matrizes */
    if (teste == 'm') {
        printf("\nRealizando o teste com matrizes\n");
        while (!(opcao == 'v' || opcao == 'a')) {
            printf("Deseja digitar valores(v) ou ler de um arquivo(a)? [v;a] ");
            scanf(" %1s", &opcao);
        }
        /* Se for ler de arquivo */
        if (opcao == 'a') {
            /* Pega o nome do arquivo e tenta abrir para leitura */
            printf("Digite o nome do arquivo: ");
            scanf("%s", &arquivo);
            file = fopen(arquivo, "r");
            /* Se o arquivo nao existe, retorna erro */
            if (file == NULL) {
                printf("Erro: Arquivo nao existente\n");
                fflush(stdin);
                printf("\n\n\aPressione qualquer tecla para fechar...");
                getchar();
                return 0;
            }
            /* Se o arquivo existe, realiza a leitura */
            fscanf(file, " %d", &tamanho);
            for (i = 0; i < tamanho && !feof(file); i++) {
                for (j = 0; j < tamanho && !feof(file); j++) {
                    fscanf(file, " %lf", &entrada[i][j]);
                }
            }
            /* Se parou de ler antes de ler tamanho^2 elementos, avisa que a matriz esta incompleta */
            if (i != tamanho || j != tamanho) {
                printf("Erro: Arquivo nao contem uma matriz %dx%d completa\n", tamanho, tamanho);
                fflush(stdin);
                printf("\n\n\aPressione qualquer tecla para fechar...");
                getchar();
                return 0;
            }
            fclose(file);
        }
        /* Se for entrar com valores */
        else {
            while (!(2 <= tamanho && tamanho <= max)) {
                printf("Qual o tamanho da matriz? [2 <= inteiro <= %d] ", max);
                scanf(" %d", &tamanho);
                fflush(stdin);
            }
            for (i = 0; i < tamanho; i++) {
                for (j = 0; j < tamanho; j++) {
                    printf("Elemento[%03d][%03d] = ", i, j);
                    scanf(" %lf", &entrada[i][j]);
                }
            }
        }

        /* Imprime a matriz de entrada e limpa a matriz Ht */
        /* A matriz Ht eh a matriz autovetores */
        printf("Matriz de entrada:\n");
        ImprimeMatriz(entrada, tamanho);
        Ident(autovetores, tamanho);

        /* Tridiagonaliza a matriz de entrada ajeitando a matriz Ht */
        Tridiagonalizacao(entrada, autovetores, tamanho);

        /* Separa as diagonais em vetores para agilizar os calculos */
        for (i = 0; i < tamanho-1; i++) {
            diagInf[i] = entrada[i+1][i  ];
            diagPri[i] = entrada[i  ][i  ];
            diagSup[i] = entrada[i  ][i+1];
        }
        diagPri[tamanho-1] = entrada[tamanho-1][tamanho-1];

        /* Chama o algoritmo QR com deslocamento */
        iter = AlgQRcomDesl(autovetores, diagInf, diagPri, diagSup, precisao, tamanho);

        /* Imprime os autovalores, os autovetores e o numero de iteracoes */
        printf("\nAutovalores:\n");
        for (i = 0; i < tamanho; i++) {
            printf("    Autovalor[%03d] = % .4f\n", i+1, diagPri[i]);
        }
        printf("Matriz de autovetores:\n");
        ImprimeMatriz(autovetores, tamanho);
        printf("Numero de iteracoes: %d\n", iter);
    }
    /* Se for testar com trelica */
    else {
        printf("\nRealizando o teste com a trelica\n");
        while (!(opcao == 'v' || opcao == 'a')) {
            printf("Deseja digitar valores(v) ou ler de um arquivo(a)? [v;a] ");
            scanf(" %1s", &opcao);
        }
        /* Se for ler de arquivo */
        if (opcao == 'a') {
            /* Pega o nome do arquivo e tenta abrir para leitura */
            printf("Digite o nome do arquivo: ");
            scanf(" %s", &arquivo);
            file = fopen(arquivo, "r");
            /* Se o arquivo nao existe, retorna erro */
            if (file == NULL) {
                printf("Erro: Arquivo nao existente\n");
                fflush(stdin);
                printf("\n\n\aPressione qualquer tecla para fechar...");
                getchar();
                return 0;
            }
            /* Se o arquivo existe, comeca a leitura */

            /* Valores de configuracao */
            fscanf(file, " %d", &totalNos);
            fscanf(file, " %d", &nosLivres);
            fscanf(file, " %d", &barras);
            fscanf(file, " %lf", &rho);
            fscanf(file, " %lf", &area);
            fscanf(file, " %lf", &elast);
            elast = elast*1E9;
            tamanho = 2*nosLivres;
            numer   = area*elast;
            
            /* Limpa a matriz K e o vetor de massas */
            for (i = 0; i < tamanho; i++) {
                for (j = 0; j < tamanho; j++) {
                    entrada[i][j] = 0;
                }
                massas[i] = 0;
            }

            /* Faz a leitura das ligacoes de barras do arquivo */
            for (k = 0; k < barras && !feof(file); k++) {
                fscanf(file,  " %d", &i);
                fscanf(file,  " %d", &j);
                fscanf(file, " %lf", &angle);
                fscanf(file, " %lf", &compr);

                posicao[0] = 2*i-1 -1;
                posicao[1] = 2*i   -1;
                posicao[2] = 2*j-1 -1;
                posicao[3] = 2*j   -1;

                mult = numer/compr;
                rhoAL = 0.5*rho*area*compr;
                coss = cos(angle*DtoR);
                seno = sin(angle*DtoR);

                for (i = 0; i < 4; i++) {
                    if (posicao[i] <= tamanho) {
                        massas[posicao[i]] = massas[posicao[i]] + rhoAL;
                    }
                    for (j = 0; j < 4; j++) {
                        somar = 1.0;
                        somar = somar*pow(coss, (i+1)%2)*pow(coss, (j+1)%2);
                        somar = somar*pow(seno, i%2)*pow(seno, j%2);
                        somar = somar*Sinal(i-1.5)*Sinal(j-1.5);
                        somar = somar*mult;

                        if (posicao[i] <= tamanho && posicao[j] <= tamanho) {
                            entrada[posicao[i]][posicao[j]] = entrada[posicao[i]][posicao[j]] + somar;
                        }
                    }
                }
            }
            fclose(file);
        }
        /* Se for entrar com valores */
        else {
            while (!(2 <= totalNos && totalNos <= max)) {
                printf("Qual o total de nos? [2 <= inteiro <= %d] ", max);
                scanf(" %d", &totalNos);
                fflush(stdin);
            }
            while (!(2 <= nosLivres && nosLivres <= totalNos)) {
                printf("Quantos nos estao livres? [2 <= inteiro <= %d] ", totalNos);
                scanf(" %d", &nosLivres);
                fflush(stdin);
            }
            while (!(2 <= barras && barras <= max)) {
                printf("Qual o total de barras? [2 <= inteiro <= %d] ", max);
                scanf(" %d", &barras);
                fflush(stdin);
            }
            while (!(0 < rho)) {
                printf("Qual a densidade volumetrica? [positivo] ");
                scanf(" %lf", &rho);
                fflush(stdin);
            }
            while (!(0 < area)) {
                printf("Qual a area da secao transversal das barras? [positivo] ");
                scanf(" %lf", &area);
                fflush(stdin);
            }
            while (!(0 < elast)) {
                printf("Qual o modulo de elasticidade das barras (em GPa)? [positivo] ");
                scanf(" %lf", &elast);
                fflush(stdin);
            }
            elast = elast*1E9;
            tamanho = 2*nosLivres;
            numer   = area*elast;
            
            /* Limpa a matriz K e o vetor de massas */
            for (i = 0; i < tamanho; i++) {
                for (j = 0; j < tamanho; j++) {
                    entrada[i][j] = 0;
                }
                massas[i] = 0;
            }

            /* Faz a leitura das ligacoes de barras do arquivo */
            i = 0;
            j = 0;
            for (k = 0; k < barras && !feof(file); k++) {
                while (!(1 <= i && i <= totalNos)) {
                    printf("Primeiro no: [1 <= inteiro <= %d] ", totalNos);
                    scanf(" %d", &i);
                    fflush(stdin);
                }
                while (!(1 <= j && j <= totalNos)) {
                    printf("Segundo no: [1 <= inteiro <= %d] ", totalNos);
                    scanf(" %d", &j);
                    fflush(stdin);
                }
                while (!(0 <= angle && angle <= 180)) {
                    printf("Angulo da barra com a horizontal (graus) [0 <= angulo <= 180] ");
                    scanf(" %lf", &angle);
                    fflush(stdin);
                }
                while (!(0 < compr)) {
                    printf("Qual o comprimento da barra? [positivo] ");
                    scanf(" %lf", &compr);
                    fflush(stdin);
                }

                posicao[0] = 2*i-1 -1;
                posicao[1] = 2*i   -1;
                posicao[2] = 2*j-1 -1;
                posicao[3] = 2*j   -1;

                mult = numer/compr;
                rhoAL = 0.5*rho*area*compr;
                coss = cos(angle*DtoR);
                seno = sin(angle*DtoR);

                for (i = 0; i < 4; i++) {
                    if (posicao[i] <= tamanho) {
                        massas[posicao[i]] = massas[posicao[i]] + rhoAL;
                    }
                    for (j = 0; j < 4; j++) {
                        somar = 1.0;
                        somar = somar*pow(coss, (i+1)%2)*pow(coss, (j+1)%2);
                        somar = somar*pow(seno, i%2)*pow(seno, j%2);
                        somar = somar*Sinal(i-1.5)*Sinal(j-1.5);
                        somar = somar*mult;

                        if (posicao[i] <= tamanho && posicao[j] <= tamanho) {
                            entrada[posicao[i]][posicao[j]] = entrada[posicao[i]][posicao[j]] + somar;
                        }
                    }
                }
            }
        }

        /* Calcula M^-1/2 */
        for (i = 0; i < tamanho; i++) {
            massas[i] = 1/sqrt(massas[i]);
        }

        /* Multiplica a matriz K pela massa */
        for (i = 0; i < tamanho; i++) {
            for (j = 0; j < tamanho; j++) {
                entrada[i][j] = massas[i]*entrada[i][j]*massas[j];
            }
        }
        /* Imprime a matriz de entrada e limpa a matriz Ht */
        /* A matriz Ht eh a matriz autovetores */
        printf("Matriz de entrada:\n");
        ImprimeMatriz(entrada, tamanho);
        Ident(autovetores, tamanho);

        /* Tridiagonaliza a matriz de entrada ajeitando a matriz Ht */
        Tridiagonalizacao(entrada, autovetores, tamanho);

        /* Separa as diagonais em vetores para agilizar os calculos */
        for (i = 0; i < tamanho-1; i++) {
            diagInf[i] = entrada[i+1][i  ];
            diagPri[i] = entrada[i  ][i  ];
            diagSup[i] = entrada[i  ][i+1];
        }
        diagPri[tamanho-1] = entrada[tamanho-1][tamanho-1];

        /* Chama o algoritmo QR com deslocamento */
        iter = AlgQRcomDesl(autovetores, diagInf, diagPri, diagSup, precisao, tamanho);

        /* Imprime os autovalores, os autovetores e o numero de iteracoes */
        printf("\nAutovalores:\n");
        for (i = 0; i < tamanho; i++) {
            printf("    Autovalor[%03d] = % .4f\n", i+1, diagPri[i]);
        }
        printf("Matriz de autovetores:\n");
        ImprimeMatriz(autovetores, tamanho);
        printf("Numero de iteracoes: %d\n", iter);

        /* Separa os 5 menores autovalores */
        for (i = 0; i < 5; i++) {
            minFreq[i] = diagPri[1];
            posicao[i] = 1;
            for (j = 1; j < tamanho; j++) {
                if (diagPri[j] < minFreq[i] && anterior < diagPri[j]) {
                    minFreq[i] = diagPri[j];
                    posicao[i] = j;
                }
            }
            anterior = minFreq[i];
        }
        /* Tira a raiz para encontrar frequencias */
        for (i = 0; i < 5; i++) {
            minFreq[i] = sqrt(minFreq[i]);
        }
        
        /* Agrupa os modos das 5 menores frequencias */
        for (i = 0; i < 5; i++) {
            for (j = 0; j < tamanho; j++) {
                minModo[j][i] = massas[j]*autovetores[j][posicao[i]];
            }
        }

        printf("As 5 menores frequencias e seus modos de vibracao estao salvos no arquivo \"saidaTesteTrelica.txt\"");

        file = fopen(".\\saidaTesteTrelica.txt", "w");
        fprintf(file, " Frequencias (rad/s):");
        for (i = 0; i < 5; i++) {
            fprintf(file, "  %08.3lf |", minFreq[i]);
        }
        fprintf(file, "\n\n Modos de vibracao:  ");
        for (i = 0; i < tamanho; i++) {
            for (j = 0; j < 5; j++) {
                fprintf(file, " % .6f |", minModo[i][j]);
            }
            fprintf(file, "\n                     ");
        }
        fclose(file);
    }
    
    /* Espera uma tecla qualquer para nao fechar o programa assim que terminar a execucao e faz um aviso sonoro */
    fflush(stdin);
    printf("\n\n\aPressione qualquer tecla para fechar...");
    getchar();
    return 0;
}