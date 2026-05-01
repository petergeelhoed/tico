#define MATMUL_M 3
#define MATMUL_N 3
#define MATMUL_T 2
#define MATMUL_S 3
#include <stdio.h>
#include <stdlib.h>

#include "mymath.h"

int main(void)
{
    unsigned int matM = MATMUL_M;
    unsigned int matN = MATMUL_N;
    // NOLINTBEGIN[readability-magic-numbers]
    double xarr[MATMUL_M * MATMUL_N] =
        {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0};

    unsigned int matT = MATMUL_T;
    unsigned int matS = MATMUL_S;
    double yvec[MATMUL_S * MATMUL_T] = {1.0, 2.0, 1.0, 2.0, 1.0, 2.0};
    // NOLINTEND[readability-magic-numbers]

    for (unsigned int j = 0; j < matN; j++)
    {
        for (unsigned int i = 0; i < matM; i++)
        {
            printf("%8.1f ", xarr[i + j * matM]);
        }
        printf("\n");
    }
    printf("\n");

    for (unsigned int j = 0; j < matS; j++)
    {
        for (unsigned int i = 0; i < matT; i++)
        {
            printf("%8.1f ", yvec[i + j * matT]);
        }
        printf("\n");
    }

    printf("\n");
    double* tmp = mulmat(xarr, matN, matM, yvec, matS, matT);
    printf("\n");

    for (unsigned int j = 0; j < matN; j++)
    {
        for (unsigned int i = 0; i < matT; i++)
        {
            printf("%8.3f ", tmp[i + j * matT]);
        }
        printf("\n");
    }
    exit(0);
}
