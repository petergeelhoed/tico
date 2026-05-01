#include <stdio.h>
#include <stdlib.h>

#include "mymath.h"



int main(void)
{
    unsigned int matM = 3;
    unsigned int matN = 3;
    double xarr[9] = // NOLINTBEGIN[readability-magic-numbers]
        {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0}; // NOLINTEND[readability-magic-numbers]

    unsigned int matT = 2;
    unsigned int matS = 3;
    double yvec[6] = // NOLINTBEGIN[readability-magic-numbers]
        {1.0, 2.0, 1.0, 2.0, 1.0, 2.0}; // NOLINTEND[readability-magic-numbers]

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
