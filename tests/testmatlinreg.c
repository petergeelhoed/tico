#include <stdio.h>
#include <stdlib.h>

#include "mymath.h"

int main(void)
{

    unsigned int matM = 1;
    unsigned int matN = 4;
    // NOLINTBEGIN[readability-magic-numbers]
    double xarr[4] = {1.0, 2.0, 3.0, 4.0};

    unsigned int matT = 1;
    unsigned int matS = 4;
    double yvec[4] = {3.0, 5.0, 8.0, 9.0};
    double wvec[4] = {1.0, 1.0, 2.0, 1.0};
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
            printf("%8.1fw%3.1f ", yvec[i + j * matT], wvec[i + j * matT]);
        }
        printf("\n");
    }

    printf("\n");
    double tmp[2];
    matlinreg(tmp, xarr, matN, matM, yvec, wvec);
    printf("\n");

    for (unsigned int i = 0; i < matT + 1; i++)
    {
        printf("%8.3f ", tmp[i]);
    }
    printf("\n");

    fastlinreg(tmp, xarr, matN, yvec, wvec);
    for (unsigned int i = 0; i < matT + 1; i++)
    {
        printf("%8.3f ", tmp[i]);
    }
    printf("\n");

    exit(0);
}
