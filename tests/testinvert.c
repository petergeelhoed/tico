#include <stdio.h>
#include <stdlib.h>

#include "mymath.h"

#define numElements 9
int main(void)
{

    const unsigned int numRows = 3;
    const unsigned int numCols = 3;
    // NOLINTBEGIN[readability-magic-numbers]
    double xarr[numElements] = {1.0, 2.0, -1.0, 2.0, 1.0, 2.0, -1.0, 2.0, 1.0};
    // NOLINTEND[readability-magic-numbers]

    for (unsigned int col = 0; col < numCols; col++)
    {
        for (unsigned int row = 0; row < numRows; row++)
        {
            printf("%8.1f ", xarr[row + col * numRows]);
        }
        printf("\n");
    }
    printf("\n");
    invert(xarr, numRows, numCols);
    printf("\n");

    for (unsigned int col = 0; col < numCols; col++)
    {
        for (unsigned int row = 0; row < numRows; row++)
        {
            printf("%8.3f ", xarr[row + col * numRows]);
        }
        printf("\n");
    }
    exit(0);
}
