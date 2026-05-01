#include <stdio.h>
#include <stdlib.h>

#include "mymath.h"

int main(void)
{


    unsigned int numRows = 3;
    unsigned int numCols = 3;
    double xarr[9] = // NOLINTBEGIN[readability-magic-numbers]
        {1.0, 2.0, -1.0, 2.0, 1.0, 2.0, -1.0, 2.0, 1.0}; // NOLINTEND[readability-magic-numbers]

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
