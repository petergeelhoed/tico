#include "myfft.h"

void crosscorint(unsigned int NN, int* array, int* ref, int* cross)
{

    fftw_complex* tmparr = fftw_alloc_complex(NN);
    fftw_complex* tmpref = fftw_alloc_complex(NN);
    for (unsigned int j = 0; j < NN; j++)
    {
        tmparr[j][0] = array[j];
        tmparr[j][1] = 0.0;
        tmpref[j][0] = ref[j];
        tmpref[j][1] = 0.0;
    }
    fftw_complex* coor = crosscor(NN, tmparr, tmpref);

    for (unsigned int j = 0; j < NN; j++)
    {
        cross[j] = (int)(coor[j][0] * NN);
    }
    fftw_free(tmpref);
    fftw_free(tmparr);
    fftw_free(coor);
}
