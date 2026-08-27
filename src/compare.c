#include "compare.h"

int compare_double(const void* first, const void* second)
{
    double dfirst = *(const double*)first;
    double dsecond = *(const double*)second;

    return (dfirst > dsecond) - (dfirst < dsecond);
}
