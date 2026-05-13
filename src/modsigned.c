#include "modsigned.h"
#include <stddef.h>

size_t modSigned(int value, size_t ArrayLength)
{
    return (size_t)(((value % (int)ArrayLength) + (int)ArrayLength) %
                    (int)ArrayLength);
}
