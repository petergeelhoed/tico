#include "mylib.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define BUF_SIZE 32

int main(void)
{
    // Part 1: Basic tests
    char test[BUF_SIZE] = " 12 ";
    printf("'%s': %d\n", test, getInt(test));
    printf("'%s': %lf\n", test, getDouble(test));

    char testd[BUF_SIZE] = "12.4";
    // Note: getInt("12.4") will return 12 and stop at the decimal point
    printf("'%s': %d\n", testd, getInt(testd));
    printf("'%s': %lf\n", testd, getDouble(testd));

    // Part 2: Working with a file stream
    char testl[] = "1.0\n3.1\n \n42.7\n 1";

    // Open a file handle to the string testl for reading
    FILE* fpInput = fmemopen(testl, strlen(testl), "r");
    if (fpInput == NULL)
    {
        perror("fmemopen failed");
        return EXIT_FAILURE;
    }

    char* line = NULL;
    size_t len = 0;

    printf("\n--- Reading from simulated file ---\n");
    while (getline(&line, &len, fpInput) != -1)
    {
        // Skip lines that are just whitespace or empty to avoid getDouble
        // exiting
        /*        if (line[0] == '\n' || line[0] == '\r' || line[0] == ' ')
                {
                    continue;
                }
        */
        printf("Line Double Read: '%s' -> Value: %lf\n", line, getDouble(line));
        printf("Line Int    Read: '%s' -> Value: %d\n", line, getInt(line));
    }

    // Cleanup
    free(line);
    (void)fclose(fpInput);

    return 0;
}

change the code to check for the correct doubles and ints
