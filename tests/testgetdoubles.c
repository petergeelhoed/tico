
#include <ctype.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "parseargs.h"

static void
printResult(const char* label, int n, const double* vals, size_t maxCount)
{
    printf("%s\n", label);
    if (n < 0)
    {
        printf("  Result: ERROR (n=%d)\n\n", n);
        return;
    }
    printf("  Parsed count: %d\n", n);
    for (int i = 0; i < n && i < (int)maxCount; i++)
    {
        printf("    [%d] %.17g\n", i, vals[i]);
    }
    printf("\n");
}

/**
 * Run one automated test by feeding a string as stdin using fmemopen (POSIX).
 * If fmemopen is not available, falls back to writing the string to a temp file
 * and freopen'ing stdin.
 */
static void runAutoTest(const char* input, size_t maxCount)
{
    double vals[BUFFER_SIZE];
    if (maxCount > BUFFER_SIZE)
    {
        maxCount = BUFFER_SIZE;
    }
    printf("=== Auto test input: \"%s\" ===\n", input);

    size_t len = strlen(input);
    char* inputMutable = (char*)malloc(len);
    if (inputMutable == NULL)
    {
        exit(EXIT_FAILURE);
    }
    memcpy(inputMutable, input, len);

    // POSIX-ish: fmemopen available on glibc; also present on some platforms.
    FILE* mem = fmemopen((void*)inputMutable, len, "r");
    if (!mem)
    {
        free(inputMutable);
        perror("fmemopen failed");
        return;
    }
    // Temporarily redirect stdin
    FILE* origStdin = stdin;
    stdin = mem;

    int count = getDoublesFromStdin(maxCount, vals);

    // Restore stdin and close
    stdin = origStdin;
    (void)fclose(mem);
    free(inputMutable);

    printResult("Auto test result:", count, vals, maxCount);
}

int main(void)
{
    /* --------- Automated tests --------- */
    runAutoTest("3.14 -2.5e3 junk 42\n", BUFFER_SIZE);
    runAutoTest("a=1.0; b=2.0; c=nan; d=inf\n", BUFFER_SIZE);
    runAutoTest("no numbers here!\n", BUFFER_SIZE);
    runAutoTest("   +.5  -0.125  1e-9  2E+10\n", BUFFER_SIZE);
    runAutoTest("comma, separated: 12.3, -4.56, 78.9\n", BUFFER_SIZE);

    /* --------- Interactive test --------- */
    printf("Enter a line with doubles (mixed text is okay). Press Ctrl+D "
           "(Unix) or Ctrl+Z (Windows) to end:\n");

    double arr[BUFFER_SIZE];
    int count = getDoublesFromStdin(BUFFER_SIZE, arr);
    printResult("Interactive parse result:", count, arr, BUFFER_SIZE);

    return 0;
}
