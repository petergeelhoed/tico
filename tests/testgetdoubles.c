
#include <ctype.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mylib.h"

static void
print_result(const char* label, int n, const double* vals, size_t max_count)
{
    printf("%s\n", label);
    if (n < 0)
    {
        printf("  Result: ERROR (n=%d)\n\n", n);
        return;
    }
    printf("  Parsed count: %d\n", n);
    for (int i = 0; i < n && i < (int)max_count; i++)
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
static void run_auto_test(const char* input, size_t max_count)
{
    double vals[BUF_SIZE];
    if (max_count > BUF_SIZE)
    {
        max_count = BUF_SIZE;
    }
    printf("=== Auto test input: \"%s\" ===\n", input);

    size_t len = strlen(input);
    char* input_mutable = (char*)malloc(len);
    if (input_mutable == NULL)
    {
        exit(EXIT_FAILURE);
    }
    memcpy(input_mutable, input, len);

    // POSIX-ish: fmemopen available on glibc; also present on some platforms.
    FILE* mem = fmemopen((void*)input_mutable, len, "r");
    if (!mem)
    {
        free(input_mutable);
        perror("fmemopen failed");
        return;
    }
    // Temporarily redirect stdin
    FILE* orig_stdin = stdin;
    stdin = mem;

    int count = getDoublesFromStdin(max_count, vals);

    // Restore stdin and close
    stdin = orig_stdin;
    (void)fclose(mem);
    free(input_mutable);

    print_result("Auto test result:", count, vals, max_count);
}

int main(void)
{
    /* --------- Automated tests --------- */
    run_auto_test("3.14 -2.5e3 junk 42\n", BUF_SIZE);
    run_auto_test("a=1.0; b=2.0; c=nan; d=inf\n", BUF_SIZE);
    run_auto_test("no numbers here!\n", BUF_SIZE);
    run_auto_test("   +.5  -0.125  1e-9  2E+10\n", BUF_SIZE);
    run_auto_test("comma, separated: 12.3, -4.56, 78.9\n", BUF_SIZE);

    /* --------- Interactive test --------- */
    printf("Enter a line with doubles (mixed text is okay). Press Ctrl+D "
           "(Unix) or Ctrl+Z (Windows) to end:\n");

    double arr[BUF_SIZE];
    int count = getDoublesFromStdin(BUF_SIZE, arr);
    print_result("Interactive parse result:", count, arr, BUF_SIZE);

    return 0;
}
