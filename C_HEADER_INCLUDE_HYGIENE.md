# C/C++ Header Include Hygiene

- The header file with the same name as the source file should always be the first include in each .c file.
- Project headers should be grouped together in a block at the top, followed by a separate block for standard library includes below.
- Only include headers in .h files if their types, macros, or function prototypes are directly used in the header itself.
- Do not include headers in .h files just because they are needed in the .c implementation—include them in the .c file instead.
- Always preserve the original formatting, order, and whitespace of includes unless explicitly requested to change them.
- When refactoring, check for and remove unnecessary includes, but do not add or remove blank lines or reorder includes unless necessary for clarity or style consistency.
- If a type like size_t is used in a header, include <stddef.h> in the header. If malloc/free/exit are only used in the .c file, include <stdlib.h> only in the .c file.
- Never introduce formatting changes (such as moving #pragma once or adding/removing blank lines) unless explicitly requested.
- When in doubt, match the style and order of includes as found in git HEAD.
