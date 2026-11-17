/* stub POSIX regex pour mingw : sfscan n'utilise pas les expressions de filtre htslib */
#ifndef SF_STUB_REGEX_H
#define SF_STUB_REGEX_H
#include <stddef.h>
typedef struct { size_t re_nsub; void *p; } regex_t;
typedef long regoff_t;
typedef struct { regoff_t rm_so, rm_eo; } regmatch_t;
#define REG_EXTENDED 1
#define REG_ICASE 2
#define REG_NOSUB 4
#define REG_NEWLINE 8
#define REG_NOMATCH 1
static inline int regcomp(regex_t *r, const char *p, int f){(void)r;(void)p;(void)f;return 2;}
static inline int regexec(const regex_t *r, const char *s, size_t n, regmatch_t m[], int f){(void)r;(void)s;(void)n;(void)m;(void)f;return REG_NOMATCH;}
static inline void regfree(regex_t *r){(void)r;}
static inline size_t regerror(int e, const regex_t *r, char *b, size_t n){(void)e;(void)r;if(n){b[0]=0;}return 0;}
#endif
