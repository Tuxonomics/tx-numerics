// tx_tests.hpp

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <setjmp.h>
#include <stdarg.h>
#include <tgmath.h>

#define TX_CONCAT_0(x, y) x ## y
#define TX_CONCAT(x, y) TX_CONCAT_0(x, y)


jmp_buf tx_text_buf;
#define TX_ERROR_BUFFER_SIZE 1024
char tx_error_buffer[TX_ERROR_BUFFER_SIZE];

#define TX_APPROX(a, b, eps) (fabs((a)-(b)) < (eps))

#define tx_approx_definition(T) \
    bool tx_approx_ ## T( T a, T b, T eps ) \
    { \
        if ( b == T(0.0) ) { \
            return TX_APPROX( a, b, eps ); \
        } \
        else { \
            return TX_APPROX( a/b, T(1.0), eps ); \
        } \
    }

tx_approx_definition(float)
tx_approx_definition(double)

#define tx_approx_array_definition(T) \
    bool tx_approx_array_ ## T( T *a, T *b, T eps, size_t n ) \
    { \
        for ( size_t i = 0; i < n; ++i ) { \
            if ( !tx_approx_ ## T(a[i], b[i], eps) ) { \
                return false; \
            } \
        } \
        return true; \
    }

tx_approx_array_definition(float)
tx_approx_array_definition(double)

#define TX_APPROX_ARRAY(a, b, eps, n, type) TX_CONCAT(tx_approx_array_,type)( a, b, eps, n)


void tx_assert_msg(
    int cond, const char *cond_str, const char *func, const char *file, const char *line, const char *msg, ...
) {
    if ( !cond ) {
        if ( msg ) {
            va_list args;
            va_start(args, msg);
            vprintf(msg, args);
            va_end(args);
        }
        snprintf(tx_error_buffer, TX_ERROR_BUFFER_SIZE, "(%s) in %s, file %s, line %s", cond_str, func, file, line);
        longjmp(tx_text_buf, 1);
    }

}


#define TX_AS_STRING2(x) #x
#define TX_AS_STRING(x) TX_AS_STRING2(x)
#define TX_LINE TX_AS_STRING(__LINE__)

#define TX_ASSERT_MSG_VA(cond, msg, ...) tx_assert_msg(cond, #cond, __func__, __FILE__, TX_LINE, msg, __VA_ARGS__)
#define TX_ASSERT_MSG(cond, msg) TX_ASSERT_MSG_VA(cond, msg, NULL)
#define TX_ASSERT(cond) TX_ASSERT_MSG(cond, NULL)

#define TX_TEST_FUN(name) void (*name) (void)


struct tx_test {
    TX_TEST_FUN(f);
    const char *name;
};

#define TX_N_TESTS 1000
#define TX_ADD_TEST(func) (tx_test) { .f = func, .name = #func }
#define TX_TEST_LIST tx_test tx_tests[TX_N_TESTS]

#define TX_TEST_OUTPUT_LEN 70
#define TX_GREEN(text) "\033[32m" text "\033[0m"
#define TX_RED(text)   "\033[31m" text "\033[0m"


#define TX_TESTS_MAIN int main( int argn, const char *argv[] ) \
    {\
        int success = 0;\
\
        for ( int i=0; i<TX_N_TESTS; i++ ) {\
            if ( tx_tests[i].f ) {\
                int failed = setjmp( tx_text_buf ); \
                char *test_result = NULL;\
                size_t name_size = strlen(tx_tests[i].name);\
     \
                if ( !failed ) {\
                    tx_tests[i].f();\
                    test_result = (char *) TX_GREEN("OK");\
                }\
                else {\
                    printf("Test failed %s\n", tx_error_buffer);\
                    test_result = (char *) TX_RED("FAILED");\
                    success += 1;\
                }\
    \
                printf("%s ", tx_tests[i].name);\
                printf("\033[90m");\
                for (size_t i = 0; i < TX_TEST_OUTPUT_LEN - name_size; i++ ) {\
                    putchar('.');\
                }\
                printf("\033[0m");\
                printf(" %s\n", test_result);\
                fflush(stdout);\
            }\
        }\
\
        return success;\
    }\

