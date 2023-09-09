// tx_tests.hpp

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <exception>
#include <assert.h>


class tx_test_exception : public std::exception {
public:
    char *message;
    char *location;

    tx_test_exception( void ) : message(NULL), location(NULL) {}
    tx_test_exception( char * msg ) : message(msg), location(NULL) {}
    tx_test_exception( char * msg, char *loc ) : message(msg), location(loc) {}

    char *what ( void ) {
        return message;
    }
};


void tx_assert_msg(
    int cond, char *msg, const char *cond_str, const char *func, const char *file, const char *line
) {
    if ( !cond ) {
        char *buff = (char *) malloc(1024);
        snprintf(buff, 1024, "(%s), function %s, file %s, line %s", cond_str, func, file, line);
        throw tx_test_exception( msg, buff );
    }
}


#define TX_AS_STRING2(x) #x
#define TX_AS_STRING(x) TX_AS_STRING2(x)
#define TX_LINE TX_AS_STRING(__LINE__)

#define TX_ASSERT_MSG(cond, msg) tx_assert_msg(cond, msg, #cond, __func__, __FILE__, TX_LINE)
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
                size_t name_size = strlen(tx_tests[i].name);\
                char *test_result = NULL;\
\
                try {\
                    tx_tests[i].f();\
                    test_result = (char *) TX_GREEN("OK");\
                }\
                catch (tx_test_exception &e) {\
                    test_result = (char *) TX_RED("FAILED");\
                    printf("Test failed %s\n", e.location);\
                    success = 1;\
                }\
                catch (...) {\
                    test_result = (char *) TX_RED("FAILED");\
                    printf("Test function could not be successfully executed!\n");\
                    success = 1;\
                }\
\
                printf("%s ", tx_tests[i].name);\
                printf("\033[90m");\
                for (size_t i = 0; i < TX_TEST_OUTPUT_LEN - name_size; i++ ) {\
                    putchar('.');\
                }\
                printf("\033[0m");\
                printf(" %s\n", test_result);\
            }\
        }\
\
        return success;\
    }\
