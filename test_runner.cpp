// test_runner.cpp

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <exception>
#include <assert.h>


class test_exception : public std::exception {
public:
    char *message;
    char *location;

    test_exception( void ) : message(NULL), location(NULL) {}
    test_exception( char * msg ) : message(msg), location(NULL) {}
    test_exception( char * msg, char *loc ) : message(msg), location(loc) {}

    char *what ( void ) {
        return message;
    }
};


void assert_msg(
    int cond, char *msg, const char *cond_str, const char *func, const char *file, const char *line
) {
    if ( !cond ) {
        char *buff = (char *) malloc(1024);
        snprintf(buff, 1024, "(%s), function %s, file %s, line %s", cond_str, func, file, line);
        throw test_exception( msg, buff );
    }
}


#define AS_STRING2(x) #x
#define AS_STRING(x) AS_STRING2(x)
#define LINE AS_STRING(__LINE__)

#define ASSERT_MSG(cond, msg) assert_msg(cond, msg, #cond, __func__, __FILE__, LINE)
#define ASSERT(cond) ASSERT_MSG(cond, NULL)

#define TEST_FUN(name) void (*name) (void)


struct test {
    TEST_FUN(f);
    const char *name;
};


#define N_TESTS 1000
#define ADD_TEST(func) (test) { .f = func, .name = #func }

// #include "tests.cpp"

void test_1( void )
{
    assert(1);
    ASSERT(1 == 0);
}


test tests[N_TESTS] = {
    ADD_TEST(test_1),
};


#define TEST_OUTPUT_LEN 70
#define GREEN(text) "\033[32m" text "\033[0m"
#define RED(text)   "\033[31m" text "\033[0m"


int main( int argn, const char *argv[] )
{
    int success = 0;

    for ( int i=0; i<N_TESTS; i++ ) {
        if ( tests[i].f ) {
            size_t name_size = strlen(tests[i].name);
            char *test_result = NULL;

            try {
                tests[i].f();
                test_result = (char *) GREEN("OK");
            }
            catch (test_exception &e) {
                test_result = (char *) RED("FAILED");
                printf("Test failed %s\n", e.location);
                success = 1;
            }
            catch (...) {
                test_result = (char *) RED("FAILED");
                printf("Test function could not be executed!\n");
                success = 1;
            }

            printf("%s ", tests[i].name);
            printf("\033[90m");
            for (size_t i = 0; i < TEST_OUTPUT_LEN - name_size; i++ ) {
                putchar('.');
            }
            printf("\033[0m");
            printf(" %s\n", test_result);
        }
    }

    return 0;
}
