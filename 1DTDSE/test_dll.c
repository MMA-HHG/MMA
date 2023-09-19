#include <stdlib.h>
#include <stdio.h>

// Compile using ```cc -fPIC -shared -o test_dll.so test_dll.c```

// Define custom structure for testing
typedef struct nested {
    int test;
} nested;

typedef struct str {
    nested nested;
    double * arr;
    double * arr2;
    int arr_size;
    int number;
    double * float_ptr;
} str;


// Declare function working with the custom structure
int test_struct(str *s) {
    s->arr2 = calloc(s->arr_size, sizeof(double));
    for (int i = 0; i < s->arr_size; i++) {
        printf("arr[%d] = %f \n", i, s->arr[i]);
        s->arr2[i] = s->arr[i];
    }
    s->arr2[0] = 0;

    printf("Number = %d \n", s->number);

    printf("Float = %f \n", *(s->float_ptr));

    printf("Int = %d\n", s->nested.test);

    

    return s->number;
}