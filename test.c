#include <stdio.h> 

typedef int interval_array_t[10][2];
/* typedef int const cinterval_array_t[10][2]; */
typedef int (*const const_interval_array_t)[2];

void func(const_interval_array_t tar) {
    printf("do something");
}

/*
void func(int (* const tar)[2]) {
    printf("do something");
}
*/

int main(void) {
    interval_array_t tar;
    func(tar);
    return 0;
}
