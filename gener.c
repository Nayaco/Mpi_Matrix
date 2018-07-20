#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#define ran(x) (rand()%x)
int n,m;
int main(int argc, char const *argv[])
{                  
    srand(time(NULL));
    scanf("%d %d", &n, &m);                                  
    FILE *file = fopen(argv[1], "w");
    fprintf(file, "%d %d\n", n, m);
    for(int i = 0; i < n; i++){
        for(int j = 0; j < m; j++){
            fprintf(file, "%d ", ran(100));
        }
        fprintf(file,"\n");
    }
    return 0;
}
