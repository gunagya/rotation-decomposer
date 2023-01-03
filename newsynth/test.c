
#include <stdio.h>

#include "include/gridsynth_ccals.h"
#include <HsFFI.h>

int main(int argc, char** argv)
{
    hs_init(&argc, &argv);

    
    const int bsize = 2048;
    char buffer[bsize];

    printf("Executing gridsynth_angle_to_seq_with_precision_hs\n");
    const int result = gridsynth_angle_to_seq_with_precision_hs(10, "3*pi/16", &buffer, bsize);
    printf("Returned %d \n", result);
    printf("Wrote %s \n", buffer);

    hs_exit();

    return 0;
}