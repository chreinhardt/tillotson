/*
 * This is a simple program to test the Tillotson EOS library.
 *
 * Author:   Christian Reinhardt
 * Date:     22.02.2021
 * Modified:  
 *
 * Test if the error messages generated in tillErrorString() are correct.
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "tillotson.h"

#define max(A,B) ((A) > (B) ? (A) : (B))
#define min(A,B) ((A) > (B) ? (B) : (A))

#define INDEX(i, j) (((i)*tillMat->nTableV) + (j))

int main(int argc, char **argv) {
    char ErrorString[256]; 
    int iError;

    /*
     * Print error messages for all error codes.
     */	
    iError = TILL_SUCCESS;
    tillErrorString(iError, ErrorString); 
    printf("Error: %s (iError= %i)\n", ErrorString, iError);

    iError = TILL_FAIL;
    tillErrorString(iError, ErrorString);
    printf("Error: %s (iError= %i)\n", ErrorString, iError);

    iError = TILL_LOOKUP_SUCCESS;
    tillErrorString(iError, ErrorString);
    printf("Error: %s (iError= %i)\n", ErrorString, iError);

    iError = TILL_LOOKUP_OUTSIDE_RHOMIN;
    tillErrorString(iError, ErrorString);
    printf("Error: %s (iError= %i)\n", ErrorString, iError);

    iError = TILL_LOOKUP_OUTSIDE_RHOMAX;
    tillErrorString(iError, ErrorString);
    printf("Error: %s (iError= %i)\n", ErrorString, iError);

    iError = TILL_LOOKUP_OUTSIDE_VMIN;
    tillErrorString(iError, ErrorString);
    printf("Error: %s (iError= %i)\n", ErrorString, iError);

    iError = TILL_LOOKUP_OUTSIDE_VMAX;
    tillErrorString(iError, ErrorString);
    printf("Error: %s (iError= %i)\n", ErrorString, iError);

    return 0;
}
