#include <stdio.h>
#include "update_state.h"

/**
* 
* Update the state vector based on the given reaction.
* @param x The state vector. Should be of length 7.
* @param reaction The index of the reaction to be performed.
* The state change matrix is a 15x7 matrix, where each row represents a reaction.
* The state change matrix is defined as follows:
*   0: 1 0 0 0 0 0 0
*   1: -1 0 0 0 0 0 0
*   2: -1 0 1 0 0 0 0
*   3: 0 1 0 0 0 0 0
*   4: 0 -1 0 0 0 0 0
*   5: 0 -1 0 1 0 0 0
*   6: 0 0 -1 0 0 0 0
*   7: 0 0 -1 0 1 0 0
*   8: 0 0 0 -1 0 0 0
*   9: 0 0 0 -1 0 1 0
*   10: 0 0 0 0 -1 0 0
*   11: 0 0 0 0 -1 0 1
*   12: 0 0 0 0 0 -1 0
*   13: 1 0 0 0 0 0 -1
*   14: 0 0 0 0 0 0 -1
*/

void update_state(int *x, int reaction)
{
    switch (reaction)
    {
    case 0:
        x[0] += 1;
        break;
    case 1:
        x[0] -= 1;
        break;
    case 2:
        x[0] -= 1;
        x[2] += 1;
        break;
    case 3:
        x[1] += 1;
        break;
    case 4:
        x[1] -= 1;
        break;
    case 5:
        x[1] -= 1;
        x[3] += 1;
        break;
    case 6:
        x[2] -= 1;
        break;
    case 7:
        x[2] -= 1;
        x[4] += 1;
        break;
    case 8:
        x[3] -= 1;
        break;
    case 9:
        x[3] -= 1;
        x[5] += 1;
        break;
    case 10:
        x[4] -= 1;
        break;
    case 11:
        x[4] -= 1;
        x[6] += 1;
        break;
    case 12:
        x[5] -= 1;
        break;
    case 13:
        x[0] += 1;
        x[6] -= 1;
        break;
    case 14:
        x[6] -= 1;
        break;
    default:
        perror("Invalid reaction index");
        break;
    }
}