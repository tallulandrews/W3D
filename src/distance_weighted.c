#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "misc_functions.h"

/* 
	Implementation by: Tallulah Andrews
	Date: 28 Sept 2015
	Untested
*/


/* Cannot pass matrices to C from R instead matrix = long array of col1,col2,col3,etc...*/
void distance_wgt (double* y, double* w, int* nrow, int* ncol, double* exponent, double* out) {
	int row1,row2, column;
	for (row1 = 0; row1 < *nrow; row1++) {
		for (row2 = 0; row2 < *nrow; row2++) {
			double dist = 0.0;
			for ( column = 0; column < *ncol; column++) {
				double difference = y[convert_2D_indices_to_1D(row1, column, nrow, ncol)]-y[convert_2D_indices_to_1D(row2, column, nrow, ncol)];
				if (*exponent == 2.0) {
					dist = dist + w[convert_2D_indices_to_1D(row1, column, nrow, ncol)]*difference*difference;
				} else if (*exponent == 1.0) {
					dist = dist + w[convert_2D_indices_to_1D(row1, column, nrow, ncol)]*abs(difference);
				} else {
					dist = dist + w[convert_2D_indices_to_1D(row1, column, nrow, ncol)]*pow(abs(difference), *exponent);
				}
			}

			if (*exponent == 2.0) {
                                dist = sqrt(dist);
                        } else if (*exponent == 1.0) {
                                dist = dist;
                        } else {
                                dist = pow(dist, 1.0/(*exponent));
                        }
			out[convert_2D_indices_to_1D(row1, row2, nrow, nrow)] = dist;
		}
	}
}
