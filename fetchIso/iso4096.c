//  Code modified by Anke Friederici
//	Original Copyright 2011 Johns Hopkins University
//
//  Licensed under the Apache License, Version 2.0 (the "License");
//  you may not use this file except in compliance with the License.
//  You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
//  Unless required by applicable law or agreed to in writing, software
//  distributed under the License is distributed on an "AS IS" BASIS,
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  See the License for the specific language governing permissions and
//  limitations under the License.

#include <stdio.h>
#include <time.h>
#include <float.h>
#include <math.h>
#include "turblib.h"

/*
 * Turbulence Database sample C client code
 */
#define N 10

int compare_threshold(const void *a_void, const void *b_void) {
    ThresholdInfo a = *((ThresholdInfo *)a_void);
    ThresholdInfo b = *((ThresholdInfo *)b_void);

    if (a.z != b.z) return a.z < b.z ? -1 : 1;
    if (a.y != b.y) return a.y < b.y ? -1 : 1;
    if (a.x != b.x) return a.x < b.x ? -1 : 1;
    return 0;
}

int main(int argc, char *argv[]) {

    // Settings //
    const char *filename = "../../Data/Percolation/Iso4096/iso64x64x64.raw";
    int x_size = 64;
    int y_size = 64;
    int z_size = 64;
    // Settings //

    char *authtoken = "please provide your own";
    char *dataset = "isotropic4096";

    char *threshold_field = "vorticity";
    ThresholdInfo *vorticity_array; /* dynamic array for the results of vorticity queries */
    int vorticity_array_size;       /* size of the vorticity array */
    float threshold = 0.0f;

    FILE *out_file = fopen(filename, "wb");  // write only

    // test for files not existing.
    if (out_file == NULL) {
        printf("Error! Could not open file\n");
        exit(-1);
    }

    /* Initialize gSOAP */
    soapinit();

    time_t rawtime;
    struct tm *timeinfo;
    time(&rawtime);
    timeinfo = localtime(&rawtime);
    printf("==> Starting at: %s", asctime(timeinfo));

    /* Enable exit on error.  See README for details. */
    turblibSetExitOnError(1);

    float average_total = 0;
    float average_block = 0;

    int max_elements = 1024 * 512;
    int y_step = max_elements / x_size;
    if (y_step > y_size) y_step = y_size;
    int z_step = y_step < y_size ? 1 : max_elements / (x_size * y_size);
    if (z_step > z_size) z_step = z_size;
    int num_requests = (y_size / y_step) * (z_size / z_step);

    printf("\n= Loading (%d, %d, %d), with requests of size (%d, %d, %d)\n", x_size, y_size, z_size,
           x_size, y_step, z_step);
    printf("= Will take %d requests.\n", num_requests);
    for (int z = 0; z < z_size; z += z_step)
        for (int y = 0; y < y_size; y += y_step) {

            // NOTE: The array storing the results is dynamically allocated inside the getThreshold
            // function, because it's size is not known. It needs to be freed after it has been used
            // to avoid leaking the memory.
            getThreshold(authtoken, dataset, threshold_field, 0, threshold, FD4NoInt, 0, y, z,
                         x_size, y_step, z_step, &vorticity_array, &vorticity_array_size);
            printf("\nGot %d vorticity values\n", vorticity_array_size);

            // Sort.
            qsort(vorticity_array, vorticity_array_size, sizeof(ThresholdInfo), compare_threshold);

            for (int p = 0; p < vorticity_array_size; p++) {
                fwrite(&vorticity_array[p].value, sizeof(float), 1, out_file);
                average_block += vorticity_array[p].value;
            }
            // Average per request.
            average_total += average_block / vorticity_array_size;
            average_block = 0;

            // Free the threshold array after using it.
            free(vorticity_array);

            // Timing.
            time(&rawtime);
            timeinfo = localtime(&rawtime);
            printf("Finished downloading at at: %s", asctime(timeinfo));
        }

    float average = average_total / num_requests;
    printf("\n>> Total average: %.6f\n", average);

    // Final time.
    time(&rawtime);
    timeinfo = localtime(&rawtime);
    printf("==> Finished at: %s\n", asctime(timeinfo));

    /* Free gSOAP resources */
    soapdestroy();
    fclose(out_file);

    // Go through all values again to compute rms.
    FILE *in_file = fopen(filename, "rb");
    float *in_data = (float *)malloc(max_elements * sizeof(float));

    float rms_total = 0;
    float rms_block = 0;

    // Go through all requests again, we know that block size makes sense.
    for (int r = 0; r < num_requests; ++r) {
        fread(in_data, sizeof(float), max_elements, in_file);

        // RMS of all errors in block.
        for (int f = 0; f < max_elements; ++f) {
            float error = in_data[f] - average;
            rms_block += error * error;
        }

        rms_total += rms_block / max_elements;
        rms_block = 0;
    }
    float rms = sqrt(rms_total / num_requests);
    printf("\n>> Total RMS: %.6f\n", rms);
    free(in_data);

    // Final time.
    time(&rawtime);
    timeinfo = localtime(&rawtime);
    printf("==> Found rms at: %s\n", asctime(timeinfo));

    fclose(in_file);
    return 0;
}
