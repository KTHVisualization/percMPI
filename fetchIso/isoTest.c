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

int compare_threshold(const void* a_void, const void* b_void) {
    ThresholdInfo a = *((ThresholdInfo*)a_void);
    ThresholdInfo b = *((ThresholdInfo*)b_void);

    if (a.z != b.z) return a.z < b.z ? -1 : 1;
    if (a.y != b.y) return a.y < b.y ? -1 : 1;
    if (a.x != b.x) return a.x < b.x ? -1 : 1;
    return 0;
}

int main() {

    // Settings //
#define x_size 50
#define y_size 50
#define z_size 50
    // Settings //

    const int x_full = 4096;
    const int y_full = 4096;
    const int z_full = 4096;

#define size (x_size * y_size * z_size)
    char* authtoken = "please provide your own";
    char* dataset = "isotropic4096";

    char* threshold_field = "vorticity";
    ThresholdInfo* vorticity_array; /* dynamic array for the results of vorticity queries */
    int vorticity_array_size;       /* size of the vorticity array */
    float threshold = 0.0f;

    // Make sample point list.
    float points[size][3];
    float gradients[size][9];

    for (int z = 0; z < z_size; ++z)
        for (int y = 0; y < y_size; ++y)
            for (int x = 0; x < x_size; ++x) {
                int idx = x + x_size * y + x_size * y_size * z;
                points[idx][0] = M_PI_2 * (float)(x / (x_full - 1));
                points[idx][1] = M_PI_2 * (float)(y / (y_full - 1));
                points[idx][2] = M_PI_2 * (float)(z / (z_full - 1));
                // printf("Index %d\n",idx);
            }

    /* Initialize gSOAP */
    soapinit();
    /* Enable exit on error.  See README for details. */
    turblibSetExitOnError(1);

    // Timing.
    time_t start_time, end_time;
    struct tm* timeinfo;
    time(&start_time);
    timeinfo = localtime(&start_time);
    printf("==> Starting at: %s", asctime(timeinfo));

    // Getting velocity.
    {
        printf("Loading %d³ raw velocity.\t", x_size);

        float* rawdata = (float*)malloc(size * sizeof(float) * 3);
        getRawVelocity(authtoken, dataset, 0.0f, 0, 0, 0, x_size, y_size, z_size, (char*)rawdata);

        // Timing.
        time(&end_time);
        printf("Finished downloading #%d in: %.2fs\n", size, difftime(end_time, start_time));
        free(rawdata);
        time(&start_time);
    }

    {
        printf("Loading %d³ block.\t\t", x_size);

        // NOTE: The array storing the results is dynamically allocated inside the getThreshold
        // function, because it's size is not known. It needs to be freed after it has been used to
        // avoid leaking the memory.
        getThreshold(authtoken, dataset, threshold_field, 0, threshold, FD4NoInt, 0, 0, 0, x_size,
                     y_size, z_size, &vorticity_array, &vorticity_array_size);

        // Sort.
        qsort(vorticity_array, size, sizeof(ThresholdInfo), compare_threshold);

        // Timing.
        time(&end_time);
        printf("Finished downloading #%d in: %.2fs\n", vorticity_array_size,
               difftime(end_time, start_time));
        time(&start_time);
    }

    {
        printf("Loading %d³ gradient block.\t", x_size);

        getVelocityGradient(authtoken, dataset, 0, FD4NoInt, NoTInt, size, points, gradients);

        // Timing.
        time(&end_time);
        printf("Finished downloading #%d in: %.2fs\n", size, difftime(end_time, start_time));
        time(&start_time);
    }

    {
        printf("Loading %d³ block in 4 slices.\t", x_size);

        for (int i = 0; i < 4; ++i) {
            // NOTE: The array storing the results is dynamically allocated inside the getThreshold
            // function, because it's size is not known. It needs to be freed after it has been used
            // to avoid leaking the memory.
            getThreshold(authtoken, dataset, threshold_field, 0, threshold, FD4NoInt, 0, 0,
                         z_size / 4 * i, x_size, y_size, z_size / 4, &vorticity_array,
                         &vorticity_array_size);

            // Sort.
            qsort(vorticity_array, vorticity_array_size, sizeof(ThresholdInfo), compare_threshold);
            free(vorticity_array);
        }
        // Timing.
        time(&end_time);
        printf("Finished downloading #%d in: %.2fs\n", vorticity_array_size * 4,
               difftime(end_time, start_time));
        time(&start_time);
    }

    {
        printf("Loading %d³ elements in x.\t", x_size);

        // NOTE: The array storing the results is dynamically allocated inside the getThreshold
        // function, because it's size is not known. It needs to be freed after it has been used to
        // avoid leaking the memory.
        getThreshold(authtoken, dataset, threshold_field, 0, 0, FD4NoInt, 0, 0, 0, x_full,
                     size / x_full, 1, &vorticity_array, &vorticity_array_size);

        // Sort.
        qsort(vorticity_array, vorticity_array_size, sizeof(ThresholdInfo), compare_threshold);
        free(vorticity_array);

        // Timing.
        time(&end_time);
        printf("Finished downloading #%d in: %.2fs\n", vorticity_array_size,
               difftime(end_time, start_time));
        time(&start_time);
    }
    {
        printf("Loading %d³ elements in z.\t", x_size);

        // NOTE: The array storing the results is dynamically allocated inside the getThreshold
        // function, because it's size is not known. It needs to be freed after it has been used to
        // avoid leaking the memory.
        getThreshold(authtoken, dataset, threshold_field, 0, -1.0f, FD4NoInt, 0, 0, 0, 1,
                     size / z_full, z_full, &vorticity_array, &vorticity_array_size);

        // Sort.
        qsort(vorticity_array, vorticity_array_size, sizeof(ThresholdInfo), compare_threshold);
        free(vorticity_array);

        // Timing.
        time(&end_time);
        printf("Finished downloading #%d in: %.2fs\n", vorticity_array_size,
               difftime(end_time, start_time));
        time(&start_time);
    }

    /* Free gSOAP resources */
    soapdestroy();
    return 0;
}
