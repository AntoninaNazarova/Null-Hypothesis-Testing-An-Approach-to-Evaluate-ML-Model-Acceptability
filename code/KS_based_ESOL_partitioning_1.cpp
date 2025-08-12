//The following is a code snippet of the Kennard-Stone (KS)-based dataset partitioning functionality

// Constants for dataset sizes
const int total_samples = 1128;
const int evaluation_set_size = 55;
const int training_set_size = total_samples - evaluation_set_size;

// Initialize sample indices
for (int i3 = 0; i3 < total_samples; i3++) {
    i_sm11[i3] = i3; // Use declared int array i_sm11
}

// Step 1: Find the most distant pair of samples
float max = 0.0; // Use float for max as per the provided type
int i1_max = 0, i2_max = 0;

for (int i1 = 0; i1 < total_samples; i1++) {
    int i11 = i_sm11[i1];
    for (int i2 = 0; i2 < total_samples; i2++) {
        int i22 = i_sm11[i2];
        float sum = 0.0; // Use float for sum as per the provided type

        // Calculate the squared distance between two samples
        for (int i = 0; i < 784 / ipr; i++) {
            float diff = s_smile[i11][i] - s_smile[i22][i];
            sum += diff * diff;
        }

        // Update maximum distance and indices
        if (max < sum) {
            max = sum;
            i1_max = i1;
            i2_max = i2;
        }
    }
}

// Initialize the sample order array
for (int i3 = 0; i3 < total_samples; i3++) {
    i_sm1[i3] = -1; // Use int array i_sm1
}

printf("\nMost distant pair: i1 = %d, i2 = %d, max distance = %f, i1_max = %d, i2_max = %d", 
       i1, i2, max, i1_max, i2_max);
getch();

// Step 2: Select samples based on the Kennard-Stone algorithm
i_sm1[0] = i1_max;
i_sm1[1] = i2_max;

for (int i1 = 2; i1 < total_samples; i1++) {
    float sum_max = 0.0; // Use float for sum_max
    int i_max = 0;

    // Iterate through all samples to find the next most distant one
    for (int i2 = 0; i2 < total_samples; i2++) {
        int i22 = i2;

        // Skip already selected samples
        int is_selected = 0; // Use int for logical flag
        for (int i3 = 0; i3 < i1; i3++) {
            if (i22 == i_sm1[i3]) {
                is_selected = 1;
                break;
            }
        }
        if (is_selected) continue;

        // Calculate the minimum distance to the current set of selected samples
        float sum_min = 1e6; // Use float for sum_min
        for (int i3 = 0; i3 < i1; i3++) {
            int i33 = i_sm1[i3];
            float sum1 = 0.0; // Use float for sum1

            for (int i = 0; i < 784 / ipr; i++) {
                float diff = s_smile[i22][i] - s_smile[i33][i];
                sum1 += diff * diff;
            }

            if (sum_min > sum1) {
                sum_min = sum1;
            }
        }

        // Update maximum minimum distance and index
        if (sum_max < sum_min) {
            sum_max = sum_min;
            i_max = i2;
        }
    }

    i_sm1[i1] = i_max;
    printf("\nStep %d: Selected sample = %d, max distance = %f", i1, i_max, sum_max);
}

// Step 3: Identify evaluation set samples
for (int i1 = 0; i1 < total_samples; i1++) {
    i_sm[i1] = 2000; // Use int array i_sm
}

int i2 = 0;
for (int i1 = 0; i1 < total_samples; i1++) {
    int is_training = 0; // Use int for logical flag

    for (int i = 0; i < training_set_size; i++) {
        if (i1 == i_sm1[i]) {
            is_training = 1;
            break;
        }
    }

    if (!is_training) {
        i_sm[i2++] = i1; // Add to evaluation set
    }
}

// Debug: Print evaluation and training sets
for (int i1 = 0; i1 < total_samples; i1++) {
    printf("\nSample %d: i_sm = %d", i1, i_sm[i1]);
}

// Write results to the output file
for (int i1 = 0; i1 < total_samples; i1++) {
    int j = i_sm1[i1]; // Use int for j
    fprintf(stream_10, "\n%d", j);
}

// Close all file streams and finalize the program

