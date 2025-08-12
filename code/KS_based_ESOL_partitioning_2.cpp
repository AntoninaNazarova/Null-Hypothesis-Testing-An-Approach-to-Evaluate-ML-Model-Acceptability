//The following is a code snippet of the Kennard-Stone (KS)-based dataset partitioning functionality

//The program randomly selects 55 samples to form the evaluation set and appends them to the 
//end of the output file (write_random_X.dat). This output file contains a list of 1128 sample numbers 
//arranged in the required order for input into the feedforward neural network. The remaining 
//1128 - 55 samples are reorganized based on the Kennard-Stone hierarchy, which was already established 
//using the KS_based_ESOL_partitioning_1.cpp

// Constants
const int total_samples = 1128;
const int evaluation_set_size = 55;
const int training_set_size = total_samples - evaluation_set_size;

// Array to store sample indices
int i_sm1[total_samples];

// Step 1: Initialize all sample indices to -1
for (int i3 = 0; i3 < total_samples; i3++) {
    i_sm1[i3] = -1;
}

// Step 2: Randomly assign unique indices to `i_sm1`
for (int i = 0; i < total_samples; i++) {
    printf("\n i = %d", i);

random_selection:
    int i10 = rand() % total_samples;

    // Check if the index is already selected
    bool is_duplicate = false;
    for (int i1 = 0; i1 < total_samples; i1++) {
        if (i_sm1[i1] == i10) {
            is_duplicate = true;
            break;
        }
    }

    if (is_duplicate) {
        goto random_selection;
    }

    // Assign unique index to `i_sm1`
    i_sm1[i] = i10;
}

// Step 3: Swap indices for evaluation and training sets
for (int i = 0; i < evaluation_set_size; i++) {
    int i1 = i * evaluation_set_size + i;
    int j1 = i_sm1[i1];
    int j2 = i_sm1[total_samples - evaluation_set_size + i];

    // Perform the swap
    i_sm1[i1] = j2;
    i_sm1[total_samples - evaluation_set_size + i] = j1;
}

// Step 4: Reinitialize the first `training_set_size` indices to -1
for (int i1 = 0; i1 < training_set_size; i1++) {
    i_sm1[i1] = -1;
}

// Step 5: Randomly assign unique indices to the training set
for (int i = 0; i < training_set_size; i++) {
    printf("\n i = %d", i);

random_training_selection:
    int i10 = rand() % total_samples;

    // Check if the index is already in the evaluation set
    bool is_in_evaluation = false;
    for (int i1 = training_set_size; i1 < total_samples; i1++) {
        if (i_sm1[i1] == i10) {
            is_in_evaluation = true;
            break;
        }
    }

    if (is_in_evaluation) {
        goto random_training_selection;
    }

    // Check if the index is already in the training set
    bool is_in_training = false;
    for (int i1 = 0; i1 < training_set_size; i1++) {
        if (i_sm1[i1] == i10) {
            is_in_training = true;
            break;
        }
    }

    if (is_in_training) {
        goto random_training_selection;
    }

    // Assign unique index to the training set
    i_sm1[i] = i10;
}

// Print the final results (for debugging purposes)
for (int i1 = 0; i1 < total_samples; i1++) {
    printf("\ni_sm1[%d] = %d", i1, i_sm1[i1]);
}
