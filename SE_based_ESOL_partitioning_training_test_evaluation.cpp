//The following is a code snippet of the sphere exclusion (SE)-based dataset partitioning functionality

// Initialization of range variables
double rmin = 10000.0; // Minimum distance (initialized to a large value)
double rmax = 0.0;     // Maximum distance (initialized to 0)

// Calculate minimum and maximum distances between all pairs of elements
for (int i1 = 0; i1 < 1128 - 1; i1++) {
    for (int i2 = i1 + 1; i2 < 1128; i2++) {
        double r = 0.0;
        for (int i = 0; i < 784 / ipr; i++) {
            double diff = s_smile[i1][i] - s_smile[i2][i];
            r += diff * diff;
        }
        if (r < rmin) rmin = r; // Update minimum distance
        if (r > rmax) rmax = r; // Update maximum distance
    }
}

// Set threshold radius for clustering
double rfor = rmin + (rmax - rmin) / 160000000.0;
rfor = 7.0; // Adjusted based on a specific use case

// ******************* CLUSTERING INITIALIZATION *******************
int i20 = 1128; // Total number of elements
int clast_l[1128], clast_t[1128], clast_cc[1128];
for (int i1 = 0; i1 < i20; i1++) {
    clast_l[i1] = -1;
    clast_t[i1] = -1;
    clast_cc[i1] = -1;
}

int i20_l = 0;      // Number of centers for spheres
int i20_t = 0;      // Test set counter
int i20_c = 0;      // Cluster counter
int i20_test = 0;   // Total test set size

// ******************* RANDOM SELECTION OF SPHERE CENTER *******************
clstt3:
int j;
bool is_duplicate;
do {
    j = rand() % i20; // Select a random index
    is_duplicate = false;

    // Check if the selected index is already in any cluster or test set
    for (int i3 = 0; i3 < i20_c; i3++) {
        if (j == clast_cc[i3]) {
            is_duplicate = true;
            break;
        }
    }
    for (int i3 = 1128 - kp; i3 < 1128; i3++) {
        if (j == i_sm1[i3]) {
            is_duplicate = true;
            break;
        }
    }
} while (is_duplicate);

// Assign selected index as a sphere center
clast_l[i20_l++] = j; // Add to sphere centers
clast_cc[i20_c++] = j; // Add to all clusters

// Initialize cluster center with selected index
for (int i = 0; i < l1; i++) {
    f_clast[i] = s_smile[j][i];
}

// ******************* CLUSTERING LOGIC *******************
clstt2:
int i_count = 1; // Count of elements in the current cluster
iss[0] = j;      // Initialize the first cluster element

for (int i1 = 0; i1 < i20; i1++) {
    bool is_in_cluster = false;

    // Check if current index is already in any cluster
    for (int i3 = 0; i3 < i20_c; i3++) {
        if (i1 == clast_cc[i3]) {
            is_in_cluster = true;
            break;
        }
    }
    for (int i3 = 1128 - kp; i3 < 1128; i3++) {
        if (i1 == i_sm1[i3]) {
            is_in_cluster = true;
            break;
        }
    }
    if (is_in_cluster) continue;

    // Calculate the distance to the cluster center
    double r = 0.0;
    for (int i = 0; i < l1; i++) {
        double diff = s_smile[i1][i] - f_clast[i];
        r += diff * diff;
    }

    // Add element to cluster if within the radius
    if (r < rfor) {
        clast_cc[i20_c++] = i1;
        iss[i_count++] = i1;
    }
}

// Add cluster elements to the test set
int split_count = static_cast<int>(0.963 * (i_count - 1));
if ((0.963 * (i_count - 1) - split_count) > 0.5) split_count++;

is_test[i20_test++] = iss[0]; // Add the cluster center to the test set
for (int i = 1; i < split_count; i++) {
    is_test[i20_test++] = iss[i];
}

printf("\n____ i20_test = %d", i20_test);

// ******************* LEARNING SET GENERATION *******************
int is_learning[1128];
int i20_learning = 0;
for (int i1 = 0; i1 < 1128; i1++) {
    bool is_in_test_or_cluster = false;

    // Check if the index is part of the test set or cluster
    for (int i3 = 0; i3 < i20_test; i3++) {
        if (i1 == is_test[i3]) {
            is_in_test_or_cluster = true;
            break;
        }
    }
    for (int i3 = 1128 - kp; i3 < 1128; i3++) {
        if (i1 == i_sm1[i3]) {
            is_in_test_or_cluster = true;
            break;
        }
    }
    if (!is_in_test_or_cluster) {
        is_learning[i20_learning++] = i1; // Add to learning set
    }
}

// Merge test and learning sets into `i_sm1`, because the first 1128-55 will be training, 
// followed by 55 for testing, and the final 55 for evaluation
for (int i = 0; i < i20_test; i++) {
    i_sm1[i] = is_test[i];
}
for (int i = 0; i < i20_learning; i++) {
    i_sm1[i20_test + i] = is_learning[i];
}

// ******************* DUPLICATE VERIFICATION *******************
for (int i1 = 0; i1 < 1128 - 1; i1++) {
    for (int i2 = i1 + 1; i2 < 1128; i2++) {
        if (i_sm1[i1] == i_sm1[i2]) {
            printf("\nDuplicate found.");
            exit(1);
        }
    }
}

// Print results for debugging
for (int i1 = 0; i1 < 1128; i1++) {
    printf("\ni1 = %d, i_sm1 = %d", i1, i_sm1[i1]);
}