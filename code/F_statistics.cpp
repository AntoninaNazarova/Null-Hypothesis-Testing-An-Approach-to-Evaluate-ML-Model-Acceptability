#include <iostream>
#include <vector>
#include <array>
#include <cmath>
#include <string>
#include <random>
#include <fstream>
#include <memory>
#include <algorithm>
#include <unordered_set>

template <typename T>
constexpr const T& clamp(const T& v, const T& lo, const T& hi) {
    return (v < lo) ? lo : (hi < v) ? hi : v;
}

int main() {
	std::cout << "Program started.\n";
    constexpr int nb = 12000000;
    constexpr int num_samples = 383;
    constexpr int feature_length = 2160;
    constexpr int kp = 20;
    constexpr int ipr = 4;
    constexpr int max_rnd = 10000;
    constexpr int max_size = 1000;
    constexpr int l1 = feature_length / ipr;
    constexpr int l2 = l1;
    constexpr int l3 = 1;
	constexpr int i200 = 210;
	
	float dd = 0.0f;
	int i1 = 0, i2 = 0;
    int ip = 0;
	int imin = 0;
	int j2 = 0;
	float fisher1 = 1000.0f;
    int p = 0, ipp = 0, jsum = 0;
    float delta0 = 0.01f, delta_max = 50.0f, delta_min = 1e-6f;
    float eta_plus = 1.2f, eta_minus = 0.5f, ap = 0.95f;
    float w_fisher = 1000.0f, w_Er2 = 0.0f, Er = 0.0f;
    float pi = 3.14159264f;
    float ap_min = 1000.0f, aak1 = 100.0f, max_eval = 10000.0f;
    float result_fisher = 0.0f, result_fisher_test = 0.0f;
    float sr = 0.0f, sum = 0.0f, amin = 10000.f, amin1 = 10000.f;
    float min = 1000.f;
    float metric_test = 0.0f, R2_min = 0.0f, R2_max = 0.0f, ER2 = 0.0f;
    float rm1 = 0.0f, b_t = 0.0f, a_t = 0.0f;
    float gg1 = 0.0f, s_a_av = 0.0f, s_eps_av = 0.0f;
	double rm = 0.0;

    int ipp_count = 0, i_counter_eval = 0;
    int count_fish = 0, count1_fish = 0, count_R2 = 0;
    int iiii = 0, i_fish1 = 0;

    std::array<char, nb> buf{};
    std::array<char, 20> c = {'0', '1', '2', '3', '4', '5', '6', '7', '8', '9'};
    std::array<float, 20> tt{}, ttt{};
    std::array<int, 5000> i_sm1{}, i_sm2{}, i_sm1B{}, i_sm2B{};
    std::array<float, 5000> s_eps{}, s_eps1{};
    std::array<float, 8000> s_a{}, s_a1{};
    std::array<float, 10000> regr{};
	std::array<float, 5000> eps_train{}, eps_evaluation{}, eps_test{};

	std::vector<std::vector<float>>
		s_smile(5000, std::vector<float>(5000)),
		s_smile1(5000, std::vector<float>(5000));
        
	std::cout << "Trying to open input/output files...\n";
    std::ifstream stream_A("s_smiles.dat", std::ios::binary);
    std::ifstream stream_Z("s_eps.dat", std::ios::binary);
    std::ifstream stream_10("write_4_random_1.dat", std::ios::binary);
    std::ifstream stream_B10("write_4.dat", std::ios::binary);

    std::ofstream stream_1("1train_4_RELU_RN_KS5.dat");
    std::ofstream stream_2("2test_4_RELU_RN_KS5.dat");
    std::ofstream stream_3("3dd_4_RELU_RN_KS5.dat");
    std::ofstream stream_4("4b_t_4_RELU_RN_KS5.dat");
    std::ofstream stream_5("5sigma_4_RELU_RN_KS5.dat");
    std::ofstream stream_6("6R2_4_RELU_RN_KS5.dat");
    std::ofstream stream_7("7a_t_4_RELU_RN_KS5.dat");
    std::ofstream stream_9("9fisher_4_RELU_RN_KS5.dat");
    std::ofstream stream_11("10Es_4_RELU_RN_KS5.dat");
    std::ofstream stream_12("11ER2_4_RELU_RN_KS5.dat");
    std::ofstream stream_13("12RMSE_4_RELU_RN_KS5.dat");
    std::ofstream stream_14("fa_st2.dat");
    std::ofstream stream_15("r2_st2.dat");

    if (!stream_A || !stream_Z || !stream_10 || !stream_B10 ||
        !stream_1 || !stream_2 || !stream_3 || !stream_4 || !stream_5 ||
        !stream_6 || !stream_7 || !stream_9 || !stream_11 || !stream_12 ||
        !stream_13 || !stream_14 || !stream_15) {
        std::cerr << "Error opening one or more files.\n";
        return 1;
    }

    std::cout << "All files opened successfully.\n";
	
    for (int i = 0; i < 20; ++i) {
        tt[i] = std::pow(10.0f, i);
        ttt[i] = std::pow(10.0f, -i - 1);
    }

        // Fill power tables
    for (int i = 0; i < 20; ++i) {
        tt[i] = std::pow(10.0f, i);
        ttt[i] = std::pow(10.0f, -i - 1);
    }

    // Read SMILES binary
    stream_A.read(buf.data(), nb);
    int ii = 0;
    for (int i = 0; i < num_samples; ++i) {
        int ii1 = 0;
        while (ii < nb) {
            char ch = buf[ii];
            if (ch == '1') s_smile[i][ii1++] = 1.0f, ++ii;
            else if (ch == '-') s_smile[i][ii1++] = -1.0f, ii += 2;
            else if (ch == '0') s_smile[i][ii1++] = 0.0f, ++ii;
            else if (ch == '\r' && buf[ii+1] == '\n') { ii += 2; break; }
        }
    }

    for (int i1 = 0; i1 < num_samples; ++i1)
        for (int i = 0; i < feature_length / ipr; ++i) {
            float r3 = s_smile[i1][i * 4 + 0];
            float r2 = s_smile[i1][i * 4 + 1];
            float r1 = s_smile[i1][i * 4 + 2];
            float r  = s_smile[i1][i * 4 + 3];

            if (r3 < 0) r3 = 2.0f;
            if (r2 < 0) r2 = 2.0f;
            if (r1 < 0) r1 = 2.0f;
            if (r  < 0) r  = 2.0f;

            float r20 = r3 * 8 + r2 * 4 + r1 * 2 + r;
            if (r20 == 0.0f) { s_smile[i1][i] = 0.0f; continue; }

            if (r3 == 2.0f) r3 = 0.0f;
            if (r2 == 2.0f) r2 = 0.0f;
            if (r1 == 2.0f) r1 = 0.0f;
            if (r  == 2.0f) r  = 0.0f;

            r20 = r3 * 8 + r2 * 4 + r1 * 2 + r - 7.5f;
            s_smile[i1][i] = r20;
        }

    for (int i = 0; i < num_samples; ++i)
        std::copy(s_smile[i].begin(), s_smile[i].end(), s_smile1[i].begin());

    // Decode EPS
    stream_Z.read(buf.data(), nb);
    i1 = 0, i2 = 0;
    while (i2 < num_samples && i1 < nb) {
        int start = i1;
        while (i1 < nb && buf[i1] != '\n' && buf[i1] != '\r') {
            ++i1;
        }

        try {
            std::string token(buf.data() + start, i1 - start);
            s_eps[i2++] = std::stof(token);
        } catch (...) {
            std::cerr << "Failed to parse EPS value at i2 = " << i2 << ": "
                      << std::string(buf.data() + start, i1 - start) << "\n";
            return 1;
        }

        // Skip newline characters
        while (i1 < nb && (buf[i1] == '\n' || buf[i1] == '\r')) ++i1;
    }

    for (int i = 0; i < num_samples; ++i) sr += s_eps[i];
    sr /= static_cast<float>(num_samples);
    for (int i = 0; i < num_samples; ++i) s_eps[i] -= sr;

    std::cout << "i2 = " << i2 << " | mean = " << sr << "\n";


    // Power table for decoding

    // Read KS input buffer
    
    if (!stream_B10) {
        std::cerr << "Error opening write_4.dat\n";
        return 1;
    }

    
    stream_B10.read(buf.data(), nb);
    if (!stream_B10) {
        std::cerr << "Error reading buffer.\n";
        return 1;
    }

    // Decode integer values into i_sm2B[]
    
    i1 = 0, i2 = 0;

    while (i2 < num_samples && i1 < nb) {
        int i_count1 = 0, j1 = 0;

        while (true) {
            if (buf[i1] == ',') j1 = i_count1;
            if (buf[i1] == '\r' && buf[i1 + 1] == '\n') {
                i1 += 2;
                break;
            }
            ++i1;
            ++i_count1;
        }

        float r = 0.0f, r2 = 1.0f;
        for (int i = 0; i < j1; ++i) {
            char ch = buf[i1 - i_count1 + j1 - i - 1];
            if (ch == '-') r2 = -1.0f;
            else if (std::isdigit(ch)) r += (ch - '0') * tt[i];
        }

        i_sm2B[i2++] = static_cast<int>(r * r2);
    }


    // Random number setup
    std::mt19937 rng(std::random_device{}());
    std::uniform_int_distribution<int> uni(0, num_samples - 1);
    std::array<int, max_rnd> rnd{};
    rnd.fill(-1);
    iiii = 0;

    // ====== BEGINNING of randomized KS sampling ======
    p = 0;
    R2_max = 0.0f;
    w_fisher = 1000.0f, w_Er2 = 0.0f;

     // <<<<<<<<<<<<<<<<<< USER-DEFINED START SEED
    rnd[iiii++] = i200;

    std::cout << "Starting KS with i200 = " << i200 << "\n";


    std::copy(i_sm2B.begin(), i_sm2B.end(), i_sm1B.begin());


    for (int i = num_samples; i < num_samples + kp; ++i) {
        i_sm1B[i] = -1;
    }

    // Fill kp random samples without duplication
    for (int i = num_samples; i < num_samples + kp; ++i) {
        int i10;
        while (true) {
            i10 = uni(rng);
            bool duplicate = false;
            for (int j = num_samples; j < i; ++j) {
                if (i_sm1B[j] == i10) {
                    duplicate = true;
                    break;
                }
            }
            if (!duplicate) break;
        }
        i_sm1B[i] = i10;
    }

    // Invalidate original training points that appear in extra randoms
    for (int i1 = 0; i1 < num_samples; ++i1) {
        for (int i = num_samples; i < num_samples + kp; ++i) {
            if (i_sm1B[i1] == i_sm1B[i]) {
                i_sm1B[i1] = -1;
            }
        }
    }

    // Build final compacted i_sm1 array
    int idx = 0;
    for (int i1 = 0; i1 < num_samples + kp; ++i1) {
        if (i_sm1B[i1] != -1) {
            i_sm1[idx++] = i_sm1B[i1];
        }
    }

    std::cout << "Final KS + randomized size = " << idx << " (should be 383)\n";
    std::cout << "KS complete.\n";


    // Copy s_smile1 → s_smile
    for (int i = 0; i < num_samples; ++i)
        std::copy(s_smile1[i].begin(), s_smile1[i].end(), s_smile[i].begin());

    // Scalars
    ip = 0, ipp = 0, jsum = 0;
    
    amin = 10000.f, amin1 = 10000.f;
    imin = 0;
    min = 1000.f;
    

    // Error and mean
    int count_sm1 = num_samples - kp - kp;
    int count_sm = kp;

    
    for (int i1 = 383 - kp - kp; i1 < 383 - kp; ++i1) {
        int j = i_sm1[i1];
        float r1 = s_eps[j];
        sum += r1;
    }
    sum /= static_cast<float>(kp);

    Er = 0.0f;
    for (int i1 = 383 - kp - kp; i1 < 383 - kp; ++i1) {
        int j = i_sm1[i1];
        float r1 = s_eps[j];
        Er += (r1 - sum) * (r1 - sum);
    }

    // Initialize deltas and derivatives
    auto make_matrix = [](int rows, int cols, float fill) {
        return std::vector<std::vector<float>>(rows, std::vector<float>(cols, fill));
    };

    auto del1x     = make_matrix(max_size, max_size, 0.0f);
    auto delta1x   = make_matrix(max_size, max_size, delta0);
    auto deltaw1x  = make_matrix(max_size, max_size, delta0);
    auto del2x     = make_matrix(max_size, max_size, 0.0f);
    auto delta2x   = make_matrix(max_size, max_size, delta0);
    auto deltaw2x  = make_matrix(max_size, max_size, delta0);
    auto del3x     = make_matrix(max_size, max_size, 0.0f);
    auto delta3x   = make_matrix(max_size, max_size, delta0);
    auto deltaw3x  = make_matrix(max_size, max_size, delta0);
    auto del1p     = make_matrix(max_size, max_size, 0.0f);
    auto delta1p   = make_matrix(max_size, max_size, delta0);
    auto deltaw1p  = make_matrix(max_size, max_size, delta0);
	auto del11 = make_matrix(max_size, max_size, 0.0f);
	auto del22 = make_matrix(max_size, max_size, 0.0f);
	auto del33 = make_matrix(max_size, max_size, 0.0f);

    // Random weight initialization in [-0.045, 0.045]
    
    std::uniform_real_distribution<float> dist(-0.045f, 0.045f);

    auto w1 = make_matrix(l1, l1, 0.0f);
    auto w2 = make_matrix(l1, l1, 0.0f);
    auto w3 = make_matrix(l3, l1, 0.0f);

    for (int j = 0; j < l1; ++j)
        for (int i = 0; i < l1; ++i)
            w1[j][i] = 3.f * 0.03f * (dist(rng));

    for (int j = 0; j < l1; ++j)
        for (int i = 0; i < l1; ++i)
            w2[j][i] = 3.f * 0.03f * (dist(rng));

    for (int j = 0; j < l3; ++j)
        for (int i = 0; i < l1; ++i)
            w3[j][i] = 3.f * 0.03f * (dist(rng));

    std::cout << "Network initialized: layers l1=" << l1 << ", l2=" << l2 << ", l3=" << l3 << "\n";
    std::cout << "Error mean = " << sum << ", Error variance = " << Er << "\n";
	
	j2 = 0;
	fisher1 = 1000.0f;
	R2_max = 0.0f;
	
	
	ip = 0;
	std::cout << "Beginning training loop...\n";
	do {
		std::cout << "Epoch " << ip << "\n";
		ipp++;

		// Zero initialize deltas
		for (int j = 0; j < 1000; ++j)
			for (int i = 0; i < 1000; ++i)
				del11[j][i] = del22[j][i] = del33[j][i] = 0.0;

		// Training loop
		for (int i111 = 0; i111 < 383 - 2 * kp; ++i111) {
			int i11 = i_sm1[i111];
			double d0 = s_eps[i11];

			std::array<double, 1000> x{}, x1{}, x11{}, x3{}, s1{}, s3{};
			for (int i = 0; i < l1; ++i)
				x[i] = s_smile[i11][i] * ap;

			// Hidden layer
			for (int j = 0; j < l1; ++j) {
				s1[j] = 0.0;
				for (int i = 0; i < l1; ++i)
					s1[j] += w1[j][i] * x[i];
				x11[j] = 0.5 * (1.0 + std::tanh(s1[j] / 2.0));
			}
			std::copy(x11.begin(), x11.end(), x1.begin());

			// Output layer
			for (int j = 0; j < l3; ++j) {
				s3[j] = 0.0;
				for (int i = 0; i < l1; ++i)
					s3[j] += w3[j][i] * x1[i];
				x3[j] = s3[j];
			}

			// del3 and gradient computation
			std::array<double, 1000> del3{};
			for (int i = 0; i < l3; ++i)
				del3[i] = x3[i] - d0;

			for (int j = 0; j < l1; ++j)
				for (int i = 0; i < l1; ++i)
					del11[j][i] += del3[0] * w3[0][j] * x1[j] * (1.0 - x1[j]) * x[i];

			for (int i = 0; i < l2; ++i)
				del33[0][i] += del3[0] * x1[i];
		}

		// RPROP weight update: W1
		for (int j = 0; j < l1; ++j) {
			for (int i = 0; i < l1; ++i) {
				double rr1 = del1x[j][i];
				double rr2 = del11[j][i];

				if (rr1 * rr2 > 0.0) {
					delta1x[j][i] = clamp(delta1x[j][i] * eta_plus, 0.0f, delta_max);
					double sign = (rr2 < 0.0) ? -1.0 : 1.0;
					deltaw1x[j][i] = -sign * delta1x[j][i];
					w1[j][i] += deltaw1x[j][i];
				} else if (rr1 * rr2 < 0.0) {
					delta1x[j][i] = clamp(delta1x[j][i] * eta_minus, delta_min, delta_max);
					w1[j][i] -= deltaw1x[j][i];
					del11[j][i] = 0.0;
				} else {
					double sign = (rr2 == 0.0) ? 0.0 : (rr2 < 0.0 ? -1.0 : 1.0);
					deltaw1x[j][i] = -sign * delta1x[j][i];
					w1[j][i] += deltaw1x[j][i];
				}

				del1x[j][i] = del11[j][i];
			}
		}

		// RPROP weight update: W3
		for (int j = 0; j < l3; ++j) {
			for (int i = 0; i < l2; ++i) {
				double rr1 = del3x[j][i];
				double rr2 = del33[j][i];

				if (rr1 * rr2 > 0.0) {
					delta3x[j][i] = clamp(delta3x[j][i] * eta_plus, 0.0f, delta_max);
					double sign = (rr2 < 0.0) ? -1.0 : 1.0;
					deltaw3x[j][i] = -sign * delta3x[j][i];
					w3[j][i] += deltaw3x[j][i];
				} else if (rr1 * rr2 < 0.0) {
					delta3x[j][i] = clamp(delta3x[j][i] * eta_minus, delta_min, delta_max);
					w3[j][i] -= deltaw3x[j][i];
					del33[j][i] = 0.0;
				} else {
					double sign = (rr2 == 0.0) ? 0.0 : (rr2 < 0.0 ? -1.0 : 1.0);
					deltaw3x[j][i] = -sign * delta3x[j][i];
					w3[j][i] += deltaw3x[j][i];
				}

				del3x[j][i] = del33[j][i];
			}
		}

		// ---------------------------- TRAINING PHASE --------------------------------
		{
			double sum = 0.0;
			for (int i1 = 0; i1 < 383 - 2 * kp; ++i1) {
				int j = i_sm1[i1];
				double d0 = s_eps[j];

				std::array<double, 1000> x{}, x11{}, x3{}, s1{}, s3{};

				for (int i = 0; i < l1; ++i)
					x[i] = s_smile[j][i] * ap;

				// HIDDEN LAYER
				for (int j_ = 0; j_ < l1; ++j_) {
					s1[j_] = 0.0;
					for (int i = 0; i < l1; ++i)
						s1[j_] += w1[j_][i] * x[i];
					x11[j_] = 0.5 * (1.0 + std::tanh(s1[j_] / 2.0));
				}

				// OUTPUT LAYER
				for (int j_ = 0; j_ < l3; ++j_) {
					s3[j_] = 0.0;
					for (int i = 0; i < l1; ++i)
						s3[j_] += w3[j_][i] * x11[i];
					x3[j_] = s3[j_];
				}

				for (int i = 0; i < l3; ++i)
					sum += std::pow(d0 - x3[0], 2.0);

				eps_train[i1] = d0;
				s_a[i1] = x3[0];
			}

			double rmse = std::sqrt(sum / (383.0 - 2 * kp));
			if (amin > rmse)
				amin = rmse;
		}

		// ---------------------------- EVALUATION PHASE -------------------------------
		{
			double sum = 0.0;
			for (int i1 = 383 - kp; i1 < 383; ++i1) {
				int j = i_sm1[i1];
				double d0 = s_eps[j];

				std::array<double, 1000> x{}, x11{}, x3{}, s1{}, s3{};

				for (int i = 0; i < l1; ++i)
					x[i] = s_smile[j][i] * ap;

				for (int j_ = 0; j_ < l1; ++j_) {
					s1[j_] = 0.0;
					for (int i = 0; i < l1; ++i)
						s1[j_] += w1[j_][i] * x[i];
					x11[j_] = 0.5 * (1.0 + std::tanh(s1[j_] / 2.0));
				}

				for (int j_ = 0; j_ < l3; ++j_) {
					s3[j_] = 0.0;
					for (int i = 0; i < l1; ++i)
						s3[j_] += w3[j_][i] * x11[i];
					x3[j_] = s3[j_];
				}

				for (int i = 0; i < l3; ++i)
					sum += std::pow(d0 - x3[0], 2.0);

				eps_evaluation[i1] = d0;
				s_a[i1] = x3[0];
			}
		}

		// ---------------------------- TEST PHASE -------------------------------------
		{
			double sum = 0.0;
			for (int i1 = 383 - 2 * kp; i1 < 383 - kp; ++i1) {
				int j = i_sm1[i1];
				double d0 = s_eps[j];

				std::array<double, 1000> x{}, x11{}, x3{}, s1{}, s3{};

				for (int i = 0; i < l1; ++i)
					x[i] = s_smile[j][i] * ap;

				for (int j_ = 0; j_ < l1; ++j_) {
					s1[j_] = 0.0;
					for (int i = 0; i < l1; ++i)
						s1[j_] += w1[j_][i] * x[i];
					x11[j_] = 0.5 * (1.0 + std::tanh(s1[j_] / 2.0));
				}

				for (int j_ = 0; j_ < l3; ++j_) {
					s3[j_] = 0.0;
					for (int i = 0; i < l1; ++i)
						s3[j_] += w3[j_][i] * x11[i];
					x3[j_] = s3[j_];
				}

				for (int i = 0; i < l3; ++i)
					sum += std::pow(d0 - x3[0], 2.0);

				s_a[i1] = x3[0];
				eps_test[i1] = d0;
			}

			double rmse = std::sqrt(sum / kp);
			if (amin1 > rmse) amin1 = rmse;
			metric_test = rmse;

			ER2 = 1.0 - sum / Er;
			if (R2_min < ER2)
				R2_min = ER2;
		}

		// ---------------------------- POST-PROCESSING -------------------------------
		for (int i1 = 0; i1 < 383; ++i1) {
			s_eps1[i1] = s_eps[i1] + sr;
			s_a1[i1] = s_a[i1] + sr;
		}

		// Hypothesis a + b*x
		double sum = 0.0;
		for (int i1 = 383 - 2 * kp; i1 < 383 - kp; ++i1)
			sum += s_a1[i1];
		sum /= kp;

		double am2 = 0.0, am3 = 0.0;
		for (int i1 = 383 - 2 * kp; i1 < 383 - kp; ++i1) {
			double r3 = s_a1[i1];
			am2 += std::pow(r3 - sum, 2.0);
			am3 += std::pow(r3 - sum, 3.0);
		}
		am2 /= kp;
		am3 /= kp;

		double gg1 = std::sqrt(am3 * am3 / (am2 * am2 * am2));

		double s_a_av = 0.0, s_eps_av = 0.0;
		for (int i1 = 383 - 2 * kp; i1 < 383 - kp; ++i1) {
			int j = i_sm1[i1];
			s_eps_av += s_eps1[j];
			s_a_av += s_a1[i1];
		}
		s_eps_av /= kp;
		s_a_av /= kp;

		double r = 0.0, r2 = 0.0;
		for (int i1 = 383 - 2 * kp; i1 < 383 - kp; ++i1) {
			int j = i_sm1[i1];
			double r1 = s_eps1[j];
			double r3 = s_a1[i1];
			r += std::pow(r1 - s_eps_av, 2.0);
			r2 += (r1 - s_eps_av) * (r3 - s_a_av);
		}
		double b_t = r2 / r;
		double a_t = s_a_av - b_t * s_eps_av;

		int i3 = 0;
		for (int i1 = 383 - 2 * kp; i1 < 383 - kp; ++i1) {
			int j = i_sm1[i1];
			double r1 = s_eps1[j];
			double r3 = s_a_av + b_t * (r1 - s_eps_av);
			regr[i3++] = r3;
		}

		double sigma = 0.0, sigma_y = 0.0, sum_y = 0.0;
		i3 = 0;
		for (int i1 = 383 - 2 * kp; i1 < 383 - kp; ++i1) {
			int j = i_sm1[i1];
			double r1 = s_eps1[j];
			double r3 = s_a1[i1];
			sigma += std::pow(regr[i3] - r3, 2.0);
			sigma_y += std::pow(r1 - s_eps_av, 2.0);
			sum_y += std::pow(r3 - s_a_av, 2.0);
			i3++;
		}
		sum = sigma;
		sigma /= (kp - 2);

		double fish = (std::pow(b_t - 1.0, 2.0) * sigma_y + kp * std::pow(a_t - s_eps_av * (1.0 - b_t), 2.0)) / (2.0 * sigma);
		if (ip > 35 && fisher1 > fish) {
			fisher1 = fish;
			i_fish1 = ip;
		}

		double denom = std::pow(1.0 - b_t, 2.0) * sigma_y - 2.0 * fish * sigma;
		if (denom != 0.0)
			r = ((1.0 - b_t) * a_t * sigma_y - 2.0 * fish * sigma * s_eps_av) / denom;

		rm = std::max(rm, (fish - dd) * 100.0f / dd);

		rm1 = (fish - dd) * 100.0 / dd;

		ER2 = r = 1.0 - sum / sum_y;
		if (R2_max < r) R2_max = r;

		ip++;
	std::cout << "Final RMSE train: " << amin << "\n";
	std::cout << "Final RMSE test: " << amin1 << "\n";
	std::cout << "Final R2 max: " << R2_max << "\n";
	} while (ip != 160);

    // Pause before exit (replacing multiple getch() calls)
    std::cout << "Press Enter to exit...\n";
    std::cin.get();

    return 0;
}