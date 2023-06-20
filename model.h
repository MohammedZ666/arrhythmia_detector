#ifndef MODEL_H
#define MODEL_H

int8_t LAYER0_KERNEL[]{-13, 69, -123, 32, 99, -32, 17, 22, -7, 8, 8, 72, -30, 41, 47, -44, 12, -1, -24, 10, 14, 38, -10, 31, 6, -49, -3, -13, -34, -12, 17, 5, -29, 20, -33, -36, -16, -14, -20, -14, 19, -14, 20, -5, -17, -31, 0, -9, -26, -27, -9, 10, 10, -11, -24, -17, -15, -1, 7, 10, -23, 8, 8, -21, -8, -14, 17, 1, 6, 12, -29, -1, -31, -8, 0, -7, 32, 3, 15, -2, -3, -9, -31, -7, 12, -7, 34, 3, 11, -31, 16, 5, -29, 9, 2, -4, 22, 5, 7, -16, 34, 12, -1, 0, -12, -7, 23, 7, -4, 3, 12, -2, 9, -6, -18, -5, 15, 10, -2, 10, -1, -13, 14, -5, -5, 7, 10, 6, -1, -2, -1, -4, 0, 1, -7, 16, 2, -3, 8, 5, 3, 8, -8, -2, -2, -1, 8, -5, 7, 19, 5, 0, 1, -2, 28, -24, 14, -1, 2, 17, 7, -7, 5, -9, 42, -26, 15, -2, -2, 3, -6, 0, -12, -3, 2, -7, -2, -5, 14, 7, 2, 4, 2, -3, -6, -1, 0, -2, 6, 3, -1, -1, 2, 2, 2, 5, -4, 0, 1, 2, -3, -2, 1, 3, 7, 2, -6, 2, -5, 0, 0, 1, -3, 2, 3, 2, -16, 4, -2, -2, -1, 2, -2, 2, 0, 0, -14, 2, -2, -3, -2, 1, 0, 2, -2, -1, -7, 1, -4, 0, -2, 0, 3, 0, 2, -3, -2, 0, -4, -1, 0, -1, 1, 0, 4, -2, -1, -3, -1, -3, 2, 1, -2, 1, 3, 0, 4, -4, -1, -4, 3, 2, -4, 0, -1, 0, 5, -3, -2, -1, 4, 0, -2, -1, -1, -2, 6, -2, -5, -1, 2, -3, 0, 1, -3, 3, 3, -2, -1, -2, 1, -2, -1, 1, -8, 8, -1, -3, 3, -1, 0, 0, -1, 2, -10, 8, -1, -2, 4, 2, 4, -2, 3, 0, -1, 4, 5, 1, -7, 1, -1, -4, -2, 2, 3, 3, 2, 1, -4, 2, -1, 0, -5, -2, 0, 0, 0, 1, -1, 2, 14, 1, -2, -3, -4, 1, 1, -1, 3, 2, 4, -10, 5, 4, -2, 1, 17, -6, -2, 4, 20, -13, 8, -3, 5, 3, 14, -7, -3, 7, 36, -8, 8, -12, 9, 6, 0, -6, 0, 1, 25, -3, 2, -1, -5, 9, -14, -3, 8, -1, 0, -1, 5, 0, -12, 4, -11, 2, -2, -5, -22, 1, 7, 0, -10, 0, -28, 7, 1, -4, -36, 4, 8, -1, -5, -2, -30, 7, 2, -4, -34, 8, 4, -7, -1, -5, -32, 8, 2, -4, -18, 12, -2, -14, 3, -2, -32, 12, -4, -3, -9, 15, -14, -10, 5, 14, -37, 14, 0, 1, -15, 19, -18, -5, 6, 27, -35, 12, 3, -2, -19, 20, -7, -12, 9, 26, -21, 10, 4, -6, -11, 18, 7, -23, 11, 12, 5, 11, -2, -7, -8, 22, 4, -4, 1, -3, 16, 11, 0, 1, 0, 28, -2, 10, -7, -23, 39, 9, -6, -7, -3, 31, 6, 11, -5, -41, 47, 10, -6, -10, -36, 17, 22, 8, 17, -55, 61, 4, -15, 7, -74, 15, 4, 15, 34, -46, 54, 5, -4, 36, -71, 5, -4, 0, 61, -29, 47, 27, -8, 26, -49, -21, -8, 4, 66, -4, 11, 68, 12, 52, -43, -58, 1, 6, 103, -25, -2, 73, -14, 37, -36, -40, -7, 5, 104, -39, -22, 69, -21, 41, -48, -30, 5, -13, 127, -54, -18, 48, -38, 25, -61, 16, -10, 16, 110, -66, -70, 53, 15, 15, -18, 63, -7, 18, 93, -111, -80, 91, 23, -3};
int8_t LAYER0_BIAS[]{-3, 1, 0, -1, 1, -1, 0, 3, 1, 1};

int8_t LAYER1_KERNEL[]{1, -53, 16, 8, 1, -9, -14, 36, -26, 19, 18, 13, 0, 24, -47, -30, -17, 13, 12, -4, 1, 5, -16, 17, -9, -9, 15, 36, 6, -23, 20, -38, 15, 8, -29, -70, 25, -2, -31, -12};
int8_t LAYER1_BIAS[]{-8, -1, 17, -9};
float SAMPLE_INPUT_N[] = {0.000000, 0.000000, 0.000000, -0.320000, -0.340000, -0.320000, -0.310000, -0.320000, -0.335000, -0.325000, -0.305000, -0.320000, -0.315000, -0.315000, -0.300000, -0.320000, -0.330000, -0.325000, -0.315000, -0.310000, -0.320000, -0.310000, -0.300000, -0.295000, -0.300000, -0.260000, -0.245000, -0.250000, -0.255000, -0.225000, -0.210000, -0.215000, -0.215000, -0.215000, -0.210000, -0.235000, -0.245000, -0.230000, -0.220000, -0.225000, -0.245000, -0.220000, -0.190000, -0.190000, -0.240000, -0.260000, -0.270000, -0.315000, -0.325000, -0.330000, -0.320000, -0.320000, -0.325000, -0.325000, -0.310000, -0.330000, -0.345000, -0.340000, -0.325000, -0.335000, -0.340000, -0.330000, -0.325000, -0.355000, -0.395000, -0.455000, -0.470000, -0.560000, -0.590000, -0.510000, -0.415000, -0.205000, -0.085000, 0.380000, 0.640000, 0.885000, 0.775000, 0.045000, -0.260000, -0.470000, -0.450000, -0.370000, -0.360000, -0.375000, -0.390000, -0.375000, -0.370000, -0.375000, -0.375000, -0.360000, -0.355000, -0.365000, -0.380000, -0.370000, -0.360000, -0.355000, -0.370000, -0.355000, -0.350000, -0.370000, -0.380000, -0.360000, -0.350000, -0.360000, -0.365000, -0.365000, -0.350000, -0.355000, -0.380000, -0.360000, -0.345000, -0.365000, -0.360000, -0.365000, -0.355000, -0.360000, -0.370000, -0.360000, -0.345000, -0.360000, -0.350000, -0.340000, -0.330000, -0.335000, -0.345000, -0.335000, -0.340000, -0.355000, -0.360000, -0.355000, -0.350000, -0.360000, -0.380000, -0.370000, -0.380000, -0.390000, -0.400000, -0.405000, -0.380000, -0.380000, -0.385000, -0.380000, -0.355000, -0.370000, -0.360000, -0.340000, -0.325000, 0.000000, 0.000000, 0.000000};
float SAMPLE_INPUT_S[] = {0.000000, 0.000000, 0.000000, -0.265000, -0.270000, -0.255000, -0.255000, -0.245000, -0.215000, -0.185000, -0.180000, -0.160000, -0.155000, -0.160000, -0.150000, -0.165000, -0.175000, -0.215000, -0.245000, -0.305000, -0.335000, -0.375000, -0.390000, -0.415000, -0.430000, -0.455000, -0.470000, -0.510000, -0.520000, -0.495000, -0.485000, -0.435000, -0.415000, -0.390000, -0.395000, -0.395000, -0.400000, -0.375000, -0.360000, -0.350000, -0.355000, -0.340000, -0.350000, -0.375000, -0.385000, -0.415000, -0.425000, -0.420000, -0.410000, -0.410000, -0.420000, -0.450000, -0.460000, -0.475000, -0.475000, -0.475000, -0.465000, -0.465000, -0.460000, -0.480000, -0.495000, -0.495000, -0.505000, -0.495000, -0.490000, -0.480000, -0.485000, -0.510000, -0.515000, -0.530000, -0.535000, -0.345000, -0.110000, 0.690000, 1.130000, 1.555000, 1.405000, 0.305000, -0.335000, -1.615000, -1.675000, -1.095000, -0.805000, -0.575000, -0.570000, -0.585000, -0.575000, -0.565000, -0.555000, -0.560000, -0.565000, -0.585000, -0.585000, -0.605000, -0.590000, -0.560000, -0.545000, -0.525000, -0.535000, -0.540000, -0.550000, -0.540000, -0.540000, -0.520000, -0.525000, -0.495000, -0.485000, -0.480000, -0.490000, -0.490000, -0.465000, -0.445000, -0.425000, -0.395000, -0.375000, -0.385000, -0.370000, -0.300000, -0.270000, -0.230000, -0.210000, -0.180000, -0.155000, -0.140000, -0.140000, -0.150000, -0.150000, -0.160000, -0.180000, -0.175000, -0.195000, -0.205000, -0.205000, -0.220000, -0.225000, -0.220000, -0.215000, -0.220000, -0.230000, -0.250000, -0.275000, -0.330000, -0.360000, -0.410000, -0.425000, -0.465000, -0.505000, 0.000000, 0.000000, 0.000000};
float SAMPLE_INPUT_V[] = {0.000000, 0.000000, 0.000000, 0.400000, 0.400000, 0.380000, 0.375000, 0.330000, 0.320000, 0.300000, 0.265000, 0.240000, 0.215000, 0.205000, 0.185000, 0.150000, 0.140000, 0.145000, 0.125000, 0.090000, 0.095000, 0.085000, 0.085000, 0.050000, 0.045000, 0.060000, 0.050000, 0.020000, 0.045000, 0.060000, 0.055000, 0.040000, 0.060000, 0.070000, 0.080000, 0.065000, 0.085000, 0.110000, 0.105000, 0.110000, 0.125000, 0.140000, 0.135000, 0.135000, 0.155000, 0.170000, 0.150000, 0.150000, 0.165000, 0.175000, 0.175000, 0.160000, 0.155000, 0.170000, 0.155000, 0.145000, 0.160000, 0.175000, 0.180000, 0.175000, 0.185000, 0.210000, 0.200000, 0.200000, 0.215000, 0.255000, 0.255000, 0.255000, 0.240000, 0.145000, 0.095000, -0.005000, -0.040000, -0.140000, -0.185000, -0.210000, -0.190000, -0.115000, -0.075000, 0.000000, 0.055000, 0.135000, 0.115000, 0.045000, 0.005000, -0.095000, -0.150000, -0.195000, -0.200000, -0.225000, -0.230000, -0.255000, -0.260000, -0.225000, -0.250000, -0.255000, -0.255000, -0.250000, -0.275000, -0.275000, -0.280000, -0.270000, -0.285000, -0.300000, -0.290000, -0.280000, -0.290000, -0.295000, -0.285000, -0.275000, -0.295000, -0.285000, -0.275000, -0.260000, -0.260000, -0.255000, -0.235000, -0.210000, -0.210000, -0.200000, -0.170000, -0.140000, -0.135000, -0.120000, -0.100000, -0.050000, -0.025000, -0.005000, 0.025000, 0.070000, 0.090000, 0.130000, 0.165000, 0.215000, 0.225000, 0.245000, 0.285000, 0.310000, 0.330000, 0.335000, 0.335000, 0.360000, 0.355000, 0.355000, 0.360000, 0.360000, 0.345000, 0.000000, 0.000000, 0.000000};
float SAMPLE_INPUT_F[] = {0.000000, 0.000000, 0.000000, 0.050000, 0.055000, 0.075000, 0.095000, 0.070000, 0.060000, 0.020000, 0.030000, 0.015000, -0.020000, -0.100000, -0.105000, -0.090000, -0.095000, -0.180000, -0.180000, -0.185000, -0.255000, -0.315000, -0.260000, -0.195000, -0.235000, -0.265000, -0.235000, -0.240000, -0.265000, -0.250000, -0.200000, -0.125000, -0.150000, -0.240000, -0.235000, -0.155000, -0.170000, -0.220000, -0.200000, -0.110000, -0.115000, -0.140000, -0.125000, -0.085000, -0.100000, -0.060000, -0.035000, -0.015000, -0.030000, 0.005000, -0.020000, -0.075000, -0.070000, -0.085000, -0.100000, -0.120000, -0.130000, -0.155000, -0.140000, -0.140000, -0.145000, -0.120000, -0.080000, -0.005000, 0.040000, 0.070000, 0.120000, 0.240000, 0.265000, 0.350000, 0.450000, 0.750000, 0.940000, 1.520000, 1.810000, 2.235000, 2.160000, 1.330000, 0.815000, 0.090000, -0.100000, -0.300000, -0.350000, -0.400000, -0.430000, -0.510000, -0.535000, -0.550000, -0.575000, -0.685000, -0.695000, -0.655000, -0.685000, -0.680000, -0.655000, -0.685000, -0.715000, -0.790000, -0.810000, -0.815000, -0.805000, -0.790000, -0.785000, -0.805000, -0.810000, -0.770000, -0.775000, -0.760000, -0.735000, -0.710000, -0.695000, -0.665000, -0.670000, -0.625000, -0.585000, -0.540000, -0.540000, -0.515000, -0.455000, -0.415000, -0.410000, -0.370000, -0.330000, -0.285000, -0.265000, -0.255000, -0.230000, -0.195000, -0.190000, -0.190000, -0.165000, -0.155000, -0.180000, -0.195000, -0.155000, -0.095000, -0.130000, -0.145000, -0.105000, -0.090000, -0.135000, -0.145000, -0.100000, -0.080000, -0.070000, -0.065000, -0.080000, 0.000000, 0.000000, 0.000000};

float dequantize(int8_t x_q)
{
    return 0.817042f * x_q;
}

#endif