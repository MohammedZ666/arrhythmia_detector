#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include "model.h"

int main()

{
    float ecg_data[150 + 15 - 1]{0.0};
    int i = detect_qrs(SAMPLE_INPUT_F, ecg_data, 150);
    printf("%d ", i);
    float layer2in[SAMPLE_IN_LEN0 * LAYER0_KERNEL_DIM1]{0.0};
    dense(layer2in, &ecg_data[i - 30], LAYER0_KERNEL, LAYER0_BIAS, SAMPLE_IN_LEN0, SAMPLE_IN_LEN1, LAYER0_KERNEL_DIM0, LAYER0_KERNEL_DIM1, 'r');
    float output[SAMPLE_IN_LEN0 * LAYER1_KERNEL_DIM1]{0.0};
    dense(output, layer2in, LAYER1_KERNEL, LAYER1_BIAS, SAMPLE_IN_LEN0, LAYER0_KERNEL_DIM1, LAYER1_KERNEL_DIM0, LAYER1_KERNEL_DIM1, 't');

    for (int i = 0; i < SAMPLE_IN_LEN0 * LAYER1_KERNEL_DIM1; i++)
    {
        printf("%f ", output[i]);
    }
    printf("\n");
}