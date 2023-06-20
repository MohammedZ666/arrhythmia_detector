#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include "model.h"
#include <sstream>
#include <sciplot/sciplot.hpp>
#include <valarray>

using namespace std;
using namespace sciplot;

void count_dim(const char *path, int *shape)
{
    int i, j = 0;
    std::ifstream file(path);
    if (file.is_open())
    {
        string line;
        bool is_first_line = true;
        while (std::getline(file, line))
        {
            if (!line.empty())
            {

                string sfloat;
                std::stringstream stream(line);
                if (is_first_line)
                {
                    is_first_line = false;
                    while (std::getline(stream, sfloat, ','))
                    {
                        shape[1]++;
                    }
                }
                shape[0]++;
            }
        }
        file.close();
    }
}
vector<vector<float>> load_testing_data(const char *path)
{
    std::ifstream file(path);
    int shape[2]{-1, -1};
    count_dim(path, shape);
    vector<vector<float>> test_data(shape[0]);
    if (file.is_open())
    {
        string line;
        bool is_first_row = true;
        int i = 0;
        while (std::getline(file, line))
        {
            if (is_first_row)
                is_first_row = false;
            else
            {
                int j = 0;
                bool is_first_col = true;
                string sfloat;
                std::stringstream stream(line);
                test_data[i] = vector<float>(shape[1]);
                while (std::getline(stream, sfloat, ','))
                {
                    if (is_first_col)
                        is_first_col = false;
                    else
                    {
                        test_data[i][j] = stof(sfloat);
                        j++;
                    }
                }

                i++;
            }
        }
        file.close();
    }
    return test_data;
}

void find_peaks(float *data, int data_len, int spacing, float limit, int *output)
{
    for (int i = 0; i < data_len; i++)
    {
        int lookup = (i - 1) - spacing < 0 ? i : spacing;
        bool is_broken = false;
        for (int j = 1; j <= lookup; j++)
        {
            if (!(data[i] > data[i - j]))
            {
                is_broken = true;
                break;
            }
        }
        if (is_broken)
            continue;

        lookup = (i + 1) + spacing > data_len - 1 ? (data_len - 1) - i : spacing;
        for (int j = 1; j <= lookup; j++)
        {
            if (!(data[i] > data[i + j]))
            {
                is_broken = true;
                break;
            }
        }
        if (!is_broken && data[i] > limit)
            output[i] = 1;
    }
}

void show_progress(float progress)
{
    if (progress <= 1.0)
    {
        int barWidth = 70;

        std::cout << "[";
        int pos = barWidth * progress;
        for (int i = 0; i < barWidth; ++i)
        {
            if (i < pos)
                std::cout << "=";
            else if (i == pos)
                std::cout << ">";
            else
                std::cout << " ";
        }
        std::cout << "] " << int(progress * 100.0) << " %\r";
    }
    if (progress == 1.0)
        std::cout << endl;
}

int detect_qrs(float *x, float *res)
{
    int len = 150;
    float y[len]{0.0};
    float x_2 = 0.0; // delayed x, y samples
    float x_1 = 0.0;
    float y_1 = 0.0;
    float y_2 = 0.0;

    const float a0 = 1.0;
    const float a1 = -1.73356294;
    const float a2 = 0.77567951;
    const float b0 = 0.11216024;
    const float b1 = 0;
    const float b2 = -0.11216024;

    for (int i = 0; i < len; ++i)
    {
        y[i] = (b0 * x[i] + b1 * x_1 + b2 * x_2 - (a1 * y_1 + a2 * y_2)) / a0;

        x_2 = x_1;
        x_1 = x[i];
        y_2 = y_1;
        y_1 = y[i];
    }

    for (int i = 0; i < 5; i++)
        y[i] = y[5];

    for (int i = 0; i + 1 < len; i++)
    {
        y[i] = y[i + 1] - y[i];
        y[i] *= y[i];
    }

    len--;

    int integration_window{15};
    int n = len + integration_window - 1;
    int max_len = len > integration_window ? len : integration_window;
    for (int i = 0; i < n; i++)
    {

        int kMax = i < max_len ? i : max_len;
        for (int k = 0; k <= kMax; k++)
        {
            if (k < len && (i - k) < integration_window)
            {
                res[i] += y[k];
            }
        }
    }

    int spacing{30};
    float limit{0.35};
    float qrs_peak_est{0.0};
    float noise_peak_est{0.0};
    float noise_peak_val{0.0};
    float qrs_filt{0.125};
    float noise_filt{0.125};
    float qrs_noise_diff_weight{0.25};
    float last_qrs_index{-1};
    uint8_t refractory_period{120};
    float threshold{0.0};
    uint8_t i = 60;
    uint8_t j = 1;
    uint8_t max_ind = i;
    bool is_peak = false;

    while (j <= spacing && i < len - 30)
    {
        if (!(res[i] > res[i - j] && res[i] > res[i + j]))
        {
            i = i + j;
            j = 1;
            continue;
        }

        else if (j == spacing)
        {
            is_peak = true;
            max_ind = res[i] > res[max_ind] ? i : max_ind;
            i += spacing + 1;
        }
        j++;
    }
    if (!is_peak)
        return -1;

    i = max_ind;
    if (last_qrs_index == -1 || i - last_qrs_index > refractory_period)
    {
        float current_peak_val = res[i];
        if (current_peak_val > threshold)
        {
            qrs_peak_est = qrs_filt * current_peak_val + (1 - qrs_filt) * qrs_peak_est;
            last_qrs_index = i;
            return i;
        }
        else
        {
            noise_peak_est = noise_filt * current_peak_val + (1 - noise_filt) * noise_peak_est;
        }
        threshold = noise_peak_est + qrs_noise_diff_weight * (qrs_peak_est - noise_peak_est);
    }
    return -1;
}

void softmax(float *x, int rows, int cols)
{
    float sum = 0;
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            x[i * cols + j] = expf(x[i * cols + j]);
            sum += x[i * cols + j];
        }
    }
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            x[i * cols + j] = x[i * cols + j] / sum;
        }
    }
}
void relu(float *x, int rows, int cols)
{
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            x[i * cols + j] = x[i * cols + j] < 0.0f ? 0.0f : x[i * cols + j];
        }
    }
}

void sigmoid(float *x, int rows, int cols)
{
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            x[i * cols + j] = 1 / (1 + expf(-1 * x[i * cols + j]));
        }
    }
}

void matmul(float *result, float *input, int8_t *kernel, uint8_t r1, uint8_t c1, uint8_t r2, uint8_t c2)
{

    for (int i = 0; i < r1; i++)
    {
        for (int j = 0; j < c2; j++)
        {
            for (int k = 0; k < c1; k++)
            {
                result[i * c1 + j] += input[i * c1 + k] * dequantize(kernel[k * c2 + j]);
            }
        }
    }
}

void dense(float *output, float *input, int8_t *kernel, int8_t *bias, uint8_t ir, uint8_t ic, uint8_t kr, uint8_t kc, char gate)
{
    matmul(output, input, kernel, ir, ic, kr, kc);
    for (int i = 0; i < ir; i++)
    {
        for (int j = 0; j < kc; j++)
        {
            output[i * kc + j] += dequantize(bias[i * kc + j]);
        }
    }

    if (gate == 's')
        sigmoid(output, ir, kc);
    else if (gate == 'r')
        relu(output, ir, kc);
    else
        softmax(output, ir, kc);
}

int argmax(float *output, int len)
{
    int max_ind = 0;
    float max = output[max_ind];

    for (int i = 1; i < len; i++)
    {
        if (output[i] > max)
        {
            max = output[i];
            max_ind = i;
        }
    }
    return max_ind;
}
void run_test()
{
    vector<vector<float>> x_test = load_testing_data("ecg_data/x_test.csv");
    vector<vector<float>> y_test = load_testing_data("ecg_data/y_test.csv");
    int corr = 0;
    for (int i = 0; i < x_test.size(); i++)
    {
        float *x_testf = &x_test[i][0];

        float layer2in[10]{};
        dense(layer2in, x_testf, LAYER0_KERNEL, LAYER0_BIAS, 1, 61, 61, 10, 's');
        float output[4]{};
        dense(output, layer2in, LAYER1_KERNEL, LAYER1_BIAS, 1, 10, 10, 4, 't');
        float *y_testf = &y_test[i][0];
        if (argmax(output, 4) == argmax(y_testf, y_test[i].size()))
            corr++;

        if (i % 50 == 0 || i == (x_test.size() - 1))
            show_progress((i + 1) * 1.0 / x_test.size());
    }

    printf("\nTesting accuracy is %.60g %%\n", (corr * 100.0 / x_test.size()));
}
void draw_beat(float *beat, int len)
{
    Vec x = linspace(0, len, len);
    // Create a Plot object
    Plot2D plot;

    // Set the x and y labels
    plot.xlabel("x");
    plot.ylabel("y");

    std::valarray<float> val(beat, len);

    // Set the x and y ranges
    plot.xrange(0.0, len);
    plot.yrange(0.0, val.max());

    // Set the legend to be on the bottom along the horizontal
    plot.legend()
        .atOutsideBottom()
        .displayHorizontal()
        .displayExpandWidthBy(2);

    // Plot sin(i*x) from i = 1 to i = 6

    Vec y = Vec(len);

    for (int i = 0; i < len; i++)
        y[i] = (double)beat[i];

    plot.drawCurve(x, y);

    // Create figure to hold plot
    Figure fig = {{plot}};
    // Create canvas to hold figure
    Canvas canvas = {{fig}};

    // Show the plot in a pop-up window
    canvas.show();
}
void test_sample()
{
    float res[149 + 15 - 1]{0.0};
    float layer2in[10]{0.0};
    float output[4]{0.0};

    int i = detect_qrs(SAMPLE_INPUT_F, res);

    // draw_beat(res, 149 + 15 - 1);

    dense(layer2in, &res[60], LAYER0_KERNEL, LAYER0_BIAS, 1, 61, 61, 10, 's');
    dense(output, layer2in, LAYER1_KERNEL, LAYER1_BIAS, 1, 10, 10, 4, 't');
    for (int i = 0; i < 4; i++)
    {
        printf("%f ", output[i]);
    }
    printf("\n");
}
int main()
{
    run_test();
    // test_sample();
}
