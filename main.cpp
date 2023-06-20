#include <iostream>
#include <fstream>
#include <vector>

using namespace std;
float last_qrs_index{-1};
int refractory_period{120}; // Change proportionally when adjusting frequency (in samples).
float threshold{0.0};
int integration_window{15};
vector<float> load_data_from_file(const char *path)
{
    std::ifstream file(path);
    vector<float> ecg_data;
    if (file.is_open())
    {
        string line;
        bool is_first_line = true;
        while (std::getline(file, line))
        {
            if (is_first_line)
                is_first_line = false;
            else
            {
                size_t c_index = line.find(",", 0);
                float ecg_value = stod(line.substr(c_index + 1, line.length()));
                ecg_data.push_back(ecg_value);
            }
        }
        file.close();
    }
    return ecg_data;
}

void detect_qrs(float *x, int len)
{
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
        // IIR difference equation
        y[i] = (b0 * x[i] + b1 * x_1 + b2 * x_2 - (a1 * y_1 + a2 * y_2)) / a0;

        // shift delayed x, y samples
        x_2 = x_1;
        x_1 = x[i];
        y_2 = y_1;
        y_1 = y[i];
    }

    // See QRSDetectorOffline.py after filtering
    for (int i = 0; i < 5; i++)
        y[i] = y[5];

    // Differentiation
    for (int i = 0; i + 1 < len; i++)
    {
        y[i] = y[i + 1] - y[i];
        y[i] *= y[i];
    }

    len--; // Reducing length since diff causes 1 element reduction

    int n = len + integration_window - 1;
    float res[n]{0.0};
    int max_len = len > integration_window ? len : integration_window;
    for (int i = 0; i < n; i++)
    {

        int kMax = i < max_len ? i : max_len;
        for (int k = 0; k <= kMax; k++)
        {
            if (k < len && (i - k) < integration_window)
                res[i] += y[k];
        }
    }

    int spacing{50};
    float limit{0.35};
    vector<int> peak_ind;

    for (int i = 0; i < n; i++)
    {
        int lookup = (i - 1) - spacing < 0 ? i : spacing;
        bool is_broken = false;
        for (int j = 1; j <= lookup; j++)
        {
            if (!(res[i] > res[i - j]))
            {
                is_broken = true;
                break;
            }
        }
        if (is_broken)
            continue;

        lookup = (i + 1) + spacing > n - 1 ? (n - 1) - i : spacing;
        for (int j = 1; j <= lookup; j++)
        {
            if (!(res[i] > res[i + j]))
            {
                is_broken = true;
                break;
            }
        }
        if (!is_broken && res[i] > limit)
            peak_ind.push_back(i);
    }
    vector<int> qrs_indices, noise_indices;
    float qrs_peak_est{0.0};
    float noise_peak_est{0.0};
    float noise_peak_val{0.0};
    float qrs_filt{0.125};
    float noise_filt{0.125};
    float qrs_noise_diff_weight{0.25};

    for (int i = 0; i < peak_ind.size(); i++)
    {
        if (last_qrs_index == -1 || i - last_qrs_index > refractory_period)
        {
            float current_peak_val = res[peak_ind[i]];
            if (current_peak_val > threshold)
            {
                qrs_indices.push_back(peak_ind[i]);
                qrs_peak_est = qrs_filt * current_peak_val + (1 - qrs_filt) * qrs_peak_est;
            }
            else
            {
                noise_indices.push_back(peak_ind[i]);
                noise_peak_est = noise_filt * current_peak_val + (1 - noise_filt) * noise_peak_est;
            }
            threshold = noise_peak_est + qrs_noise_diff_weight * (qrs_peak_est - noise_peak_est);
        }
    }
    for (int i = 0; i < qrs_indices.size(); i++)
    {
        printf("%d, ", qrs_indices[i]);
    }
    for (int i = 0; i < noise_indices.size(); i++)
    {
        printf("%d, ", noise_indices[i]);
    }
    printf("\n");
}

void print_output(int *qrs_indices, int *noise_indices, int qrs_len, int noise_len)
{
    for (int i = 0; i < qrs_len; i++)
    {
        printf("%d, ", qrs_indices[i]);
    }
    for (int i = 0; i < noise_len; i++)
    {
        printf("%d, ", noise_indices[i]);
    }
    printf("\n");
}

void detect_qrs_a(float *x, float *res, int len, int *qrs_indices, int *noise_indices, int *qcount, int *ncount)
{
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
        // IIR difference equation
        y[i] = (b0 * x[i] + b1 * x_1 + b2 * x_2 - (a1 * y_1 + a2 * y_2)) / a0;

        // shift delayed x, y samples
        x_2 = x_1;
        x_1 = x[i];
        y_2 = y_1;
        y_1 = y[i];
    }

    // See QRSDetectorOffline.py after filtering
    for (int i = 0; i < 5; i++)
        y[i] = y[5];

    // Differentiation
    for (int i = 0; i + 1 < len; i++)
    {
        y[i] = y[i + 1] - y[i];
        y[i] *= y[i];
    }

    len--; // Reducing length since diff causes 1 element reduction

    int n = len + integration_window - 1;
    // float res[n]{0.0};
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

    int spacing{50};
    float limit{0.35};
    // int qcount = 0;
    // int ncount = 0;
    float qrs_peak_est{0.0};
    float noise_peak_est{0.0};
    float noise_peak_val{0.0};
    float qrs_filt{0.125};
    float noise_filt{0.125};
    float qrs_noise_diff_weight{0.25};

    for (int i = 0; i < n; i++)
    {
        int lookup = (i - 1) - spacing < 0 ? i : spacing;
        bool is_broken = false;
        for (int j = 1; j <= lookup; j++)
        {
            if (!(res[i] > res[i - j]))
            {
                is_broken = true;
                break;
            }
        }
        if (is_broken)
            continue;

        lookup = (i + 1) + spacing > n - 1 ? (n - 1) - i : spacing;
        for (int j = 1; j <= lookup; j++)
        {
            if (!(res[i] > res[i + j]))
            {
                is_broken = true;
                break;
            }
        }
        if (!is_broken && res[i] > limit)
        {
            if (last_qrs_index == -1 || i - last_qrs_index > refractory_period)
            {
                float current_peak_val = res[i];
                if (current_peak_val > threshold)
                {
                    last_qrs_index = i;
                    qrs_indices[*qcount] = i;
                    *qcount = *qcount + 1;
                    qrs_peak_est = qrs_filt * current_peak_val + (1 - qrs_filt) * qrs_peak_est;
                }
                else
                {
                    noise_indices[*ncount] = i;
                    *ncount = *ncount + 1;
                    noise_peak_est = noise_filt * current_peak_val + (1 - noise_filt) * noise_peak_est;
                }
                threshold = noise_peak_est + qrs_noise_diff_weight * (qrs_peak_est - noise_peak_est);
            }
        }
    }

    // print_output(qrs_indices, noise_indices, *qcount, *ncount);
}
int main()
{
    vector<float> ecg_data_raw = load_data_from_file("ecg_data/ecg_data_2.csv");
    float x[ecg_data_raw.size()];
    copy(ecg_data_raw.begin(), ecg_data_raw.end(), x);
    int qrs_indices[ecg_data_raw.size() + integration_window - 2]{0};
    int noise_indices[ecg_data_raw.size() + integration_window - 2]{0};
    float res[ecg_data_raw.size() + integration_window - 2]{0.0};

    int *qrs_len = NULL;
    int *noise_len = NULL;

    qrs_len = (int *)malloc(sizeof(int));
    noise_len = (int *)malloc(sizeof(int));

    *qrs_len = 0;
    *noise_len = 0;

    // detect_qrs(x, ecg_data_raw.size());
    detect_qrs_a(x, res, ecg_data_raw.size(), qrs_indices, noise_indices, qrs_len, noise_len);
    // print_output(qrs_indices, noise_indices, *qrs_len, *noise_len);

    // not freeing since program is ending
    // free(qrs_len);
    // free(noise_len);
    int count = 0;
    for (int i = qrs_indices[0] - 30; i <= qrs_indices[0] + 30; i++, count++)
    {
        printf("%f\n", res[i]);
    }
}