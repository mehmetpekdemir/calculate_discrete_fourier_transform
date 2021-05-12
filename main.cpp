#include <iostream>
#include "AudioFile.h"
#include <math.h>
using std::vector;
#define PI 3.14159265

AudioFile<double> load_audio_file();
void print_audio_file_summary(AudioFile<double> audioFile);
vector<double> get_signal(AudioFile<double> audioFile);
int get_num_samples_per_channel(AudioFile<double> audioFile);
void calculate_discrete_fourier_transform();

int main()
{

    /*
        load audio file
    */
    AudioFile<double> auidoFile = load_audio_file();

    /*
        print audio file summary
    */
    print_audio_file_summary(auidoFile);

    /*
        get signal
    */
    vector<double> signal = get_signal(auidoFile);

    /*
        get num samples per channel
    */
    int num_samples_per_channel = get_num_samples_per_channel(auidoFile);

    /*
        calculate discrete fourier transform
    */
    calculate_discrete_fourier_transform(num_samples_per_channel, signal);

    return 0;
}

AudioFile<double> load_audio_file()
{
    AudioFile<double> audioFile;
    audioFile.load("/home/mehmet/Workspace/DFT/track.wav"); // try catch
}

void print_audio_file_summary(AudioFile<double> audioFile)
{
    audioFile.printSummary();
}

vector<double> get_signal(AudioFile<double> audioFile)
{
    const int channel = 0;
    vector<double> signal = audioFile.samples[channel];
    return signal;
}

int get_num_samples_per_channel(AudioFile<double> audioFile)
{
    int num_samples_per_channel = audioFile.getNumSamplesPerChannel();
    return num_samples_per_channel;
}

void calculate_discrete_fourier_transform(int num_samples_per_channel, vector<double> signal)
{
    const int half_of_num_samples_per_channel = num_samples_per_channel / 2;
    double *real = new double[half_of_num_samples_per_channel];
    double *imag = new double[half_of_num_samples_per_channel];

    for (int i = 0; i < half_of_num_samples_per_channel; i++)
    {
        double sumReal = 0;
        double sumImaginary = 0;
        for (int j = 0; j < num_samples_per_channel; j++)
        {
            sumReal += signal[j] * cos((2 * PI * i * j) / num_samples_per_channel);
            sumImaginary += signal[j] * sin((2 * PI * i * j) / num_samples_per_channel);
        }
        real[i] = sumReal;
        imag[i] = -sumImaginary;
        std::cout << real[i] << "\t" << imag[i] << std::endl;
    }
}