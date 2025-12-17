#include <iostream>
#include "AudioFile.h"
#include <map>
#define _USE_MATH_DEFINES
#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

struct MidSide {
    AudioFile<double> mid;
    AudioFile<double> side;
};

double w(int n, int N)
{
    return std::sin(M_PI / (2.0 * N) * (n + 0.5));
}

MidSide splitToMidAndSide(const AudioFile<double>& audio) {
    MidSide ms;

    if (audio.getNumChannels() != 2) {
        throw std::runtime_error("Audio ni stereo!");
    }

    int numSamplesPerChannel = audio.getNumSamplesPerChannel();
    
    ms.mid.setNumChannels(1);
    ms.side.setNumChannels(1);
    ms.mid.setNumSamplesPerChannel(numSamplesPerChannel);
    ms.side.setNumSamplesPerChannel(numSamplesPerChannel);
    ms.mid.setSampleRate(audio.getSampleRate());
    ms.side.setSampleRate(audio.getSampleRate());

    for (int i = 0; i < numSamplesPerChannel; ++i) {
        double L = audio.samples[0][i];
        double R = audio.samples[1][i];

        ms.mid.samples[0][i] = (L + R) / 2;
        ms.side.samples[0][i] = (L - R) / 2;
    }

    return ms;
}

void buildBlocks(int N, std::vector<std::vector<double>>& blocks) {
    int blockSize = N;
    int half = N / 2;
    int numBlocks = N - 1;

    blocks.resize(numBlocks, std::vector<double>(blockSize));

    int globalValue = 0;

    for (int b = 0; b < numBlocks; ++b) {
        for (int j = 0; j < blockSize; ++j) {
            if (b == 0 && j < half) {
                blocks[b][j] = 0;
            }
            else if (b == numBlocks - 1 && j >= half) {
                blocks[b][j] = 0;
            }
            else {
                blocks[b][j] = ++globalValue;
            }
        }

        if (b < numBlocks - 1) {
            globalValue -= half;
        }
    }
}

void windowFunction(std::vector<std::vector<double>>& blocks) {
    for (int b = 0; b < blocks.size() - 1; ++b) {
        for (int j = 0; j < blocks.size(); ++j) {
            blocks[b][j] *= w(j, blocks[b].size());
        }
    }
}

void MDCT(const std::vector<double>& block) {
    int N = block.size() / 2;
}

void compressAudio(const AudioFile<double>& audio, int M, int N) {
    MidSide ms = splitToMidAndSide(audio);
    AudioFile<double> Mid = ms.mid;
    AudioFile<double> Side = ms.side;

    std::vector<std::vector<double>> blocks;
    buildBlocks(N, blocks);
    windowFunction(blocks);


}

void dekompresijaZvoka() {
    std::cout << "Izbrali ste dekompresijo zvoka!";
}

void predvajajDekompresiranZvok() {
    std::cout << "Izbrali ste predvajanje dekompresiranjega zvoka!";
}

int main()
{
    AudioFile<double> audio;
    int userChoice = 0, M, N; //M = Faktor kompresije    - N = BLOCK SIZE

    while (true) {
        std::cout << "Izberite nacin:\n\n"
            << "1 - Kompresija zvoka\n"
            << "2 - Izvedi dekompresijo\n"
            << "3 - Predvajaj dekompresiran zvok\n"
            << "0 - IZHOD\n\n";

        std::cin >> userChoice;

        switch (userChoice)
        {
        case 1: {
            std::cout << "Izberite faktor stiskanja M:\n";
            std::cin >> M;
            std::cout << "Izberite velikost bloka N:\n";
            std::cin >> N;
            std::cout << "\n\n";

            compressAudio(audio, M, N * 2);

            break;
        }
        case 2: {
            dekompresijaZvoka();
            break;
        }
        case 3: {
            predvajajDekompresiranZvok();
            break;
        }
        case 0: {
            std::cout << "\nIzhod iz programa.\n";
            return 0;
        }
        default:
            std::cout << "Neveljavna izbira!\n";
        }
    }


    if (!audio.load("test.wav"))
    {
        std::cout << "FAIL\n";
        return 1;
    }

    audio.printSummary();
    return 0;
}
