#include <iostream>
#include "AudioFile.h"
#include <vector>
#include <cmath>
#include <iomanip>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

struct MidSide {
    AudioFile<double> mid;
    AudioFile<double> side;
};

double w(int n, int N) {
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

        ms.mid.samples[0][i] = (L + R) / 2.0;
        ms.side.samples[0][i] = (L - R) / 2.0;
    }

    return ms;
}

void buildBlocks(const std::vector<double>& signal, int N, std::vector<std::vector<double>>& blocks) {
    int blockSize = 2 * N;
    int half = N;
    int numSamples = signal.size();

    std::vector<double> paddedSignal;
    paddedSignal.insert(paddedSignal.end(), N, 0.0);
    paddedSignal.insert(paddedSignal.end(), signal.begin(), signal.end());
    paddedSignal.insert(paddedSignal.end(), N, 0.0);

    int numBlocks = (paddedSignal.size() - N) / N;
    blocks.clear();
    blocks.resize(numBlocks, std::vector<double>(blockSize, 0.0));

    for (int b = 0; b < numBlocks; ++b) {
        int startIdx = b * N;
        for (int j = 0; j < blockSize && (startIdx + j) < paddedSignal.size(); ++j) {
            blocks[b][j] = paddedSignal[startIdx + j];
        }
    }
}

void applyWindow(std::vector<std::vector<double>>& blocks) {
    for (int b = 0; b < blocks.size(); ++b) {
        int blockSize = blocks[b].size();
        for (int j = 0; j < blockSize; ++j) {
            blocks[b][j] *= w(j, blockSize);
        }
    }
}

std::vector<double> MDCT(const std::vector<double>& block) {
    int twoN = block.size();
    int N = twoN / 2;
    std::vector<double> mdctBlock(N, 0.0);

    for (int k = 0; k < N; ++k) {
        double sum = 0.0;
        for (int n = 0; n < twoN; ++n) {
            sum += block[n] * std::cos(M_PI / N * (n + 0.5 + N / 2.0) * (k + 0.5));
        }
        mdctBlock[k] = sum;
    }

    return mdctBlock;
}

std::vector<std::vector<double>> MDCTAllBlocks(const std::vector<std::vector<double>>& blocks) {
    std::vector<std::vector<double>> mdctBlocks;
    mdctBlocks.reserve(blocks.size());

    for (const auto& block : blocks) {
        mdctBlocks.push_back(MDCT(block));
    }

    return mdctBlocks;
}

void quantizeAndCompress(std::vector<std::vector<double>>& mdctBlocks, int M) {
    for (auto& block : mdctBlocks) {
        int N = block.size();

        for (int i = 0; i < N; ++i) {
            block[i] = std::round(block[i]);
        }

        for (int i = N - M; i < N; ++i) {
            block[i] = 0.0;
        }
    }
}

void compressAudio(const AudioFile<double>& audio, int M, int N) {
    std::cout << "=== KOMPRESIJA ZVOKA ===\n";
    std::cout << "N = " << N << ", M = " << M << "\n";
    std::cout << "Velikost bloka: " << 2 * N << "\n\n";

    MidSide ms = splitToMidAndSide(audio);

    // Procesiranje Mid kanala
    std::vector<std::vector<double>> midBlocks;
    buildBlocks(ms.mid.samples[0], N, midBlocks);
    std::cout << "Mid kanal: " << midBlocks.size() << " blokov\n";

    applyWindow(midBlocks);
    std::vector<std::vector<double>> midMDCT = MDCTAllBlocks(midBlocks);
    quantizeAndCompress(midMDCT, M);

    // Procesiranje Side kanala
    std::vector<std::vector<double>> sideBlocks;
    buildBlocks(ms.side.samples[0], N, sideBlocks);
    std::cout << "Side kanal: " << sideBlocks.size() << " blokov\n";

    applyWindow(sideBlocks);
    std::vector<std::vector<double>> sideMDCT = MDCTAllBlocks(sideBlocks);
    quantizeAndCompress(sideMDCT, M);

    std::cout << "\nKompresija zakljucena!\n";

    // TODO: Zapis v binarno datoteko
}

void dekompresijaZvoka() {
    std::cout << "Izbrali ste dekompresijo zvoka!\n";
}

void predvajajDekompresiranZvok() {
    std::cout << "Izbrali ste predvajanje dekompresiranjega zvoka!\n";
}

int main() {
    AudioFile<double> audio;

    std::cout << "Nalagam test.wav...\n";
    if (!audio.load("test.wav")) {
        std::cout << "NAPAKA: Ne morem naloziti test.wav\n";
        std::cout << "Program bo deloval v demo nacinu brez audio datoteke.\n\n";
    }
    else {
        std::cout << "Audio datoteka nalozena uspesno!\n";
        audio.printSummary();
        std::cout << "\n";
    }

    int userChoice = 0, M, N;

    while (true) {
        std::cout << "\nIzberite nacin:\n\n"
            << "1 - Kompresija zvoka\n"
            << "2 - Izvedi dekompresijo\n"
            << "3 - Predvajaj dekompresiran zvok\n"
            << "0 - IZHOD\n\n"
            << "Izbira: ";

        std::cin >> userChoice;

        switch (userChoice) {
        case 1: {
            if (audio.getNumChannels() == 0) {
                std::cout << "NAPAKA: Audio datoteka ni nalozena!\n";
                break;
            }

            std::cout << "\nIzberite velikost bloka N (64, 128, 256): ";
            std::cin >> N;
            std::cout << "Izberite faktor stiskanja M (0-" << N << "): ";
            std::cin >> M;

            if (M > N) {
                std::cout << "NAPAKA: M ne sme biti vecji od N!\n";
                break;
            }

            std::cout << "\n";
            compressAudio(audio, M, N);
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

    return 0;
}