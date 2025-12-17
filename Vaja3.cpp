#include <iostream>
#include "AudioFile.h"
#include <map>

using namespace std;

// - Izbira vhodne nekompresirane zvoène datoteke tipa WAV.
// - Faktor stiskanja(parameter M) in velikost bloka(parameter N, npr. 64, 128, itd.).
// - Shranjevanje kompresiranega zvoka in izvedba dekompresije nad prebranim zvokom.
// - Predvajanje dekompresiranega zvoka.Alternativno lahko hranite dekompresiran zvok na disk v formatu WAV.

struct MidSide {
    AudioFile<uint32_t> mid;
    AudioFile<uint32_t> side;
};

MidSide splitToMidAndSide(const AudioFile<uint32_t>& audio) {
    MidSide ms;

    if (audio.getNumChannels() != 2) {
        throw runtime_error("Audio ni stereo!");
    }

    int numSamplesPerChannel = audio.getNumSamplesPerChannel();
    
    ms.mid.setNumChannels(1);
    ms.side.setNumChannels(1);
    ms.mid.setNumSamplesPerChannel(numSamplesPerChannel);
    ms.side.setNumSamplesPerChannel(numSamplesPerChannel);
    ms.mid.setSampleRate(audio.getSampleRate());
    ms.side.setSampleRate(audio.getSampleRate());

    for (int i = 0; i < numSamplesPerChannel; ++i) {
        uint32_t L = audio.samples[0][i];
        uint32_t R = audio.samples[1][i];

        ms.mid.samples[0][i] = (L + R) / 2;
        ms.side.samples[0][i] = (L - R) / 2;
    }

    return ms;
}

void compressAudio(const AudioFile<uint32_t>& audio, int compressionFactor, int blockSize) {
    MidSide ms = splitToMidAndSide(audio);
    AudioFile<uint32_t> M = ms.mid;
    AudioFile<uint32_t> S = ms.side;

    //BUILD BLOCKS
    vector<vector<uint32_t>> blocks(blockSize, vector<uint32_t>(blockSize - 1, 0));

    int half = blockSize / 2;
    int lastCol = blockSize - 1;

    for (int row = 0; row < half; ++row) {
        blocks[row][0] = 0;
    }

    for (int row = blockSize - half; row < blockSize; ++row) {
        blocks[row][lastCol] = 0;
    }
}

void dekompresijaZvoka() {
    cout << "Izbrali ste dekompresijo zvoka!";
}

void predvajajDekompresiranZvok() {
    cout << "Izbrali ste predvajanje dekompresiranjega zvoka!";
}

int main()
{
    AudioFile<uint32_t> audio;
    int userChoice = 0, compressionFactor, blockSize;

    while (true) {
        cout << "Izberite nacin:\n\n"
            << "1 - Kompresija zvoka\n"
            << "2 - Izvedi dekompresijo\n"
            << "3 - Predvajaj dekompresiran zvok\n"
            << "0 - IZHOD\n\n";

        cin >> userChoice;

        switch (userChoice)
        {
        case 1:
            cout << "Izberite faktor stiskanja M:\n";
            cin >> compressionFactor;
            cout << "Izberite velikost bloka N:\n";
            cin >> blockSize;
            cout << "\n\n";

            compressAudio(audio, compressionFactor, blockSize * 2);
            break;

        case 2:
            dekompresijaZvoka();
            break;

        case 3:
            predvajajDekompresiranZvok();
            break;

        case 0:
            cout << "\nIzhod iz programa.\n";
            return 0;

        default:
            cout << "Neveljavna izbira!\n";
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
