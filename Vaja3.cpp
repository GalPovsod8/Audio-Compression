#include <iostream>
#include "AudioFile.h"
#include <vector>
#include <cmath>
#include <iomanip>
#include <fstream>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

class BitWriter {
private:
    std::ofstream& out;
    uint8_t buffer;
    int bitsFilled;

public:
    BitWriter(std::ofstream& file) : out(file), buffer(0), bitsFilled(0) {}

    void writeBit(bool bit) {
        buffer = (buffer << 1) | (bit ? 1 : 0);
        bitsFilled++;
        if (bitsFilled == 8) {
            out.put(static_cast<char>(buffer));
            buffer = 0;
            bitsFilled = 0;
        }
    }

    void writeBits(uint32_t value, int nBits) {
        for (int i = nBits - 1; i >= 0; --i) {
            writeBit((value >> i) & 1);
        }
    }

    void flush() {
        if (bitsFilled > 0) {
            buffer <<= (8 - bitsFilled);
            out.put(static_cast<char>(buffer));
            buffer = 0;
            bitsFilled = 0;
        }
    }
};


class BitReader {
private:
    std::ifstream& in;
    uint8_t buffer;
    int bitsRemaining;

public:
    BitReader(std::ifstream& file) : in(file), buffer(0), bitsRemaining(0) {}

    void fillBuffer() {
        if (in.eof()) {
            bitsRemaining = 0;
            return;
        }
        char byte;
        if (in.get(byte)) {
            buffer = static_cast<uint8_t>(byte);
            bitsRemaining = 8;
        }
        else {
            bitsRemaining = 0;
        }
    }

    bool readBit() {
        if (bitsRemaining == 0) fillBuffer();
        if (bitsRemaining == 0) return false;

        bool bit = (buffer & 0x80) != 0;
        buffer <<= 1;
        bitsRemaining--;
        return bit;
    }

    uint32_t readBits(int nBits) {
        uint32_t value = 0;
        for (int i = 0; i < nBits; i++) {
            value = (value << 1) | (readBit() ? 1 : 0);
        }
        return value;
    }

    bool hasMoreData() {
        if (bitsRemaining > 0) return true;
        fillBuffer();
        return bitsRemaining > 0;
    }
};

int getBitLength(int value) {
    if (value == 0) return 0;
    int absValue = abs(value);
    int bits = 0;
    while (absValue > 0) {
        absValue >>= 1;
        bits++;
    }
    return bits;
}

uint32_t encodeValue(int value, int bitLength) {
    if (value >= 0) {
        return static_cast<uint32_t>(value);
    }
    else {
        return (1u << bitLength) + value;
    }
}

int decodeValue(uint32_t encodedValue, int bitLength) {
    if (bitLength == 0) return 0;

    uint32_t signBit = 1u << (bitLength - 1);
    if (encodedValue >= signBit) {
        return static_cast<int>(encodedValue) - (1 << bitLength);
    }
    else {
        return static_cast<int>(encodedValue);
    }
}

void writeMDCTBlock(BitWriter& writer, const std::vector<int>& mdctCoeffs) {
    for (int coeff : mdctCoeffs) {
        int bitLength = getBitLength(coeff);

        writer.writeBits(bitLength, 6);

        if (bitLength > 0) {
            uint32_t encodedCoeff = encodeValue(coeff, bitLength);
            writer.writeBits(encodedCoeff, bitLength);
        }
    }
}

void writeCompressedAudio(const std::string& filename,
    const std::vector<std::vector<int>>& midBlocks,
    const std::vector<std::vector<int>>& sideBlocks,
    int numSamples,
    int N,
    int sampleRate,
    int M) {

    std::ofstream outFile(filename, std::ios::binary);
    if (!outFile) {
        throw std::runtime_error("Ne morem odpreti datoteke za pisanje!");
    }

    std::cout << "\n=== ZAPIS V DATOTEKO ===\n";

    uint32_t numSamples32 = static_cast<uint32_t>(numSamples);
    uint16_t N16 = static_cast<uint16_t>(N);
    uint32_t sampleRate32 = static_cast<uint32_t>(sampleRate);
    uint16_t M16 = static_cast<uint16_t>(M);

    outFile.write(reinterpret_cast<const char*>(&numSamples32), sizeof(uint32_t));
    outFile.write(reinterpret_cast<const char*>(&N16), sizeof(uint16_t));
    outFile.write(reinterpret_cast<const char*>(&sampleRate32), sizeof(uint32_t));
    outFile.write(reinterpret_cast<const char*>(&M16), sizeof(uint16_t));

    std::cout << "Glava datoteke:\n";
    std::cout << "  Stevilo vzorcev: " << numSamples << "\n";
    std::cout << "  N: " << N << "\n";
    std::cout << "  Sample rate: " << sampleRate << " Hz\n";
    std::cout << "  M (faktor stiskanja): " << M << "\n\n";

    BitWriter writer(outFile);

    int totalBlocks = midBlocks.size() + sideBlocks.size();
    int written = 0;

    std::cout << "Pisanje Mid blokov...\n";
    for (const auto& block : midBlocks) {
        writeMDCTBlock(writer, block);
        written++;
        if (written % 100 == 0) {
            std::cout << "\rZapisano: " << written << "/" << totalBlocks << " blokov" << std::flush;
        }
    }

    std::cout << "\nPisanje Side blokov...\n";
    for (const auto& block : sideBlocks) {
        writeMDCTBlock(writer, block);
        written++;
        if (written % 100 == 0) {
            std::cout << "\rZapisano: " << written << "/" << totalBlocks << " blokov" << std::flush;
        }
    }

    writer.flush();
    outFile.close();

    std::cout << "\n\nZapis zakljucen!\n";
    std::cout << "Velikost datoteke: " << written << " blokov\n\n";
}

std::vector<int> readMDCTBlock(BitReader& reader, int N) {
    std::vector<int> mdctCoeffs(N, 0);

    for (int i = 0; i < N; ++i) {
        if (!reader.hasMoreData()) break;

        int bitLength = reader.readBits(6);

        if (bitLength > 0) {
            uint32_t encodedCoeff = reader.readBits(bitLength);
            mdctCoeffs[i] = decodeValue(encodedCoeff, bitLength);
        }
        else {
            mdctCoeffs[i] = 0;
        }
    }

    return mdctCoeffs;
}

bool readCompressedAudio(const std::string& filename,
    std::vector<std::vector<int>>& midBlocks,
    std::vector<std::vector<int>>& sideBlocks,
    int& numSamples,
    int& N,
    int& sampleRate,
    int& M) {

    std::ifstream inFile(filename, std::ios::binary);
    if (!inFile) {
        std::cerr << "Ne morem odpreti datoteke za branje!\n";
        return false;
    }

    std::cout << "\n=== BRANJE IZ DATOTEKE ===\n";

    uint32_t numSamples32;
    uint16_t N16;
    uint32_t sampleRate32;
    uint16_t M16;

    inFile.read(reinterpret_cast<char*>(&numSamples32), sizeof(uint32_t));
    inFile.read(reinterpret_cast<char*>(&N16), sizeof(uint16_t));
    inFile.read(reinterpret_cast<char*>(&sampleRate32), sizeof(uint32_t));
    inFile.read(reinterpret_cast<char*>(&M16), sizeof(uint16_t));

    numSamples = numSamples32;
    N = N16;
    sampleRate = sampleRate32;
    M = M16;

    std::cout << "Glava datoteke:\n";
    std::cout << "  Stevilo vzorcev: " << numSamples << "\n";
    std::cout << "  N: " << N << "\n";
    std::cout << "  Sample rate: " << sampleRate << " Hz\n";
    std::cout << "  M: " << M << "\n\n";

    int paddedSize = numSamples + 2 * N;
    int numBlocks = (paddedSize - N) / N;

    std::cout << "Pricakovano stevilo blokov na kanal: " << numBlocks << "\n\n";

    BitReader reader(inFile);

    std::cout << "Branje Mid blokov...\n";
    midBlocks.clear();
    for (int i = 0; i < numBlocks; ++i) {
        midBlocks.push_back(readMDCTBlock(reader, N));
        if ((i + 1) % 100 == 0) {
            std::cout << "\rPrebrano: " << (i + 1) << "/" << numBlocks << " blokov" << std::flush;
        }
    }

    std::cout << "\nBranje Side blokov...\n";
    sideBlocks.clear();
    for (int i = 0; i < numBlocks; ++i) {
        sideBlocks.push_back(readMDCTBlock(reader, N));
        if ((i + 1) % 100 == 0) {
            std::cout << "\rPrebrano: " << (i + 1) << "/" << numBlocks << " blokov" << std::flush;
        }
    }

    inFile.close();

    std::cout << "\n\nBranje zakljuceno!\n";
    std::cout << "Mid blokov: " << midBlocks.size() << "\n";
    std::cout << "Side blokov: " << sideBlocks.size() << "\n\n";

    return true;
}


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