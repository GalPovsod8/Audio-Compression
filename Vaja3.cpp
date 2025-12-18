#include <iostream>
#include "AudioFile.h"
#include <vector>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <string>
#include <chrono>

#ifdef _WIN32
#include <windows.h>
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

std::string openFileDialog(const std::string& title, const std::string& filter) {
#ifdef _WIN32
    char filename[MAX_PATH] = { 0 };
    OPENFILENAMEA ofn;
    ZeroMemory(&ofn, sizeof(ofn));
    ofn.lStructSize = sizeof(ofn);
    ofn.hwndOwner = NULL;
    ofn.lpstrFile = filename;
    ofn.nMaxFile = sizeof(filename);
    ofn.lpstrFilter = filter.c_str();
    ofn.nFilterIndex = 1;
    ofn.lpstrFileTitle = NULL;
    ofn.nMaxFileTitle = 0;
    ofn.lpstrInitialDir = NULL;
    ofn.lpstrTitle = title.c_str();
    ofn.Flags = OFN_PATHMUSTEXIST | OFN_FILEMUSTEXIST;

    if (GetOpenFileNameA(&ofn) == TRUE) {
        return std::string(filename);
    }
#else
    std::cout << title << "\nVnesite pot do datoteke: ";
    std::string path;
    std::getline(std::cin, path);
    return path;
#endif
    return "";
}

std::string saveFileDialog(const std::string& title, const std::string& filter, const std::string& defaultExt) {
#ifdef _WIN32
    char filename[MAX_PATH] = { 0 };
    OPENFILENAMEA ofn;
    ZeroMemory(&ofn, sizeof(ofn));
    ofn.lStructSize = sizeof(ofn);
    ofn.hwndOwner = NULL;
    ofn.lpstrFile = filename;
    ofn.nMaxFile = sizeof(filename);
    ofn.lpstrFilter = filter.c_str();
    ofn.nFilterIndex = 1;
    ofn.lpstrFileTitle = NULL;
    ofn.nMaxFileTitle = 0;
    ofn.lpstrInitialDir = NULL;
    ofn.lpstrTitle = title.c_str();
    ofn.lpstrDefExt = defaultExt.c_str();
    ofn.Flags = OFN_PATHMUSTEXIST | OFN_OVERWRITEPROMPT;

    if (GetSaveFileNameA(&ofn) == TRUE) {
        return std::string(filename);
    }
#else
    std::cout << title << "\nVnesite pot za shranjevanje: ";
    std::string path;
    std::getline(std::cin, path);
    return path;
#endif
    return "";
}

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
    int absValue = std::abs(value);
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

//////////////////////////////// KOMPRESIJA //////////////////////

void writeMDCTBlock(BitWriter& writer, const std::vector<int>& mdctCoeffs, int M) {
    int N = mdctCoeffs.size();
    for (int i = 0; i < (N - M); ++i) {
        int coeff = mdctCoeffs[i];
        int bitLength = getBitLength(coeff);
        writer.writeBits(bitLength, 6);
        if (bitLength > 0) {
            uint32_t encodedCoeff = encodeValue(coeff, bitLength);
            writer.writeBits(encodedCoeff, bitLength);
        }
    }
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

std::vector<double> MDCT(const std::vector<double>& x) {
    int twoN = x.size();
    int N = twoN / 2;
    std::vector<double> X(N, 0.0);

    for (int k = 0; k < N; ++k) {
        double sum = 0.0;
        for (int n = 0; n < twoN; ++n) {
            sum += x[n] * std::cos(M_PI / N * (n + 0.5 + N / 2.0) * (k + 0.5));
        }
        X[k] = sum;
    }
    return X;
}

std::vector<std::vector<double>> MDCTAllBlocks(const std::vector<std::vector<double>>& blocks) {
    std::vector<std::vector<double>> mdctBlocks;
    mdctBlocks.reserve(blocks.size());
    for (const auto& block : blocks) {
        mdctBlocks.push_back(MDCT(block));
    }
    return mdctBlocks;
}

void quantizeAndCompress(std::vector<std::vector<double>>& mdctBlocks, int M,
    std::vector<std::vector<int>>& intBlocks) {
    intBlocks.clear();
    intBlocks.reserve(mdctBlocks.size());

    for (auto& block : mdctBlocks) {
        int N = block.size();
        std::vector<int> intBlock(N);

        for (int i = 0; i < N; ++i) {
            intBlock[i] = static_cast<int>(std::round(block[i]));
        }

        for (int i = N - M; i < N; ++i) {
            intBlock[i] = 0;
        }

        intBlocks.push_back(intBlock);
    }
}

bool writeCompressedAudio(const std::string& filename,
    const std::vector<std::vector<int>>& midBlocks,
    const std::vector<std::vector<int>>& sideBlocks,
    int numSamples, int N, int sampleRate, int M) {

    std::ofstream outFile(filename, std::ios::binary);
    if (!outFile) {
        std::cerr << "Napaka: Ne morem odpreti datoteke za pisanje!\n";
        return false;
    }

    std::cout << "\n=== ZAPIS V DATOTEKO ===\n";
    std::cout << "Datoteka: " << filename << "\n";

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
    std::cout << "  M: " << M << "\n\n";

    BitWriter writer(outFile);
    int totalBlocks = midBlocks.size() + sideBlocks.size();
    int written = 0;

    std::cout << "Pisanje Mid blokov...\n";
    for (const auto& block : midBlocks) {
        writeMDCTBlock(writer, block);
        written++;
        if (written % 100 == 0) {
            std::cout << "\r  Zapisano: " << written << "/" << midBlocks.size() << std::flush;
        }
    }
    std::cout << "\n";

    std::cout << "Pisanje Side blokov...\n";
    for (const auto& block : sideBlocks) {
        writeMDCTBlock(writer, block);
        written++;
        if (written % 100 == 0) {
            std::cout << "\r  Zapisano: " << (written - midBlocks.size())
                << "/" << sideBlocks.size() << std::flush;
        }
    }
    std::cout << "\n";

    writer.flush();
    outFile.close();

    std::cout << "Zapis zakljucen!\n";
    std::cout << "Skupaj blokov: " << totalBlocks << "\n\n";
    return true;
}

bool compressAudio(const std::string& inputPath, const std::string& outputPath,
    int M, int N, double& compressionTime) {

    auto startTime = std::chrono::high_resolution_clock::now();

    std::cout << "\n========================================\n";
    std::cout << "    KOMPRESIJA ZVOKA\n";
    std::cout << "========================================\n";
    std::cout << "Vhodna datoteka: " << inputPath << "\n";
    std::cout << "Izhodna datoteka: " << outputPath << "\n";
    std::cout << "N = " << N << ", M = " << M << "\n";
    std::cout << "Velikost bloka: " << 2 * N << "\n\n";

    AudioFile<double> audio;
    if (!audio.load(inputPath)) {
        std::cerr << "Napaka: Ne morem naloziti audio datoteke!\n";
        return false;
    }

    std::cout << "Audio datoteka nalozena:\n";
    std::cout << "  Kanali: " << audio.getNumChannels() << "\n";
    std::cout << "  Vzorci: " << audio.getNumSamplesPerChannel() << "\n";
    std::cout << "  Sample rate: " << audio.getSampleRate() << " Hz\n";
    std::cout << "  Trajanje: " << std::fixed << std::setprecision(2)
        << (double)audio.getNumSamplesPerChannel() / audio.getSampleRate()
        << " sekund\n\n";

    std::cout << "1. Mid-Side razdelitev...\n";
    MidSide ms = splitToMidAndSide(audio);

    std::cout << "2. Gradnja blokov...\n";
    std::vector<std::vector<double>> midBlocks, sideBlocks;
    buildBlocks(ms.mid.samples[0], N, midBlocks);
    buildBlocks(ms.side.samples[0], N, sideBlocks);
    std::cout << "   Mid: " << midBlocks.size() << " blokov\n";
    std::cout << "   Side: " << sideBlocks.size() << " blokov\n";

    std::cout << "3. Aplikacija okenske funkcije...\n";
    applyWindow(midBlocks);
    applyWindow(sideBlocks);

    std::cout << "4. MDCT transformacija...\n";
    std::vector<std::vector<double>> midMDCT = MDCTAllBlocks(midBlocks);
    std::vector<std::vector<double>> sideMDCT = MDCTAllBlocks(sideBlocks);

    std::cout << "5. Kvantizacija in kompresija...\n";
    std::vector<std::vector<int>> midInt, sideInt;
    quantizeAndCompress(midMDCT, M, midInt);
    quantizeAndCompress(sideMDCT, M, sideInt);

    std::cout << "6. Zapis v datoteko...\n";
    if (!writeCompressedAudio(outputPath, midInt, sideInt,
        audio.getNumSamplesPerChannel(),
        N, audio.getSampleRate(), M)) {
        return false;
    }

    auto endTime = std::chrono::high_resolution_clock::now();
    compressionTime = std::chrono::duration<double>(endTime - startTime).count();

    std::cout << "========================================\n";
    std::cout << "KOMPRESIJA USPESNA!\n";
    std::cout << "Cas: " << std::fixed << std::setprecision(3)
        << compressionTime << " sekund\n";
    std::cout << "========================================\n\n";

    return true;
}

///////////////////////   DEKOMPRESIJA   /////////////////////////

void combineMidSide(const std::vector<double>& mid, const std::vector<double>& side, AudioFile<double>& audio) {
    int numSamples = mid.size();
    audio.setNumChannels(2);
    audio.setNumSamplesPerChannel(numSamples);

    for (int i = 0; i < numSamples; ++i) {
        audio.samples[0][i] = mid[i] + side[i];
        audio.samples[1][i] = mid[i] - side[i];
    }
}

std::vector<double> IMDCT(const std::vector<int>& X) {
    int N = X.size();
    int twoN = 2 * N;
    std::vector<double> y(twoN, 0.0);

    for (int n = 0; n < twoN; ++n) {
        double sum = 0.0;
        for (int k = 0; k < N; ++k) {
            sum += X[k] * std::cos(
                M_PI / N * (n + 0.5 + N / 2.0) * (k + 0.5)
            );
        }
        y[n] = sum;
    }
    return y;
}

std::vector<std::vector<double>> IMDCTAllBlocks(const std::vector<std::vector<int>>& mdctBlocks)
{
    std::vector<std::vector<double>> timeBlocks;
    timeBlocks.reserve(mdctBlocks.size());

    for (const auto& block : mdctBlocks) {
        timeBlocks.push_back(IMDCT(block));
    }

    return timeBlocks;
}

void applyWindow1D(std::vector<double>& block) {
    int N = block.size();
    for (int n = 0; n < N; ++n) {
        block[n] *= w(n, N);
    }
}

std::vector<double> overlapAdd(const std::vector<std::vector<double>>& blocks, int N) {
    int numBlocks = blocks.size();
    int outSize = (numBlocks + 1) * N;
    std::vector<double> output(outSize, 0.0);

    for (int b = 0; b < numBlocks; ++b) {
        for (int n = 0; n < 2 * N; ++n) {
            output[b * N + n] += blocks[b][n];
        }
    }

    output.erase(output.begin(), output.begin() + N);
    output.resize(output.size() - N);

    return output;
}

bool decompressAudio(const std::string& inputPath, const std::string& outputPath) {
    std::ifstream inFile(inputPath, std::ios::binary);
    if (!inFile) return false;

    uint32_t numSamples;
    uint16_t N;
    uint32_t sampleRate;
    uint16_t M;

    inFile.read(reinterpret_cast<char*>(&numSamples), sizeof(uint32_t));
    inFile.read(reinterpret_cast<char*>(&N), sizeof(uint16_t));
    inFile.read(reinterpret_cast<char*>(&sampleRate), sizeof(uint32_t));
    inFile.read(reinterpret_cast<char*>(&M), sizeof(uint16_t));

    BitReader reader(inFile);

    int paddedSize = numSamples + 2 * N;
    int numBlocks = (paddedSize - N) / N;

    auto readChannelBlocks = [&](int count) {
        std::vector<std::vector<int>> blocks(count, std::vector<int>(N, 0));
        for (int i = 0; i < count; ++i) {
            for (int j = 0; j < (N - M); ++j) {
                int bitLength = reader.readBits(6);
                if (bitLength > 0) {
                    uint32_t val = reader.readBits(bitLength);
                    blocks[i][j] = decodeValue(val, bitLength);
                }
            }
        }
        return blocks;
        };

    std::cout << "Berem bloke...\n";
    std::vector<std::vector<int>> midIntBlocks = readChannelBlocks(numBlocks);
    std::vector<std::vector<int>> sideIntBlocks = readChannelBlocks(numBlocks);

    auto processBlocks = [&](std::vector<std::vector<int>>& intBlocks) {
        std::vector<std::vector<double>> timeBlocks;
        for (const auto& b : intBlocks) {
            std::vector<double> ib = IMDCT(b);
            applyWindow1D(ib);
            timeBlocks.push_back(ib);
        }
        return overlapAdd(timeBlocks, N);
        };

    std::cout << "Izvajam IMDCT in Overlap-Add...\n";
    std::vector<double> midSignal = processBlocks(midIntBlocks);
    std::vector<double> sideSignal = processBlocks(sideIntBlocks);

    AudioFile<double> outputAudio;
    outputAudio.setSampleRate(sampleRate);
    combineMidSide(midSignal, sideSignal, outputAudio);

    if (outputAudio.getNumSamplesPerChannel() > numSamples) {
        outputAudio.setNumSamplesPerChannel(numSamples);
    }

    return outputAudio.save(outputPath);
}

int main() {
    std::cout << "========================================\n";
    std::cout << "    MDCT KOMPRESIJA ZVOKA\n";
    std::cout << "========================================\n\n";

    int userChoice = 0, M, N;

    while (true) {
        std::cout << "\nIzbira:\n\n"
            << "1 - Kompresija zvoka\n"
            << "2 - Dekompresija zvoka\n"
            << "3 - Testiranje (vsi N in M)\n"
            << "0 - IZHOD\n\n"
            << "Vnesite izbiro: ";

        std::cin >> userChoice;

        switch (userChoice) {
        case 1: {
            std::cout << "\n=== KOMPRESIJA ===\n";

            std::string inputPath = openFileDialog(
                "Izberite vhodno WAV datoteko",
                "WAV Files\0*.wav\0All Files\0*.*\0"
            );

            if (inputPath.empty()) {
                std::cout << "Izbira preklicana.\n";
                break;
            }

            std::cout << "Izbrana vhodna datoteka: " << inputPath << "\n";

            std::string outputPath = saveFileDialog(
                "Shrani kompresirano datoteko",
                "MDCT Files\0*.mdct\0All Files\0*.*\0",
                "mdct"
            );

            if (outputPath.empty()) {
                std::cout << "Izbira preklicana.\n";
                break;
            }

            std::cout << "Izbrana izhodna datoteka: " << outputPath << "\n\n";

            std::cout << "Vnesite N (64, 128, 256): ";
            std::cin >> N;
            std::cout << "Vnesite M (0 do " << N << "): ";
            std::cin >> M;

            if (M > N || M < 0) {
                std::cout << "Napaka: M mora biti med 0 in " << N << "!\n";
                break;
            }

            double compressionTime;
            if (compressAudio(inputPath, outputPath, M, N, compressionTime)) {
                std::ifstream inFile(inputPath, std::ios::binary | std::ios::ate);
                std::ifstream outFile(outputPath, std::ios::binary | std::ios::ate);

                if (inFile && outFile) {
                    size_t inSize = inFile.tellg();
                    size_t outSize = outFile.tellg();
                    double ratio = (double)inSize / outSize;

                    std::cout << "\n--- STATISTIKA ---\n";
                    std::cout << "Originalna velikost: " << inSize << " B\n";
                    std::cout << "Kompresirana velikost: " << outSize << " B\n";
                    std::cout << "Kompresijsko razmerje: " << std::fixed
                        << std::setprecision(2) << ratio << ":1\n";
                    std::cout << "Prihranek: " << std::fixed << std::setprecision(1)
                        << (1.0 - (double)outSize / inSize) * 100 << "%\n";
                }
            }
            break;
        }

        case 2: {
            std::string inputPath = openFileDialog("Izberite .mdct datoteko", "MDCT Files\0*.mdct\0");
            if (inputPath.empty()) break;

            std::string outputPath = saveFileDialog("Shrani kot WAV", "WAV Files\0*.wav\0", "wav");
            if (outputPath.empty()) break;

            if (decompressAudio(inputPath, outputPath)) {
                std::cout << "Dekompresija uspešna: " << outputPath << "\n";
            }
            else {
                std::cout << "Napaka pri dekompresiji!\n";
            }
            break;
        }

        case 3: {
            std::vector<int> N_vals = { 64, 128, 256 };
            std::vector<int> M_vals = { 0, 32, 60 };
            std::string testFile = openFileDialog("Izberi WAV za test", "WAV Files\0*.wav\0");

            std::cout << "N\tM\tRatio\tTime(s)\n";
            for (int n : N_vals) {
                for (int m : M_vals) {
                    if (m >= n) continue;
                    double cTime;
                    std::string outName = "test_N" + std::to_string(n) + "_M" + std::to_string(m) + ".mdct";
                    compressAudio(testFile, outName, m, n, cTime);
                }
            }
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