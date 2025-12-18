#ifndef BITWRITER_H
#define BITWRITER_H

#include <fstream>
#include <cstdint>

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

#endif