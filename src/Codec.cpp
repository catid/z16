// Copyright (c) 2020, Christopher A. Taylor.  All rights reserved.

#include "Codec.hpp"

#include <functional>
#include <iostream>
using namespace std;

namespace z16 {


//------------------------------------------------------------------------------
// Tools

static const unsigned kAlignmentBytes = 32;

uint8_t* SIMDSafeAllocate(size_t size)
{
    uint8_t* data = (uint8_t*)calloc(1, kAlignmentBytes + size);
    if (!data) {
        return nullptr;
    }
    unsigned offset = (unsigned)((uintptr_t)data % kAlignmentBytes);
    data += kAlignmentBytes - offset;
    data[-1] = (uint8_t)offset;
    return data;
}

void SIMDSafeFree(void* ptr)
{
    if (!ptr) {
        return;
    }
    uint8_t* data = (uint8_t*)ptr;
    unsigned offset = data[-1];
    if (offset >= kAlignmentBytes) {
        return;
    }
    data -= kAlignmentBytes - offset;
    free(data);
}

static uint32_t ReadU32(const uint8_t* data)
{
    uint32_t x = data[0];
    x |= (uint32_t)data[1] << 8;
    x |= (uint32_t)data[2] << 16;
    x |= (uint32_t)data[3] << 24;
    return x;
}

static uint16_t ReadU16(const uint8_t* data)
{
    uint16_t x = data[0];
    x |= (uint16_t)data[1] << 8;
    return x;
}

static void WriteU32(uint32_t x, uint8_t* output_data)
{
    output_data[0] = (uint8_t)(x);
    output_data[1] = (uint8_t)(x >> 8);
    output_data[2] = (uint8_t)(x >> 16);
    output_data[3] = (uint8_t)(x >> 24);
}

static void WriteU16(uint16_t x, uint8_t* output_data)
{
    output_data[0] = (uint8_t)(x);
    output_data[1] = (uint8_t)(x >> 8);
}

/// Calls the provided (lambda) function at the end of the current scope
class ScopedFunction
{
public:
    ScopedFunction(std::function<void()> func) {
        Func = func;
    }
    ~ScopedFunction() {
        if (Func) {
            Func();
        }
    }

    std::function<void()> Func;

    // Prevent function from getting called
    void Cancel() {
        Func = std::function<void()>();
    }
};


//------------------------------------------------------------------------------
// Codec

bool Codec::Initialize()
{
    if (!ZstdResiduals.Initialize()) {
        return false;
    }

    return true;
}

void Codec::Shutdown()
{
    ZstdResiduals.Shutdown();
    FreeAll();
}

z16_Frame* Codec::Decode(const uint8_t* data, int32_t bytes, int32_t& error)
{
    if (bytes < kFileHeaderBytes || data[0] != kHeaderByte0 || data[1] != kHeaderByte1) {
        cerr << "Invalid header" << endl;
        error = z16_Error_InvalidInput;
        return nullptr;
    }
    if (data[2] != kHeaderVersion) {
        cerr << "Wrong version" << endl;
        error = z16_Error_InvalidInput;
        return nullptr;
    }

    error = z16_Success;

    uint32_t format = data[3];
    if (format == z16_Format_Invalid || format >= z16_Format_Count) {
        cerr << "Invalid format" << endl;
        error = z16_Error_InvalidInput;
        return nullptr;
    }

    uint32_t width = ReadU16(data + 4);
    uint32_t height = ReadU16(data + 6);
    uint32_t compressed_choices_bytes = ReadU32(data + 8);
    uint32_t compressed_residuals_bytes = ReadU32(data + 12);

    z16_Frame* frame = Allocate(format, width, height);
    if (!frame) {
        cerr << "OOM" << endl;
        error = z16_Error_OOM;
        return nullptr;
    }
    ScopedFunction frame_scope([&]() {
        Free(frame);
    });

    data += kFileHeaderBytes;
    bytes -= kFileHeaderBytes;

    // If compression was disabled:
    if (compressed_choices_bytes == 0) {
        const int image_bytes = width * 2 * height;
        if (bytes < image_bytes) {
            cerr << "Truncated" << endl;
            error = z16_Error_InvalidInput;
            return nullptr;
        }
        memcpy(frame->Image, data, image_bytes);
        frame_scope.Cancel();
        return frame;
    }

    if ((int)compressed_choices_bytes > bytes) {
        cerr << "Truncated" << endl;
        error = z16_Error_InvalidInput;
        return nullptr;
    }

    const uint8_t* compressed_choices_data = data;
    int choices_bytes = (int)ZSTD_getFrameContentSize(compressed_choices_data, bytes);
    if (choices_bytes <= 0) {
        cerr << "Invalid choices header: " << choices_bytes << endl;
        error = z16_Error_InvalidInput;
        return nullptr;
    }

    Choices.resize(choices_bytes);
    size_t choices_result = ZSTD_decompress(
        Choices.data(), choices_bytes,
        compressed_choices_data, compressed_choices_bytes);
    if (choices_result != (size_t)choices_bytes) {
        cerr << "Invalid choices data: " << choices_result << endl;
        error = z16_Error_Zstd;
        return nullptr;
    }

    const int width_blocks = (width + kBlockSize - 1) / kBlockSize;
    const int height_blocks = (height + kBlockSize - 1) / kBlockSize;
    const int expected_choices_bytes = (width_blocks * height_blocks + 7) / 8;
    if (choices_bytes != expected_choices_bytes) {
        cerr << "Invalid choices size: " << choices_bytes << " != " << expected_choices_bytes << endl;
        error = z16_Error_InvalidInput;
        return nullptr;
    }

    data += compressed_choices_bytes;
    bytes -= compressed_choices_bytes;

    if ((int)compressed_residuals_bytes > bytes) {
        cerr << "Truncated" << endl;
        error = z16_Error_InvalidInput;
        return nullptr;
    }

    const uint8_t* compressed_residuals_data = data;
    int residuals_bytes = (int)ZSTD_getFrameContentSize(compressed_residuals_data, bytes);
    if (residuals_bytes <= 0) {
        cerr << "Invalid residuals header: " << residuals_bytes << endl;
        error = z16_Error_InvalidInput;
        return nullptr;
    }

    Residuals.resize(residuals_bytes);
    size_t residuals_result = ZSTD_decompress(
        Residuals.data(), residuals_bytes,
        compressed_residuals_data, compressed_residuals_bytes);
    if (residuals_result != (size_t)residuals_bytes) {
        cerr << "Invalid residuals data: " << residuals_result << endl;
        error = z16_Error_Zstd;
        return nullptr;
    }

    const int expected_residuals_bytes = width * 2 * height;
    if (residuals_bytes != expected_residuals_bytes) {
        cerr << "Invalid residuals size: " << residuals_bytes << " != " << expected_residuals_bytes << endl;
        error = z16_Error_InvalidInput;
        return nullptr;
    }

    uint8_t* row_data = frame->Image;
    uint8_t* choice_data = Choices.data();
    int choice_bit = 0;

    uint8_t* residuals_data = Residuals.data();
    const int plane_bytes = width * 2 * height / kSliceCount;

    for (int y = 0; y < height_blocks; ++y)
    {
        for (int pixel_y = 0; pixel_y < kBlockSize; ++pixel_y)
        {
            for (int x = 0; x < width_blocks; ++x)
            {
                const int choice_i = x + choice_bit;
                const int choice = (choice_data[choice_i / 8] >> (choice_i % 8)) & 1;

                uint8_t* block_data = row_data + x * kBlockSize * 2;
                uint8_t* row_residuals = residuals_data + x * kBlockSize * kBlockSize * 2 / kSliceCount;

                if (x > 0 && y > 0 && x + 1 < width_blocks && y + 1 < height_blocks)
                {
                    if (choice == Predictor_GAP) {
                        DecodeGapBlockInner(
                            block_data, frame->StrideBytes,
                            row_residuals, plane_bytes);
                    } else {
                        DecodeMedBlockInner(
                            block_data, frame->StrideBytes,
                            row_residuals, plane_bytes);
                    }
                }
                else
                {
                    if (choice == Predictor_GAP) {
                        DecodeGapBlockOuter(
                            block_data, frame->WidthPixels,
                            x * kBlockSize, y * kBlockSize + pixel_y, frame->StrideBytes,
                            row_residuals, plane_bytes);
                    } else {
                        DecodeMedBlockOuter(
                            block_data,
                            x * kBlockSize, y * kBlockSize + pixel_y, frame->StrideBytes,
                            row_residuals, plane_bytes);
                    }
                }
            } // next block column

            residuals_data += kBlockSize * 2 / kSliceCount;
            row_data += frame->StrideBytes;
        } // next pixel row

        residuals_data += (width_blocks - 1) * kBlockSize * kBlockSize * 2 / kSliceCount;
        choice_bit += width_blocks;
    } // next block row

    frame_scope.Cancel();
    return frame;
}

z16_Frame* Codec::Allocate(int32_t format, int width, int height)
{
    if (!FreedFrames.empty()) {
        auto& old = FreedFrames.front();
        if (old->Format == format && old->WidthPixels == width && old->HeightPixels == height) {
            return old;
        } else {
            FreeAll();
        }
    }

    const int image_bytes = width * height * 2;
    uint8_t* image_data = SIMDSafeAllocate(width * height * 2 + sizeof(z16_Frame));
    z16_Frame* frame = reinterpret_cast<z16_Frame*>( image_data + image_bytes );
    frame->Format = format;
    frame->WidthPixels = width;
    frame->HeightPixels = height;
    frame->StrideBytes = width * 2;
    frame->Image = image_data;

    return frame;
}

void Codec::FreeAll()
{
    for (z16_Frame* frame : FreedFrames) {
        SIMDSafeFree(frame->Image);
    }
    FreedFrames.clear();
}

void Codec::Free(z16_Frame* frame)
{
    FreedFrames.push_back(frame);
}

uint8_t* Codec::Encode(const z16_Frame& frame, int32_t& bytes, int32_t& error)
{
    error = z16_Success;

    CandidateBuffer candidates[Predictor_Count];

    if (frame.WidthPixels % kBlockSize != 0 ||
        frame.HeightPixels % kBlockSize != 0) {
        error = z16_Error_InvalidInput;
        cerr << "Image (" << frame.WidthPixels << "x" << frame.HeightPixels << ") not a multiple of 16x16 pixels: Truncating extra pixels" << endl;
    }

    const int width_blocks = frame.WidthPixels / kBlockSize;
    const int height_blocks = frame.HeightPixels / kBlockSize;

    const int width = width_blocks * kBlockSize;
    const int height = height_blocks * kBlockSize;

    EncodedResiduals.Initialize(width, height);
    Picker.ImageReset();
    EncodedChoices.Initialize(width, height);

    // For each block (row-first):
    for (int y = 0; y < height_blocks; ++y)
    {
        for (int x = 0; x < width_blocks; ++x)
        {
            if (x > 0 && y > 0 && x + 1 < width_blocks && y + 1 < height_blocks) {
                EncodeInner(frame.Image, x * kBlockSize, y * kBlockSize, frame.StrideBytes, candidates);
            } else {
                EncodeOuter(frame.Image, width, x * kBlockSize, y * kBlockSize, frame.StrideBytes, candidates);
            }

            const int best_i = Picker.PickPredictor(candidates);

            EncodedChoices.Append(best_i);
            EncodedResiduals.Append(candidates[best_i]);
        }
    }

    const size_t choices_bound = ZSTD_compressBound(EncodedChoices.PredictorChoice.size());
    const size_t residuals_bound = ZSTD_compressBound(EncodedResiduals.Output.size());
    const size_t output_max_bytes = kFileHeaderBytes + choices_bound + residuals_bound;

    Output.resize(output_max_bytes);
    uint8_t* output_data = Output.data();

    uint8_t* choice_data = output_data + kFileHeaderBytes;
    const size_t choices_compressed_bytes = ZSTD_compress(
        choice_data, choices_bound,
        EncodedChoices.PredictorChoice.data(), EncodedChoices.PredictorChoice.size(), kZstdCompressLevel);
    if (ZSTD_isError(choices_compressed_bytes)) {
        error = z16_Error_Zstd;
        cerr << "Zstd failed (choices)" << endl;
        return nullptr;
    }

    uint8_t* residual_data = choice_data + choices_compressed_bytes;
    const size_t residual_compressed_bytes = ZstdResiduals.Compress(
        residual_data, residuals_bound,
        EncodedResiduals.Output.data(), EncodedResiduals.Output.size());
    if (ZSTD_isError(residual_compressed_bytes)) {
        error = z16_Error_Zstd;
        cerr << "Zstd failed (residuals)" << endl;
        return nullptr;
    }

    output_data[0] = kHeaderByte0;
    output_data[1] = kHeaderByte1;
    output_data[2] = kHeaderVersion;
    output_data[3] = (uint8_t)frame.Format;

    WriteU16((uint16_t)width, output_data + 4);
    WriteU16((uint16_t)height, output_data + 6);

    const size_t compressed_bytes = choices_compressed_bytes + residual_compressed_bytes;
    const size_t image_row_bytes = width * 2;
    const size_t input_bytes = image_row_bytes * height;

    if (compressed_bytes > input_bytes) {
        bytes = kFileHeaderBytes + (int32_t)input_bytes;

        // Original data is smaller so copy that over
        const uint8_t* src = frame.Image;
        uint8_t* dst = output_data + kFileHeaderBytes;
        for (int y = 0; y < height; ++y) {
            memcpy(dst, src, image_row_bytes);
            dst += image_row_bytes;
            src += frame.StrideBytes;
        }

        WriteU32(0, output_data + 8);
        WriteU32((uint32_t)input_bytes, output_data + 12);
    } else {
        bytes = kFileHeaderBytes + (int32_t)compressed_bytes;

        WriteU32((uint32_t)choices_compressed_bytes, output_data + 8);
        WriteU32((uint32_t)residual_compressed_bytes, output_data + 12);
    }

    return Output.data();
}


} // namespace z16
