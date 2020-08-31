// Copyright (c) 2020, Christopher A. Taylor.  All rights reserved.

#include "Tools.hpp"

namespace z16 {


//------------------------------------------------------------------------------
// Logarithms

/// Base-2 logarithm, rounding down
/// Precondition: x != 0
inline unsigned Log2Floor32(uint32_t x)
{
#ifdef _MSC_VER
    unsigned long index;
    // Note: Ignoring result because x != 0
    _BitScanReverse(&index, x);
    return (unsigned)index;
#else
    // Note: Ignoring return value of 0 because x != 0
    return 31 - (unsigned)__builtin_clz(x);
#endif
}

// Base-2 logarithm, rounding up
inline unsigned Log2Ceil32(uint32_t x)
{
    if (x <= 1) {
        return 0;
    }
    return Log2Floor32(x - 1) + 1;
}


//------------------------------------------------------------------------------
// Tools

static inline uint16_t ZigZag(int16_t x)
{
    return (x << 1) ^ (x >> 15);
}

static inline int16_t UnZigZag(uint16_t x)
{
    return (x >> 1) ^ -(x & 1);
}

static inline int Clamp(int x, int maxval)
{
    if (x < 0) {
        return 0;
    }
    if (x > maxval) {
        return maxval;
    }
    return x;
}

static inline int Gradient(int n, int w, int nw, int maxval)
{
    return Clamp(n + w - nw, maxval);
}

static inline int Median(int a, int b, int c)
{ 
    const int x = a - b;
    const int y = b - c;
    const int z = a - c;

    if (x * y > 0) {
        return b;
    }
    if (x * z > 0)  {
        return c;
    }
    return a;
}


//------------------------------------------------------------------------------
// GAP Predictor [1]

static inline uint16_t GapPredict16(
    int n2, int ne2,
    int nw, int n, int ne,
    int w2, int w)
{
    const int thresh = 80;

    const int dh = std::abs(w2 - w) + std::abs(nw - n) + std::abs(n - ne);
    const int dv = std::abs(w - nw) + std::abs(n - n2) + std::abs(ne2 - ne);

    const int dvh = dv - dh;
    if (dvh > thresh) {
        return static_cast<uint16_t>( w );
    }
    if (dvh < -thresh) {
        return static_cast<uint16_t>( n );
    }

    const int gap = (w + n) / 2 + (ne - nw) / 4;
    return static_cast<uint16_t>( Clamp(gap, 65535) );
}


//------------------------------------------------------------------------------
// MED Predictor [2]

static inline uint16_t MedPredict16(uint16_t n, uint16_t w, uint16_t nw)
{
    return static_cast<uint16_t>( Median(n, w, Gradient(n, w, nw, 65535)) );
}


//------------------------------------------------------------------------------
// PredictorPicker

void PredictorPicker::ImageReset()
{
    SliceTotal = 0;

    for (int j = 0; j < kSliceCount; ++j) {
        for (int k = 0; k < kByteLevels; ++k) {
            ImageCounters[j][k] = 0;
        }
    }
}

int PredictorPicker::PickPredictor(const CandidateBuffer buffers[Predictor_Count])
{
    return 0;

    int best_i = 0;

    // If there is exemplar data:
    if (SliceTotal > 0)
    {
        for (int i = 0; i < Predictor_Count; ++i) {
            // Zero histogram
            for (int j = 0; j < kSliceCount; ++j) {
                for (int k = 0; k < kByteLevels; ++k) {
                    BlockCounters[i][j][k] = 0;
                }
            }

            const uint8_t* buffer = buffers[i].Buffer;
            const int slice_bytes = kBlockSize * kBlockSize * 2 / kSliceCount;

            // Fill histogram
            for (int j = 0; j < kSliceCount; ++j) {
                uint32_t* hist = BlockCounters[i][j];

                for (int k = 0; k < slice_bytes; ++k) {
                    const int value = buffer[k];
                    hist[value]++;
                }

                buffer += slice_bytes;
            }
        }

        const unsigned scale = 1;
        const unsigned scaled_total = SliceTotal * scale;

        unsigned predictor_bits[Predictor_Count] = {0};

        for (int j = 0; j < kSliceCount; ++j) {
            for (int k = 0; k < kByteLevels; ++k) {
                const uint32_t image_count = ImageCounters[j][k];
                if (image_count <= 0) {
                    continue;
                }

                const unsigned bits_per_pixel = Log2Ceil32(scaled_total / image_count);

                for (int i = 0; i < Predictor_Count; ++i) {
                    predictor_bits[i] += bits_per_pixel * BlockCounters[i][j][k];
                }
            }
        }

        // Pick the one with the lowest total expected bit representation
        unsigned best_bits = predictor_bits[0];
        for (int i = 1; i < Predictor_Count; ++i) {
            if (best_bits > predictor_bits[i]) {
                best_bits = predictor_bits[i];
                best_i = i;
            }
        }
    } // end if

    // Incorporate into the image counters
    for (int j = 0; j < kSliceCount; ++j) {
        for (int k = 0; k < kByteLevels; ++k) {
            ImageCounters[j][k] += BlockCounters[best_i][j][k];
        }
    }

    SliceTotal += kBlockSize * kBlockSize / kSliceCount;

    return best_i;
}


//------------------------------------------------------------------------------
// EncodeBuffer

void EncodeBuffer::Initialize(int width, int height)
{
    Output.resize(width * height * 2);
    PlaneBytes = width * height * 2 / kSliceCount;
    WriteOffset = 0;
}

void EncodeBuffer::Append(const CandidateBuffer& buffer)
{
    const uint8_t* src = buffer.Buffer;
    uint8_t* dest = Output.data() + WriteOffset;

    const int slice_bytes = kBlockSize * kBlockSize * 2 / kSliceCount;

    for (int i = 0; i < kSliceCount; ++i) {
        memcpy(dest, src, slice_bytes);
        dest += PlaneBytes;
        src += slice_bytes;
    }

    WriteOffset += slice_bytes;
}


//------------------------------------------------------------------------------
// Encoding

// Input size: kBlockSize * kBlockSize * 2
void Shuffle(const CandidateBuffer& source_buffer, CandidateBuffer& dest_buffer)
{
    const uint16_t* words = reinterpret_cast<const uint16_t*>( source_buffer.Buffer );
    uint8_t* dest = dest_buffer.Buffer;

    static_assert(kSliceCount == 4, "FIXME");

    const int slice_bytes = kBlockSize * kBlockSize * 2 / kSliceCount;
    for (int i = 0; i < slice_bytes; ++i)
    {
        const uint16_t w0 = words[i * 2 + 0];
        const uint16_t w1 = words[i * 2 + 1];

        // Spread every 4 bits into a different plane
        dest[0] = (w0 & 15) | ((w1 & 15) << 4);
        dest[slice_bytes] = ((w0 >> 4) & 15) | (((w1 >> 4) & 15) << 4);
        dest[slice_bytes * 2] = ((w0 >> 8) & 15) | (((w1 >> 8) & 15) << 4);
        dest[slice_bytes * 3] = ((w0 >> 12) & 15) | (((w1 >> 12) & 15) << 4);

        ++dest;
    }
}

void EncodeInner(const uint8_t* image_data, int pixel_offset_x, int pixel_offset_y, int stride_bytes, CandidateBuffer buffers[Predictor_Count])
{
    const uint8_t* row_data = image_data + pixel_offset_x * 2 + pixel_offset_y * stride_bytes;

    const uint16_t* prev_row = reinterpret_cast<const uint16_t*>( row_data - stride_bytes );
    const uint16_t* prev_row2 = reinterpret_cast<const uint16_t*>( row_data - stride_bytes * 2 );

    CandidateBuffer candidates[Predictor_Count];

    for (int y = 0; y < kBlockSize; ++y)
    {
        const uint16_t* row = reinterpret_cast<const uint16_t*>( row_data );
        uint16_t* gap_dest = reinterpret_cast<uint16_t*>( candidates[Predictor_GAP].Buffer + y * kBlockSize * 2 );
        uint16_t* med_dest = reinterpret_cast<uint16_t*>( candidates[Predictor_MED].Buffer + y * kBlockSize * 2 );

        for (int x = 0; x < kBlockSize; ++x)
        {
            /*
                        n2  ne2
                    nw  n   ne
                w2  w   ?
            */
            const uint16_t n2 = prev_row2[x];
            const uint16_t ne2 = prev_row2[x + 1];
            const uint16_t nw = prev_row[x - 1];
            const uint16_t n = prev_row[x];
            const uint16_t ne = prev_row[x + 1];
            const uint16_t w2 = row[x - 2];
            const uint16_t w = row[x - 1];
            const uint16_t value = row[x];

            //gap_dest[x] = ZigZag(value - GapPredict16(n2, ne2, nw, n, ne, w2, w));
            med_dest[x] = ZigZag(value - MedPredict16(n, w, nw));
        }

        row_data += stride_bytes;
        prev_row2 = prev_row;
        prev_row = row;
    }

    for (int i = 0; i < 1; ++i) {
        Shuffle(candidates[i], buffers[i]);
    }
}

// Must be kept in sync with above
void EncodeOuter(const uint8_t* image_data, int width, int pixel_offset_x, int pixel_offset_y, int stride_bytes, CandidateBuffer buffers[Predictor_Count])
{
    const uint8_t* row_data = image_data + pixel_offset_x * 2 + pixel_offset_y * stride_bytes;

    const uint16_t* prev_row = reinterpret_cast<const uint16_t*>( row_data - stride_bytes );
    const uint16_t* prev_row2 = reinterpret_cast<const uint16_t*>( row_data - stride_bytes * 2 );

    CandidateBuffer candidates[Predictor_Count];

    for (int y = 0; y < kBlockSize; ++y)
    {
        const uint16_t* row = reinterpret_cast<const uint16_t*>( row_data );
        uint16_t* gap_dest = reinterpret_cast<uint16_t*>( candidates[Predictor_GAP].Buffer + y * kBlockSize * 2 );
        uint16_t* med_dest = reinterpret_cast<uint16_t*>( candidates[Predictor_MED].Buffer + y * kBlockSize * 2 );

        for (int x = 0; x < kBlockSize; ++x)
        {
            /*
                        n2  ne2
                    nw  n   ne
                w2  w   ?
            */
            const uint16_t value = row[x];

            //------------------------------------------
            // This part is unique to the Outer version
            //------------------------------------------
            uint16_t n2 = 0, ne2 = 0, nw = 0, n = 0, ne = 0, w2 = 0, w = 0;

            const int image_x = x + pixel_offset_x;
            const int image_y = y + pixel_offset_y;

            if (image_x >= 2) {
                w2 = row[x - 2];
            }
            if (image_x >= 1) {
                w = row[x - 1];
            }
            if (image_y >= 1) {
                if (image_x >= 1) {
                    nw = prev_row[x - 1];
                }
                n = prev_row[x];
                if (image_x + 1 < width) {
                    ne = prev_row[x + 1];
                }
            }
            if (image_y >= 2) {
                n2 = prev_row2[x];
                if (image_x + 1 < width) {
                    ne2 = prev_row2[x + 1];
                }
            }
            //------------------------------------------
            // This part is unique to the Outer version
            //------------------------------------------

            //gap_dest[x] = ZigZag(value - GapPredict16(n2, ne2, nw, n, ne, w2, w));
            med_dest[x] = ZigZag(value - MedPredict16(n, w, nw));
        }

        row_data += stride_bytes;
        prev_row2 = prev_row;
        prev_row = row;
    }

    for (int i = 0; i < 1; ++i) {
        Shuffle(candidates[i], buffers[i]);
    }
}


//------------------------------------------------------------------------------
// Bit Accumulator

void BitAccumulator::Initialize(int width, int height)
{
    const int width_blocks = (width + kBlockSize - 1) / kBlockSize;
    const int height_blocks = (height + kBlockSize - 1) / kBlockSize;

    int bit_count = width_blocks * height_blocks;
    int byte_count = (bit_count + 7) / 8;
    PredictorChoice.resize(byte_count);
    BitOffset = 0;
}

void BitAccumulator::Append(int choice)
{
    PredictorChoice[BitOffset / 8] |= (uint8_t)choice << (BitOffset % 8);
    ++BitOffset;
}


//------------------------------------------------------------------------------
// Zstd Residual Context

bool ZstdResidualContext::Initialize()
{
    Context = ZSTD_createCCtx();
    if (!Context) {
        return false;
    }

    ZSTD_CCtx_setParameter(Context, ZSTD_c_compressionLevel, kZstdCompressLevel);
    ZSTD_CCtx_setParameter(Context, ZSTD_c_strategy, ZSTD_dfast);
    ZSTD_CCtx_setParameter(Context, ZSTD_c_enableLongDistanceMatching, 1);
    ZSTD_CCtx_setParameter(Context, ZSTD_c_hashLog, 6);
    ZSTD_CCtx_setParameter(Context, ZSTD_c_chainLog, 6);
    ZSTD_CCtx_setParameter(Context, ZSTD_c_windowLog, 16);
    //ZSTD_CCtx_setParameter(Context, ZSTD_c_searchLog, 4);
    //ZSTD_CCtx_setParameter(Context, ZSTD_c_targetLength, 8);
    //ZSTD_CCtx_setParameter(Context, ZSTD_c_minMatch, 9);
    ZSTD_CCtx_setParameter(Context, ZSTD_c_nbWorkers, 4);

    return true;
}

void ZstdResidualContext::Shutdown()
{
    if (Context) {
        ZSTD_freeCCtx(Context);
        Context = nullptr;
    }
}

size_t ZstdResidualContext::Compress(
    void* dst, size_t dstCapacity,
    const void* src, size_t srcSize)
{
    size_t result = ZSTD_compress2(Context, dst, dstCapacity, src, srcSize);

    if (ZSTD_isError(result)) {
        return 0;
    }

    return result;
}


//------------------------------------------------------------------------------
// Decoding

void DecodeGapBlockInner(
    uint8_t* dest_row_data, int dest_stride_bytes,
    const uint8_t* residuals_data, int plane_bytes)
{
    static_assert(kSliceCount == 4, "FIXME");

    uint16_t* row = reinterpret_cast<uint16_t*>( dest_row_data );
    const uint16_t* prev_row = reinterpret_cast<const uint16_t*>( dest_row_data - dest_stride_bytes );
    const uint16_t* prev_row2 = reinterpret_cast<const uint16_t*>( dest_row_data - dest_stride_bytes * 2 );

    for (int x = 0; x < kBlockSize; x++)
    {
        {
            uint16_t w0 = residuals_data[0] & 15;
            w0 |= ((uint16_t)residuals_data[plane_bytes] & 15) << 4;
            w0 |= ((uint16_t)residuals_data[plane_bytes * 2] & 15) << 8;
            w0 |= ((uint16_t)residuals_data[plane_bytes * 3] & 15) << 12;
            int16_t value = UnZigZag(w0);

            const uint16_t n2 = prev_row2[x];
            const uint16_t ne2 = prev_row2[x + 1];
            const uint16_t nw = prev_row[x - 1];
            const uint16_t n = prev_row[x];
            const uint16_t ne = prev_row[x + 1];
            const uint16_t w2 = row[x - 2];
            const uint16_t w = row[x - 1];

            row[x] = value + GapPredict16(n2, ne2, nw, n, ne, w2, w);
        }

        ++x;

        {
            uint16_t w1 = residuals_data[0] >> 4;
            w1 |= ((uint16_t)residuals_data[plane_bytes] >> 4) << 4;
            w1 |= ((uint16_t)residuals_data[plane_bytes * 2] >> 4) << 8;
            w1 |= ((uint16_t)residuals_data[plane_bytes * 3] >> 4) << 12;
            int16_t value = UnZigZag(w1);

            const uint16_t n2 = prev_row2[x];
            const uint16_t ne2 = prev_row2[x + 1];
            const uint16_t nw = prev_row[x - 1];
            const uint16_t n = prev_row[x];
            const uint16_t ne = prev_row[x + 1];
            const uint16_t w2 = row[x - 2];
            const uint16_t w = row[x - 1];

            row[x] = value + GapPredict16(n2, ne2, nw, n, ne, w2, w);
        }

        ++residuals_data;
    }
}

void DecodeGapBlockOuter(
    uint8_t* dest_row_data,
    int width, int pixel_offset_x, int pixel_offset_y, int dest_stride_bytes,
    const uint8_t* residuals_data, int plane_bytes)
{
    static_assert(kSliceCount == 4, "FIXME");

    uint16_t* row = reinterpret_cast<uint16_t*>( dest_row_data );
    const uint16_t* prev_row = reinterpret_cast<const uint16_t*>( dest_row_data - dest_stride_bytes );
    const uint16_t* prev_row2 = reinterpret_cast<const uint16_t*>( dest_row_data - dest_stride_bytes * 2 );

    for (int x = 0; x < kBlockSize; ++x)
    {
        {
            uint16_t w0 = residuals_data[0] & 15;
            w0 |= ((uint16_t)residuals_data[plane_bytes] & 15) << 4;
            w0 |= ((uint16_t)residuals_data[plane_bytes * 2] & 15) << 8;
            w0 |= ((uint16_t)residuals_data[plane_bytes * 3] & 15) << 12;
            int16_t value = UnZigZag(w0);

            uint16_t n2 = 0, ne2 = 0, nw = 0, n = 0, ne = 0, w2 = 0, w = 0;

            const int image_x = x + pixel_offset_x;
            const int image_y = pixel_offset_y;

            if (image_x > 1) {
                w2 = row[x - 2];
            }
            if (image_x > 0) {
                w = row[x - 1];
            }
            if (image_y > 0) {
                if (image_x > 0) {
                    nw = prev_row[x - 1];
                }
                n = prev_row[x];
                if (image_x + 1 < width) {
                    ne = prev_row[x + 1];
                }
            }
            if (image_y > 1) {
                n2 = prev_row2[x];
                if (image_x + 1 < width) {
                    ne2 = prev_row2[x + 1];
                }
            }

            row[x] = value + GapPredict16(n2, ne2, nw, n, ne, w2, w);
        }

        ++x;

        {
            uint16_t w1 = residuals_data[0] >> 4;
            w1 |= ((uint16_t)residuals_data[plane_bytes] >> 4) << 4;
            w1 |= ((uint16_t)residuals_data[plane_bytes * 2] >> 4) << 8;
            w1 |= ((uint16_t)residuals_data[plane_bytes * 3] >> 4) << 12;
            int16_t value = UnZigZag(w1);

            uint16_t n2 = 0, ne2 = 0, nw = 0, n = 0, ne = 0, w2 = 0, w = 0;

            const int image_x = x + pixel_offset_x;
            const int image_y = pixel_offset_y;

            if (image_x > 1) {
                w2 = row[x - 2];
            }
            if (image_x > 0) {
                w = row[x - 1];
            }
            if (image_y > 0) {
                if (image_x > 0) {
                    nw = prev_row[x - 1];
                }
                n = prev_row[x];
                if (image_x + 1 < width) {
                    ne = prev_row[x + 1];
                }
            }
            if (image_y > 1) {
                n2 = prev_row2[x];
                if (image_x + 1 < width) {
                    ne2 = prev_row2[x + 1];
                }
            }

            row[x] = value + GapPredict16(n2, ne2, nw, n, ne, w2, w);
        }

        ++residuals_data;
    }
}

void DecodeMedBlockInner(
    uint8_t* dest_row_data, int dest_stride_bytes,
    const uint8_t* residuals_data, int plane_bytes)
{
    static_assert(kSliceCount == 4, "FIXME");

    uint16_t* row = reinterpret_cast<uint16_t*>( dest_row_data );
    const uint16_t* prev_row = reinterpret_cast<const uint16_t*>( dest_row_data - dest_stride_bytes );

    for (int x = 0; x < kBlockSize; ++x)
    {
        {
            uint16_t w0 = residuals_data[0] & 15;
            w0 |= ((uint16_t)residuals_data[plane_bytes] & 15) << 4;
            w0 |= ((uint16_t)residuals_data[plane_bytes * 2] & 15) << 8;
            w0 |= ((uint16_t)residuals_data[plane_bytes * 3] & 15) << 12;
            int16_t value = UnZigZag(w0);

            const uint16_t nw = prev_row[x - 1];
            const uint16_t n = prev_row[x];
            const uint16_t w = row[x - 1];

            row[x] = value + MedPredict16(n, w, nw);
        }

        ++x;

        {
            uint16_t w1 = residuals_data[0] >> 4;
            w1 |= ((uint16_t)residuals_data[plane_bytes] >> 4) << 4;
            w1 |= ((uint16_t)residuals_data[plane_bytes * 2] >> 4) << 8;
            w1 |= ((uint16_t)residuals_data[plane_bytes * 3] >> 4) << 12;
            int16_t value = UnZigZag(w1);

            const uint16_t nw = prev_row[x - 1];
            const uint16_t n = prev_row[x];
            const uint16_t w = row[x - 1];

            row[x] = value + MedPredict16(n, w, nw);
        }

        ++residuals_data;
    }
}

void DecodeMedBlockOuter(
    uint8_t* dest_row_data,
    int pixel_offset_x, int pixel_offset_y, int dest_stride_bytes,
    const uint8_t* residuals_data, int plane_bytes)
{
    static_assert(kSliceCount == 4, "FIXME");

    uint16_t* row = reinterpret_cast<uint16_t*>( dest_row_data );
    const uint16_t* prev_row = reinterpret_cast<const uint16_t*>( dest_row_data - dest_stride_bytes );

    for (int x = 0; x < kBlockSize; x++)
    {
        {
            uint16_t w0 = residuals_data[0] & 15;
            w0 |= ((uint16_t)residuals_data[plane_bytes] & 15) << 4;
            w0 |= ((uint16_t)residuals_data[plane_bytes * 2] & 15) << 8;
            w0 |= ((uint16_t)residuals_data[plane_bytes * 3] & 15) << 12;
            int16_t value = UnZigZag(w0);

            uint16_t nw = 0, n = 0, w = 0;
            int image_x = x + pixel_offset_x;
            int image_y = pixel_offset_y;
            if (image_y > 0) {
                if (image_x > 0) {
                    nw = prev_row[x - 1];
                }
                n = prev_row[x];
            }
            if (image_x > 0) {
                w = row[x - 1];
            }

            row[x] = value + MedPredict16(n, w, nw);
        }

        ++x;

        {
            uint16_t w1 = residuals_data[0] >> 4;
            w1 |= ((uint16_t)residuals_data[plane_bytes] >> 4) << 4;
            w1 |= ((uint16_t)residuals_data[plane_bytes * 2] >> 4) << 8;
            w1 |= ((uint16_t)residuals_data[plane_bytes * 3] >> 4) << 12;
            int16_t value = UnZigZag(w1);

            uint16_t nw = 0, n = 0, w = 0;
            int image_x = x + pixel_offset_x;
            int image_y = pixel_offset_y;
            if (image_y > 0) {
                nw = prev_row[x - 1];
                n = prev_row[x];
            }
            if (image_x > 0) {
                w = row[x - 1];
            }

            row[x] = value + MedPredict16(n, w, nw);
        }

        ++residuals_data;
    }
}


} // namespace z16
