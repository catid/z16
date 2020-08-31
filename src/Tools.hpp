// Copyright (c) 2020, Christopher A. Taylor.  All rights reserved.

/*
    References:

    [1] "CALIC - A Context Based Adaptive Lossless Image Codec" (1996)
    Xiaolin Wu, Nassr Memon

        This provides the GAP (Gradient Adjusted Predictor).

    [2] "The LOCO-I Lossless Image Compression Algorithm: Principles and Standardization into JPEG-LS"
    https://www.hpl.hp.com/loco/HPL-98-193R1.pdf

        Median Edge Detector (MED) Predictor

        The formulation is different but it's identical to this one:
            FFmpeg algorithm https://en.wikipedia.org/wiki/FFV1
            FFV1 Lossless Video Codec Predictor (also used in FLIF)
*/

#pragma once

#include <zstd.h>

#include <stdint.h>
#include <cstring>
#include <algorithm>
#include <vector>

#include "z16.h"

namespace z16 {


//------------------------------------------------------------------------------
// Constants

static const int kBlockSize = 16; // 16x16 blocks

static const int kSliceCount = 4; // 16/4 bits per slice

static const int kZstdCompressLevel = 1;

enum Predictor_t
{
    Predictor_MED, // [1]
    Predictor_GAP, // [2]

    Predictor_Count
};

static const int kFileHeaderBytes = 16;

static const uint8_t kHeaderByte0 = 0xfe; // magic
static const uint8_t kHeaderByte1 = 0xad; // magic
static const uint8_t kHeaderVersion = 0x01;


//------------------------------------------------------------------------------
// PredictorPicker

struct CandidateBuffer
{
    uint8_t Buffer[kBlockSize * kBlockSize * 2];
};

struct PredictorPicker
{
    static const int kByteLevels = 256;

    // Total pixels accumulated so far by PickPredictor()
    uint32_t SliceTotal = 0;

    // Count of pixels accumulated so far at each level
    uint32_t ImageCounters[kSliceCount][kByteLevels];

    // Accumulator for PixelAccumulate()
    uint32_t BlockCounters[Predictor_Count][kSliceCount][kByteLevels];


    // Resets the ImageCounters and SliceTotal
    void ImageReset();

    // Returns enum AicPredictor value
    int PickPredictor(const CandidateBuffer buffers[Predictor_Count]);
};


//------------------------------------------------------------------------------
// EncodeBuffer

struct EncodeBuffer
{
    std::vector<uint8_t> Output;

    // Current write offset
    int WriteOffset = 0;

    // Bytes for each plane
    int PlaneBytes = 0;

    void Initialize(int width, int height);
    void Append(const CandidateBuffer& buffer);
};


//------------------------------------------------------------------------------
// Encoder

// Input size: kBlockSize * kBlockSize * 2
void Shuffle(const CandidateBuffer& source_buffer, CandidateBuffer& dest_buffer);

// Encode a block to a CandidateBuffer
void EncodeInner(const uint8_t* image_data, int pixel_offset_x, int pixel_offset_y, int stride_bytes, CandidateBuffer buffers[Predictor_Count]);
void EncodeOuter(const uint8_t* image_data, int width, int pixel_offset_x, int pixel_offset_y, int stride_bytes, CandidateBuffer buffers[Predictor_Count]);


//------------------------------------------------------------------------------
// Bit Accumulator

struct BitAccumulator
{
    std::vector<uint8_t> PredictorChoice;
    int BitOffset = 0;

    void Initialize(int width, int height);
    void Append(int choice);
};


//------------------------------------------------------------------------------
// Zstd Residual Context

struct ZstdResidualContext
{
    ZSTD_CCtx* Context = nullptr;

    ~ZstdResidualContext()
    {
        Shutdown();
    }

    bool Initialize();
    void Shutdown();

    size_t Compress(
        void* dst, size_t dstCapacity,
        const void* src, size_t srcSize);
};


//------------------------------------------------------------------------------
// Decoding

// GAP (Gradient Adjusted Predictor) from [1]
void DecodeGapBlockInner(
    uint8_t* dest_row_data, int dest_stride_bytes,
    const uint8_t* residuals_data, int plane_bytes);
void DecodeGapBlockOuter(
    uint8_t* dest_row_data,
    int width, int pixel_offset_x, int pixel_offset_y, int dest_stride_bytes,
    const uint8_t* residuals_data, int plane_bytes);

// Median Edge Detector (MED) Predictor from [2]
void DecodeMedBlockInner(
    uint8_t* dest_row_data, int dest_stride_bytes,
    const uint8_t* residuals_data, int plane_bytes);
void DecodeMedBlockOuter(
    uint8_t* dest_row_data,
    int pixel_offset_x, int pixel_offset_y, int dest_stride_bytes,
    const uint8_t* residuals_data, int plane_bytes);


} // namespace z16
