// Copyright (c) 2020, Christopher A. Taylor.  All rights reserved.

#pragma once

#include "Tools.hpp"

namespace z16 {


//------------------------------------------------------------------------------
// Codec

class Codec
{
public:
    ~Codec()
    {
        Shutdown();
    }

    bool Initialize();
    void Shutdown();

    z16_Frame* Decode(const uint8_t* data, int32_t bytes, int32_t& error);
    void Free(z16_Frame* frame);

    uint8_t* Encode(const z16_Frame& frame, int32_t& bytes, int32_t& error);

protected:
    // Encoder:
    PredictorPicker Picker;
    EncodeBuffer EncodedResiduals;
    BitAccumulator EncodedChoices;
    ZstdResidualContext ZstdResiduals;
    std::vector<uint8_t> Output;

protected:
    // Decoder:
    std::vector<z16_Frame*> FreedFrames;
    std::vector<uint8_t> Choices;
    std::vector<uint8_t> Residuals;

    z16_Frame* Allocate(int32_t format, int width, int height);
    void FreeAll();
};


} // namespace z16
