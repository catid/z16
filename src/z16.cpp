// Copyright (c) 2020, Christopher A. Taylor.  All rights reserved.

#include "z16.h"

#include <new> // std::nothrow
#include <iostream>
using namespace std;

#include "Codec.hpp"
using namespace z16;


//------------------------------------------------------------------------------
// C API Functions

extern "C" {


z16_Codec z16_Create()
{
    Codec* codec = new (std::nothrow) Codec;
    if (!codec) {
        return nullptr;
    }

    if (!codec->Initialize()) {
        codec->Shutdown();
        delete codec;
        return nullptr;
    }

    return reinterpret_cast<z16_Codec>( codec );
}

z16_Frame* z16_Decode(
    z16_Codec codec_,
    const void* data,
    int32_t bytes,
    int32_t* error)
{
    Codec* codec = reinterpret_cast<Codec*>( codec_ );
    if (!codec || !data || bytes <= 0) {
        cerr << "Null parameter";
        return nullptr;
    }

    int32_t error_result = z16_Success;
    z16_Frame* frame = codec->Decode((const uint8_t*)data, bytes, error_result);
    if (error) {
        *error = error_result;
    }
    return frame;
}

void z16_Free(z16_Codec codec_, z16_Frame** frame)
{
    Codec* codec = reinterpret_cast<Codec*>( codec_ );
    if (!codec || !frame || !*frame) {
        return; // Ignore invalid input
    }

    codec->Free(*frame);

    *frame = nullptr;
}

uint8_t* z16_Encode(
    z16_Codec codec_,
    const z16_Frame* frame,
    int32_t* bytes,
    int32_t* error)
{
    if (!frame || !bytes || !codec_) {
        cerr << "Null parameter";
        return nullptr;
    }

    Codec* codec = reinterpret_cast<Codec*>( codec_ );
    int32_t error_result = z16_Success;
    uint8_t* result = codec->Encode(*frame, *bytes, error_result);
    if (error) {
        *error = error_result;
    }
    return result;
}

void z16_Destroy(z16_Codec* codec_)
{
    if (!codec_ || !*codec_) {
        return;
    }

    Codec* codec = reinterpret_cast<Codec*>( *codec_ );

    codec->Shutdown();
    delete codec;

    *codec_ = nullptr;
}


} // extern "C"
