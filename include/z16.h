// Copyright (c) 2020, Christopher A. Taylor.  All rights reserved.

/**
 * Z16 :: 16-bit Monochrome Image Compressor
 *
 * Our goals are to provide real-time (fast) compression for 16-bit RAW images
 * with a faster codec and higher compression than existing formats like TIFF.
 */

#ifndef Z16_H
#define Z16_H

#include <stdint.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#define Z16_CPP(x) x
#else
#define Z16_CPP(x)
#endif


//------------------------------------------------------------------------------
// Constants

// Version of the header
#define Z16_VERSION 100 /* 1.0.0 */

/*
    Release Notes:

    1.0: First release
*/

/*
    File Format:

    Header (16 bytes):

        0xFE 0xAD
        <version=1 (8 bits)>
        <enum z16_Format (8 bits)>
        <width-1 (16 bits)>
        <height-1 (16 bits)>
        <choices bytes (32 bits)> = 0 if just storing original data
        <residuals bytes (32 bits)>

    Body:

        <compressed choices (x)>
        <compressed residuals (x)>
*/


//------------------------------------------------------------------------------
// Datatypes

typedef struct z16_Codec_t { int32_t impl; }* z16_Codec;

typedef enum z16_Error_t
{
    z16_Success = 0,

    z16_Error_InvalidInput = 1,
    z16_Error_Zstd = 2,
    z16_Error_OOM = 3,
} z16_Error;

typedef enum z16_Format_t
{
    z16_Format_Invalid,

    /*
        Monochrome input:
        The input is considered a single plane.
    */
    z16_Format_Mono16,

    z16_Format_Count
} z16_Format;

typedef struct z16_Frame_t
{
    /// Format of the image. One of the enum z16_Format constants
    int32_t Format Z16_CPP(= z16_Format_Invalid);

    /// Pointer to image
    uint8_t* Image Z16_CPP(= nullptr);

    /// Width of each plane in pixels
    int32_t WidthPixels Z16_CPP(= 0);

    /// Number of bytes per row for each plane
    int32_t StrideBytes Z16_CPP(= 0);

    /// Height of each plane in pixels
    int32_t HeightPixels Z16_CPP(= 0);
} z16_Frame;


//------------------------------------------------------------------------------
// API

/// Create a codec context.
/// Returns 0 on failure.
extern z16_Codec z16_Create();

/// Shutdown and free codec.
/// All decoded frames must be freed using z16_Free() prior to this call.
/// Also sets the codec pointer to null.
extern void z16_Destroy(z16_Codec* codec);

/*
    Encode a frame.

    Returns a buffer pointer on encode success, and sets `bytes`
    to the number of bytes provided.
    The buffer pointer will be valid until the next call to Encode.

    Returns 0 on encode failure; error will be set to an error code.
*/
extern uint8_t* z16_Encode(
    z16_Codec codec,
    const z16_Frame* frame,
    int32_t* bytes,
    int32_t* error);

/*
    Decode a frame from compressed buffer, including metadata.

    The returned frame must be freed with z16_Free(),
    prior to calling z16_Destroy().
    Sets the number of bytes used in the `used_bytes` parameter.

    Returns 0 on decode failure; error will be set to an error code.
*/
extern z16_Frame* z16_Decode(
    z16_Codec codec,
    const void* data,
    int32_t bytes,
    int32_t* error);

/// Must be called by application for each successful z16_Decode call.
/// Sets the given frame pointer to null.
extern void z16_Free(z16_Codec codec, z16_Frame** frame);


#ifdef __cplusplus
} // extern "C"
#endif


#endif // Z16_H
