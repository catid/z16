# Z16: 16-bit Monochrome Image Compressor

Roughly 2x faster to encode/decode than PNG, while providing smaller file sizes.

## Quick Guide

Build and run the tool from source.

Compressing `test.png` to `test.z16`:

```
$ ./z16_tool.exe -i test.png -o test.z16 -v -f
Image read in 474.741 msec: test.png
Image encoded in 779.719 msec: test.z16
Image written to file in 31.449 msec: test.z16

$ ./z16_tool.exe -i test.z16 -o test2.png -v -f
Decode in 256.454 msec: test.z16
Write PNG output in 1784.75 msec: test2.png
```

The resulting file is 10% smaller than the original PNG, while [de]compressing faster.

## SDK Guide

Link to the z16 library.
Add to your source code: `#include <z16.h>`
Take a look at the app/z16_tool.cpp software to see how to use the simple API.

Simplified example:

```
    z16_Codec codec = z16_Create();
    if (!codec) {
        // error handling
    }
    ScopedFunction codec_scope([&]() {
        z16_Destroy(&codec);
    });

    z16_Frame frame;
    frame.Format = z16_Format_Mono16;
    frame.HeightPixels = height;
    frame.Image = data;
    frame.StrideBytes = stride;
    frame.WidthPixels = width;

    int32_t bytes = 0, error;
    uint8_t* data = z16_Encode(codec, &frame, &bytes, &error);
    if (!data) {
        // error handling
    }

    z16_Frame* decoded = z16_Decode(codec, data, bytes, &error);
    if (!decoded) {
        // error handling
    }
    ScopedFunction frame_scope([&]() {
        z16_Free(codec, &decoded);
    });
```
