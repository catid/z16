// Copyright (c) 2020, Christopher A. Taylor.  All rights reserved.

#include "z16.h"

#include "argagg.hpp"
#include "teenypath.h"
#include "lodepng.h"

#include <chrono>
#include <functional>
#include <iostream>
#include <fstream>
using namespace std;


//------------------------------------------------------------------------------
// Tools

static bool m_Verbose = false;

// Convert string to upper case
static std::string ToUpperCase(const std::string& s)
{
    std::string result = s;
    std::transform(s.begin(), s.end(), result.begin(),
        [](unsigned char c) { return (char)std::toupper(c); });
    return result;
}

static bool IsZ16Filename(const std::string& filename)
{
    TeenyPath::path filepath(filename);
    std::string ext = ToUpperCase(filepath.extension());
    return ext == ".Z16";
}

uint64_t GetUnixTimeNowUsec()
{
    auto nowTsMicrosDuration =
        std::chrono::duration_cast<std::chrono::microseconds>(
            std::chrono::system_clock::now().time_since_epoch());
    int64_t imageTsMicrosNow = nowTsMicrosDuration.count();
    return imageTsMicrosNow;
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
// Encode

static bool Encode(const std::string& input, const std::string& output)
{
    z16_Codec codec = z16_Create();
    if (!codec) {
        cerr << "z16_Create failed: " << input << endl;
        return false;
    }
    ScopedFunction codec_scope([&]() {
        z16_Destroy(&codec);
    });

    uint64_t t0 = GetUnixTimeNowUsec();

    std::vector<uint8_t> image;
    unsigned w = 0, h = 0;
    unsigned png_error = lodepng::decode(image, w, h, input, LCT_GREY, 16);

    if (png_error) {
        cerr << "Failed to load image: " << input << " error = " << lodepng_error_text(png_error) << endl;
        return false;
    }

    for (int i = 0; i < image.size(); i += 2) {
        uint8_t temp = image[i];
        image[i] = image[i + 1];
        image[i + 1] = temp;
    }

    const int32_t format = z16_Format_Mono16;

    uint64_t t1 = GetUnixTimeNowUsec();
    if (m_Verbose) {
        cout << "Image read in " << (t1 - t0) / 1000.f << " msec: " << input << endl;
    }

    z16_Frame frame{};

    // Set frame format
    frame.Format = format;
    switch (format)
    {
    case z16_Format_Mono16:
        frame.WidthPixels = w;
        frame.StrideBytes = w * 2;
        frame.HeightPixels = h;
        break;
    default:
        cerr << "FIXME: Add support for format" << endl;
        return false;
    }

    frame.Image = image.data();

    int32_t bytes = 0, error;
    uint8_t* data = z16_Encode(codec, &frame, &bytes, &error);
    if (!data) {
        cerr << "z16_Encode failed: " << input << " error = " << error << " w=" << w << " h=" << h << endl;
        return false;
    }

    uint64_t t2 = GetUnixTimeNowUsec();
    if (m_Verbose) {
        cout << "Image encoded in " << (t2 - t1) / 1000.f << " msec: " << output << endl;
    }

    std::ofstream file(output, std::ifstream::binary);
    if (!file) {
        cerr << "Failed to write output file: " << output << endl;
        return false;
    }
    file.write((const char*)data, bytes);

    uint64_t t3 = GetUnixTimeNowUsec();
    if (m_Verbose) {
        cout << "Image written to file in " << (t3 - t2) / 1000.f << " msec: " << output << endl;
    }

    return true;
}


//------------------------------------------------------------------------------
// Decode

static bool Decode(const std::string& input, const std::string& output)
{
    uint64_t t0 = GetUnixTimeNowUsec();

    std::vector<uint8_t> file_data;
    std::ifstream file(input, std::ifstream::ate | std::ifstream::binary);
    if (!file) {
        cerr << "Failed to read input file: " << input << endl;
        return false;
    }
    file_data.resize(file.tellg());
    file.seekg(0);
    file.read((char*)file_data.data(), file_data.size());

    z16_Codec codec = z16_Create();
    if (!codec) {
        cerr << "z16_Create failed: " << input << endl;
        return false;
    }
    ScopedFunction codec_scope([&]() {
        z16_Destroy(&codec);
    });

    int32_t error;
    z16_Frame* frame = z16_Decode(
        codec,
        file_data.data(),
        (int32_t)file_data.size(),
        &error);
    if (!frame) {
        cerr << "z16_Decode failed err=" << error << ": " << input << endl;
        return false;
    }

    uint64_t t1 = GetUnixTimeNowUsec();
    if (m_Verbose) {
        cout << "Decode in " << (t1 - t0) / 1000.f << " msec: " << input << endl;
    }

    // Only supports 16-bit monochrome format
    if (frame->Format != z16_Format_Mono16) {
        cerr << "FIXME: Unsupported format" << endl;
        return false;
    }

    const int w = frame->WidthPixels;
    const int h = frame->HeightPixels;

    for (int i = 0; i < frame->WidthPixels * frame->HeightPixels * 2; i += 2) {
        uint8_t temp = frame->Image[i];
        frame->Image[i] = frame->Image[i + 1];
        frame->Image[i + 1] = temp;
    }

    unsigned png_error = lodepng::encode(
        output, frame->Image,
        frame->WidthPixels, frame->HeightPixels,
        LCT_GREY, 16);
    if (png_error) {
        cerr << "Failed to write 16-bit PNG output: " << output << endl;
        return false;
    }

    uint64_t t2 = GetUnixTimeNowUsec();
    if (m_Verbose) {
        cout << "Write PNG output in " << (t2 - t1) / 1000.f << " msec: " << output << endl;
    }

    return true;
}


//------------------------------------------------------------------------------
// Entrypoint

static void banner()
{
    cout << "Z16: 16-bit monochrome image compressor V" << Z16_VERSION << endl;
    cout << "Compress example:" << endl;
    cout << "    z16_tool -i source.png -o dest.z16 -v -f" << endl;
    cout << "Decompress example:" << endl;
    cout << "    z16_tool -i source.z16 -o dest.png -v -f" << endl;
}

int main(int argc, char *argv[])
{
    // Define arguments
    argagg::parser argparser{{
        {"help", {"-h", "--help"},
            "Shows this help message.", 0},
        {"input", {"-i", "--input"},
            "Input image", 1},
        {"output", {"-o", "--output"},
            "Output image", 1},
        {"verbose", {"-v", "--verbose"},
            "Log extra information", 0},
        {"force", {"-f", "--force"},
            "Force overwrite existing output file", 0},
    }};

    // Parse arguments
    argagg::parser_results args;
    try {
        args = argparser.parse(argc, argv);
    } catch (const std::exception &e) {
        cout << argparser << endl;
        cerr << "Invalid arguments: " << e.what() << endl;
        return -1;
    }

    if (args["help"]) {
        banner();
        cout << argparser << endl;
        return 0;
    }

    std::string input = args["input"].as<std::string>("");
    std::string output = args["output"].as<std::string>("");
    m_Verbose = args["verbose"];
    const bool force = args["force"];

    if (input.empty()) {
        banner();
        cout << argparser << endl;
        cerr << "Please specify an input file: -i filename" << endl;
        return -1;
    }
    const bool is_z16_input = IsZ16Filename(input);

    if (output.empty()) {
        TeenyPath::path input_path(input);
        if (is_z16_input) {
            // The PNG format is a good default.
            // It supports compression for 16-bit images.
            // This is very inefficient, so the goal here is to convert to
            // an interchange format for previewing and not for heavy lifting.
            input_path.replace_extension(".png");
        } else {
            input_path.replace_extension(".z16");
        }
        if (input_path.exists()) {
            if (!force) {
                cerr << "Default output path already exists: " << input_path.native_string() << endl;
                return -1;
            } else if (m_Verbose) {
                cerr << "Overwriting output file that already exists: " << input_path.native_string() << endl;
            }
        }
        output = input_path.native_string();
    } else {
        TeenyPath::path input_path(output);
        if (input_path.exists()) {
            if (!force) {
                cerr << "Output path already exists: " << input_path.native_string() << endl;
                return -1;
            } else if (m_Verbose) {
                cerr << "Overwriting output file that already exists: " << input_path.native_string() << endl;
            }
        }
    }

    const bool is_z16_output = IsZ16Filename(output);

    if (!is_z16_input && !is_z16_output) {
        cerr << "Either input or output file must have extension: .z16" << endl;
        return -1;
    }
    if (is_z16_input && is_z16_output) {
        cerr << "Either input or output file must not be an .z16 file" << endl;
        return -1;
    }

    const bool is_writing_z16 = is_z16_output;

    bool success = false;
    if (is_writing_z16) {
        success = Encode(input, output);
    } else {
        success = Decode(input, output);
    }

    if (!success) {
        cerr << "Operation failed" << endl;
        return -1;
    }

    return 0;
}
