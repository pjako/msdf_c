#include "svpng.h"
#define STB_TRUETYPE_IMPLEMENTATION
#include "../stb_truetype.h"
#define MSDF_IMPLEMENTATION
#include "../msdf.h"

#include "stdio.h"


#ifndef SAMPLE_ROOT
#define SAMPLE_ROOT ""
#endif

typedef struct Range {
    void* content;
    size_t size;
} Range;

Range g_fileRead(const char* fileName) {
    Range range = {NULL, 0};

    FILE* fileHandle = fopen(fileName, "rb");
    if (fileHandle == NULL) {
        assert(!"File not fount");
        return range;
    }

    fseek(fileHandle, 0, SEEK_END); // seek to end of file
    size_t size = ftell(fileHandle); // get current file pointer
    fseek(fileHandle, 0, SEEK_SET); // seek back to beginning of file

    if (size == 0) {
        assert(!"File size zero");
        return range;
    }

    void* mem = malloc(size);
    size_t readSize = fread(mem, 1, size, fileHandle);
    if (readSize == size) {
        range.content = mem;
        range.size = size;
    } else {
        free(mem);
    }
    fclose(fileHandle);
    return range;
}

static void* g_alloc(size_t size, void* ctx) {
    return malloc(size);
}

static void g_free(void* ptr, void* ctx) {
    return free(ptr);
}

int main() {
    stbtt_fontinfo stbttInfo;

    Range fontFile = g_fileRead(SAMPLE_ROOT "/fonts/Roboto-Regular.ttf");

    if (fontFile.content == NULL) {
        assert(!"could not load font from disk");
        return 1;
    }

    if (!stbtt_InitFont(&stbttInfo, fontFile.content, stbtt_GetFontOffsetForIndex(fontFile.content, 0))) {
        assert(!"Couldn't parse font");
		return 1;
    }


    int ascent, descent;
    stbtt_GetFontVMetrics(&stbttInfo, &ascent, &descent, 0);

    int genSize = 128;
    float genScale = stbtt_ScaleForPixelHeight(&stbttInfo, genSize);

    int glyphIdx = stbtt_FindGlyphIndex(&stbttInfo, 'A');

    msdf_AllocCtx allocCtx = {g_alloc, g_free, NULL};

    int borderSize = 1;

    msdf_Result result;
    int success = msdf_genGlyph(&result, &stbttInfo, glyphIdx, borderSize, genScale, 2.0f / genSize, &allocCtx);

    if (success == 0) {
        assert(!"Failed to generate msdf glyph");
        return 1;
    }

    FILE* fp = fopen(SAMPLE_ROOT "a.png", "wb");

    uint8_t* pixels = malloc(sizeof(uint8_t) * result.width * result.height * 3);
    float scale = genSize;
    float maxValue = 1.0 * (scale);
    float transistionWidth = ((((((float) genSize) * 0.3f) + scale) / (scale * 2.0f)));
    //float ff = expf(0.01);
    for (int y = 0; y < result.height; y++) {
        int yPos = result.width * 3 * y;
        uint8_t* pixelRow = pixels + (y * result.width * 3);
        for (int x = 0; x < result.width; x++) {
            int indexSdf = yPos + (x * 3);
            float r = result.rgb[indexSdf + 0];
            float g = result.rgb[indexSdf + 1];
            float b = result.rgb[indexSdf + 2];


            r = ((((r) + scale) / (scale * 2.0f)));
            g = ((((g) + scale) / (scale * 2.0f)));
            b = ((((b) + scale) / (scale * 2.0f)));

            if (r > 0.5f) {
                if (r > (transistionWidth)) {
                    r = 1.0f;
                } else {
                    r = 0.0f + (r - 0.5f) * (transistionWidth) * 10.0f;
                }
            } else {
                r = 0.0f;
            }
            if (g > 0.5f) {
                if (g > (transistionWidth)) {
                    g = 1.0f;
                } else {
                    g = 0.0f + (g - 0.5f) * (transistionWidth) * 10.0f;
                }
            } else {
                g = 0.0f;
            }
            if (b > 0.5f) {
                if (b > (transistionWidth)) {
                    b = 1.0f;
                } else {
                    b = 0.0f + (b - 0.5f) * (transistionWidth) * 10.0f;
                }
            } else {
                b = 0.0f;
            }

            pixelRow[x * 3 + 0] = r * 255.0f; // (r > 0.5f) ? 255.0f : r * 255.0f;
            pixelRow[x * 3 + 1] = g * 255.0f; // (g > 0.5f) ? 255.0f : g * 255.0f;
            pixelRow[x * 3 + 2] = b * 255.0f; // (b > 0.5f) ? 255.0f : b * 255.0f;
        }
    }
    
    svpng(fp, result.width, result.height, pixels, 0);
    fclose(fp);
    return 0;
}