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

float g_min(float a, float b) {
    return a > b ? b : a;
}

float g_max(float a, float b) {
    return a > b ? a : b;
}

float g_clamp(float val, float min, float max) {
    return g_min(g_max(val, min), max);
}

float g_median(float r, float g, float b) {
    return g_max(g_min(r, g), g_min(g_max(r, g), b));
}

float g_map(float min, float max, float v) {
    return (v - min) / (max - min);
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

    int glyph = 'Y';

    int glyphIdx = stbtt_FindGlyphIndex(&stbttInfo, glyph);

    msdf_AllocCtx allocCtx = {g_alloc, g_free, NULL};

    int borderSize = 4;

    msdf_Result result;
    int success = msdf_genGlyph(&result, &stbttInfo, glyphIdx, borderSize, genScale, 2.0f / genSize, &allocCtx);

    if (success == 0) {
        assert(!"Failed to generate msdf glyph");
        return 1;
    }

    FILE* fp = fopen(SAMPLE_ROOT "/sdf.png", "wb");

    uint8_t* pixels = malloc(sizeof(uint8_t) * result.width * result.height * 3);
    float scale = genSize;
    float maxValue = 1.0 * (scale);
    float transistionWidth = ((((((float) genSize) * 0.7f) + scale) / (scale * 2.0f)));
    float transistionAbs = (transistionWidth - 0.5);
    float transistStart = 0.5 - transistionAbs / 2;
    float transistEnd = transistStart + transistionAbs;
   
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

            if (r > transistStart) {
                if (r > (transistEnd)) {
                    r = 1.0f;
                } else {
                    r = 0.0f + (r - transistStart) / (transistionAbs);
                }
            } else {
                r = 0.0f;
            }
            if (g > transistStart) {
                if (g > (transistEnd)) {
                    g = 1.0f;
                } else {
                    g = 0.0f + (g - transistStart) / (transistionAbs);
                }
            } else {
                g = 0.0f;
            }
            if (b > transistStart) {
                if (b > (transistEnd)) {
                    b = 1.0f;
                } else {
                    b = 0.0f + (b - transistStart) / (transistionAbs);
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

    // apply sdf

    FILE* finalImage = fopen(SAMPLE_ROOT "/final.png", "wb");

    uint8_t* rgbaFinal = malloc(sizeof(uint8_t) * result.width * result.height * 4);

    for (int y = 0; y < result.height; y++) {
        uint8_t* pixelRow = pixels + (y * result.width * 3);
        uint8_t* rgbaFinalRow = rgbaFinal + (y * result.width * 4);
        for (int x = 0; x < result.width; x++) {
            float r = ((float)pixelRow[x * 3 + 0]) / 255.0f;// = r * 255.0f; // (r > 0.5f) ? 255.0f : r * 255.0f;
            float g = ((float)pixelRow[x * 3 + 1]) / 255.0f;// = g * 255.0f; // (g > 0.5f) ? 255.0f : g * 255.0f;
            float b = ((float)pixelRow[x * 3 + 2]) / 255.0f;// = b * 255.0f; // (b > 0.5f) ? 255.0f : b * 255.0f;
            float dist = g_median(r, g, b);
            float opacity = g_clamp(dist - 0.5f, 0, 1) * 2.0f;
            rgbaFinalRow[x * 4 + 0] = 0;//255.0f;
            rgbaFinalRow[x * 4 + 1] = 0;//255.0f;
            rgbaFinalRow[x * 4 + 2] = 0;//255.0f;
            rgbaFinalRow[x * 4 + 3] = opacity * 255.0f;

        }
    }
    svpng(finalImage, result.width, result.height, rgbaFinal, 1);
    fclose(finalImage);
    return 0;
}