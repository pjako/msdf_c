/* msdf
  Handles multi-channel signed distance field bitmap
  generation from given ttf (stb_truetype.h) font.
  https://github.com/exezin/msdf-c
  Depends on stb_truetype.h to load the ttf file.
  This is in an unstable state, ymmv.
  Based on the C++ implementation by Viktor Chlumský.
  https://github.com/Chlumsky/msdfgen
*/

#ifndef MSDF_H
#define MSDF_H

#include <inttypes.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif
	
typedef struct {
  int glyphIdx;
  int left_bearing;
  int advance;
  int ix0, ix1;
  int iy0, iy1;
} ex_metrics_t;

typedef struct msdf_Result {
  int glyphIdx;
  int left_bearing;
  int advance;
  float* rgb;
  int width;
  int height;
  int yOffset;
} msdf_Result;

typedef struct msdf_AllocCtx {
  void* (*alloc)(size_t size, void* ctx);
  // free is optional and will not be called if it is null (useful for area allocators that free everything at once)
  void (*free)(void* ptr, void* ctx);
  void* ctx;
} msdf_AllocCtx;

/*
  Generates a bitmap from the specified glyph index of a stbtt font

  Returned result is 1 for success or 0 in case of an error
 */
int msdf_genGlyph(msdf_Result* result, stbtt_fontinfo *font, int stbttGlyphIndex, uint32_t borderWidth, float scale, float range, msdf_AllocCtx* alloc);

#ifdef __cplusplus
}
#endif

#ifdef MSDF_IMPLEMENTATION
// pixel at (x, y) in bitmap (arr)
#define msdf_pixelAt(x, y, w, arr) ((msdf_Vec3){arr[(3 * (((y)*w) + x))], arr[(3 * (((y)*w) + x)) + 1], arr[(3 * (((y)*w) + x)) + 2]})

#define msdf_max(x, y) (((x) > (y)) ? (x) : (y))
#define msdf_min(x, y) (((x) < (y)) ? (x) : (y))

#define MSDF_INF -1e24
#define MSDF_EDGE_THRESHOLD 0.02

#ifndef MSDF_PI
#define MSDF_PI 3.14159265358979323846
#endif

typedef float msdf_Vec2[2];
typedef float msdf_Vec3[3];

typedef struct
{
    double dist;
    double d;
} msdf_signedDistance;

// the possible types:
// STBTT_vmove  = start of a contour
// STBTT_vline  = linear segment
// STBTT_vcurve = quadratic segment
// STBTT_vcubic = cubic segment
typedef struct
{
    int color;
    msdf_Vec2 p[4];
    int type;
} msdf_EdgeSegment;

// defines what color channel an edge belongs to
typedef enum
{
    msdf_edgeColor_black = 0,
    msdf_edgeColor_red = 1,
    msdf_edgeColor_green = 2,
    msdf_edgeColor_yellow = 3,
    msdf_edgeColor_blue = 4,
    msdf_edgeColor_magenta = 5,
    msdf_edgeColor_cyan = 6,
    msdf_edgeColor_white = 7
} msdf_edgeColor;

static double msdf_median(double a, double b, double c)
{
    return msdf_max(msdf_min(a, b), msdf_min(msdf_max(a, b), c));
}

static int msdf_nonZeroSign(double n)
{
    return 2 * (n > 0) - 1;
}

static double msdf_cross(msdf_Vec2 a, msdf_Vec2 b)
{
    return a[0] * b[1] - a[1] * b[0];
}

static void msdf_v2Scale(msdf_Vec2 r, msdf_Vec2 const v, float const s)
{
    int i;
    for (i = 0; i < 2; ++i)
        r[i] = v[i] * s;
}

static float msdf_v2MulInner(msdf_Vec2 const a, msdf_Vec2 const b)
{
    float p = 0.;
    int i;
    for (i = 0; i < 2; ++i)
        p += b[i] * a[i];
    return p;
}

static float msdf_v2Leng(msdf_Vec2 const v)
{
    return sqrtf(msdf_v2MulInner(v, v));
}

static void msdf_v2Norm(msdf_Vec2 r, msdf_Vec2 const v)
{
    float k = 1.0 / msdf_v2Leng(v);
    msdf_v2Scale(r, v, k);
}

static void msdf_v2Sub(msdf_Vec2 r, msdf_Vec2 const a, msdf_Vec2 const b)
{
    int i;
    for (i = 0; i < 2; ++i)
        r[i] = a[i] - b[i];
}

int msdf_solveQuadratic(double x[2], double a, double b, double c)
{
    if (fabs(a) < 1e-14)
    {
        if (fabs(b) < 1e-14)
        {
            if (c == 0)
                return -1;
            return 0;
        }
        x[0] = -c / b;
        return 1;
    }

    double dscr = b * b - 4 * a * c;
    if (dscr > 0)
    {
        dscr = sqrt(dscr);
        x[0] = (-b + dscr) / (2 * a);
        x[1] = (-b - dscr) / (2 * a);
        return 2;
    }
    else if (dscr == 0)
    {
        x[0] = -b / (2 * a);
        return 1;
    }
    else
    {
        return 0;
    }
}

int msdf_solveCubicNormed(double *x, double a, double b, double c)
{
    double a2 = a * a;
    double q = (a2 - 3 * b) / 9;
    double r = (a * (2 * a2 - 9 * b) + 27 * c) / 54;
    double r2 = r * r;
    double q3 = q * q * q;
    double A, B;
    if (r2 < q3)
    {
        double t = r / sqrt(q3);
        if (t < -1)
            t = -1;
        if (t > 1)
            t = 1;
        t = acos(t);
        a /= 3;
        q = -2 * sqrt(q);
        x[0] = q * cos(t / 3) - a;
        x[1] = q * cos((t + 2 * MSDF_PI) / 3) - a;
        x[2] = q * cos((t - 2 * MSDF_PI) / 3) - a;
        return 3;
    }
    else
    {
        A = -pow(fabs(r) + sqrt(r2 - q3), 1 / 3.);
        if (r < 0)
            A = -A;
        B = A == 0 ? 0 : q / A;
        a /= 3;
        x[0] = (A + B) - a;
        x[1] = -0.5 * (A + B) - a;
        x[2] = 0.5 * sqrt(3.) * (A - B);
        if (fabs(x[2]) < 1e-14)
            return 2;
        return 1;
    }
}

int msdf_solveCubic(double x[3], double a, double b, double c, double d)
{
    if (fabs(a) < 1e-14)
        return msdf_solveQuadratic(x, b, c, d);

    return msdf_solveCubicNormed(x, b / a, c / a, d / a);
}

void msdf_getOrtho(msdf_Vec2 r, msdf_Vec2 const v, int polarity, int allow_zero)
{
    double len = msdf_v2Leng(v);

    if (len == 0)
    {
        if (polarity)
        {
            r[0] = 0;
            r[1] = !allow_zero;
        }
        else
        {
            r[0] = 0;
            r[1] = -!allow_zero;
        }
        return;
    }

    if (polarity)
    {
        r[0] = -v[1] / len;
        r[1] = v[0] / len;
    }
    else
    {
        r[0] = v[1] / len;
        r[1] = -v[0] / len;
    }
}

int msdf_pixelClash(const msdf_Vec3 a, const msdf_Vec3 b, double threshold)
{
    int aIn = (a[0] > .5f) + (a[1] > .5f) + (a[2] > .5f) >= 2;
    int bIn = (b[0] > .5f) + (b[1] > .5f) + (b[2] > .5f) >= 2;
    if (aIn != bIn)
        return 0;
    if ((a[0] > .5f && a[1] > .5f && a[2] > .5f) || (a[0] < .5f && a[1] < .5f && a[2] < .5f) || (b[0] > .5f && b[1] > .5f && b[2] > .5f) || (b[0] < .5f && b[1] < .5f && b[2] < .5f))
        return 0;
    float aa, ab, ba, bb, ac, bc;
    if ((a[0] > .5f) != (b[0] > .5f) && (a[0] < .5f) != (b[0] < .5f))
    {
        aa = a[0], ba = b[0];
        if ((a[1] > .5f) != (b[1] > .5f) && (a[1] < .5f) != (b[1] < .5f))
        {
            ab = a[1], bb = b[1];
            ac = a[2], bc = b[2];
        }
        else if ((a[2] > .5f) != (b[2] > .5f) && (a[2] < .5f) != (b[2] < .5f))
        {
            ab = a[2], bb = b[2];
            ac = a[1], bc = b[1];
        }
        else
        {
            return 0;
        }
    }
    else if ((a[1] > .5f) != (b[1] > .5f) && (a[1] < .5f) != (b[1] < .5f) && (a[2] > .5f) != (b[2] > .5f) && (a[2] < .5f) != (b[2] < .5f))
    {
        aa = a[1], ba = b[1];
        ab = a[2], bb = b[2];
        ac = a[0], bc = b[0];
    }
    else
    {
        return 0;
    }
    return (fabsf(aa - ba) >= threshold) && (fabsf(ab - bb) >= threshold) && fabsf(ac - .5f) >= fabsf(bc - .5f);
}

void msdf_mix(msdf_Vec2 r, msdf_Vec2 a, msdf_Vec2 b, double weight)
{
    r[0] = (1 - weight) * a[0] + weight * b[0];
    r[1] = (1 - weight) * a[1] + weight * b[1];
}

void msdf_linearDirection(msdf_Vec2 r, msdf_EdgeSegment *e, double param)
{
    r[0] = e->p[1][0] - e->p[0][0];
    r[1] = e->p[1][1] - e->p[0][1];
}

void msdf_quadraticDirection(msdf_Vec2 r, msdf_EdgeSegment *e, double param)
{
    msdf_Vec2 a, b;
    msdf_v2Sub(a, e->p[1], e->p[0]);
    msdf_v2Sub(b, e->p[2], e->p[1]);
    msdf_mix(r, a, b, param);
}

void msdf_cubicDirection(msdf_Vec2 r, msdf_EdgeSegment *e, double param)
{
    msdf_Vec2 a, b, c, d, t;
    msdf_v2Sub(a, e->p[1], e->p[0]);
    msdf_v2Sub(b, e->p[2], e->p[1]);
    msdf_mix(c, a, b, param);
    msdf_v2Sub(a, e->p[3], e->p[2]);
    msdf_mix(d, b, a, param);
    msdf_mix(t, c, d, param);

    if (!t[0] && !t[1])
    {
        if (param == 0)
        {
            r[0] = e->p[2][0] - e->p[0][0];
            r[1] = e->p[2][1] - e->p[0][1];
            return;
        }
        if (param == 1)
        {
            r[0] = e->p[3][0] - e->p[1][0];
            r[1] = e->p[3][1] - e->p[1][1];
            return;
        }
    }

    r[0] = t[0];
    r[1] = t[1];
}

void msdf_direction(msdf_Vec2 r, msdf_EdgeSegment *e, double param)
{
    switch (e->type)
    {
    case STBTT_vline:
    {
        msdf_linearDirection(r, e, param);
        break;
    }
    case STBTT_vcurve:
    {
        msdf_quadraticDirection(r, e, param);
        break;
    }
    case STBTT_vcubic:
    {
        msdf_cubicDirection(r, e, param);
        break;
    }
    }
}

void msdf_linearPoint(msdf_Vec2 r, msdf_EdgeSegment *e, double param)
{
    msdf_mix(r, e->p[0], e->p[1], param);
}

void msdf_quadraticPoint(msdf_Vec2 r, msdf_EdgeSegment *e, double param)
{
    msdf_Vec2 a, b;
    msdf_mix(a, e->p[0], e->p[1], param);
    msdf_mix(b, e->p[1], e->p[2], param);
    msdf_mix(r, a, b, param);
}

void msdf_cubicPoint(msdf_Vec2 r, msdf_EdgeSegment *e, double param)
{
    msdf_Vec2 p12, a, b, c, d;
    msdf_mix(p12, e->p[1], e->p[2], param);

    msdf_mix(a, e->p[0], e->p[1], param);
    msdf_mix(b, a, p12, param);

    msdf_mix(c, e->p[2], e->p[3], param);
    msdf_mix(d, p12, c, param);

    msdf_mix(r, b, d, param);
}

void msdf_point(msdf_Vec2 r, msdf_EdgeSegment *e, double param)
{
    switch (e->type)
    {
    case STBTT_vline:
    {
        msdf_linearPoint(r, e, param);
        break;
    }
    case STBTT_vcurve:
    {
        msdf_quadraticPoint(r, e, param);
        break;
    }
    case STBTT_vcubic:
    {
        msdf_cubicPoint(r, e, param);
        break;
    }
    }
}

// linear edge signed distance
msdf_signedDistance msdf_linearDist(msdf_EdgeSegment *e, msdf_Vec2 origin, double *param)
{
    msdf_Vec2 aq, ab, eq;
    msdf_v2Sub(aq, origin, e->p[0]);
    msdf_v2Sub(ab, e->p[1], e->p[0]);
    *param = msdf_v2MulInner(aq, ab) / msdf_v2MulInner(ab, ab);
    msdf_v2Sub(eq, e->p[*param > .5], origin);

    double endpoint_distance = msdf_v2Leng(eq);
    if (*param > 0 && *param < 1)
    {
        msdf_Vec2 ab_ortho;
        msdf_getOrtho(ab_ortho, ab, 0, 0);
        double ortho_dist = msdf_v2MulInner(ab_ortho, aq);
        if (fabs(ortho_dist) < endpoint_distance)
            return (msdf_signedDistance){ortho_dist, 0};
    }

    msdf_v2Norm(ab, ab);
    msdf_v2Norm(eq, eq);
    double dist = msdf_nonZeroSign(msdf_cross(aq, ab)) * endpoint_distance;
    double d = fabs(msdf_v2MulInner(ab, eq));
    return (msdf_signedDistance){dist, d};
}

// quadratic edge signed distance
msdf_signedDistance msdf_quadraticDist(msdf_EdgeSegment *e, msdf_Vec2 origin, double *param)
{
    msdf_Vec2 qa, ab, br;
    msdf_v2Sub(qa, e->p[0], origin);
    msdf_v2Sub(ab, e->p[1], e->p[0]);
    br[0] = e->p[0][0] + e->p[2][0] - e->p[1][0] - e->p[1][0];
    br[1] = e->p[0][1] + e->p[2][1] - e->p[1][1] - e->p[1][1];

    double a = msdf_v2MulInner(br, br);
    double b = 3 * msdf_v2MulInner(ab, br);
    double c = 2 * msdf_v2MulInner(ab, ab) + msdf_v2MulInner(qa, br);
    double d = msdf_v2MulInner(qa, ab);
    double t[3];
    int solutions = msdf_solveCubic(t, a, b, c, d);

    // distance from a
    double min_distance = msdf_nonZeroSign(msdf_cross(ab, qa)) * msdf_v2Leng(qa);
    *param = -msdf_v2MulInner(qa, ab) / msdf_v2MulInner(ab, ab);
    {
        msdf_Vec2 a, b;
        msdf_v2Sub(a, e->p[2], e->p[1]);
        msdf_v2Sub(b, e->p[2], origin);

        // distance from b
        double distance = msdf_nonZeroSign(msdf_cross(a, b)) * msdf_v2Leng(b);
        if (fabs(distance) < fabs(min_distance))
        {
            min_distance = distance;

            msdf_v2Sub(a, origin, e->p[1]);
            msdf_v2Sub(b, e->p[2], e->p[1]);
            *param = msdf_v2MulInner(a, b) / msdf_v2MulInner(b, b);
        }
    }

    for (int i = 0; i < solutions; ++i)
    {
        if (t[i] > 0 && t[i] < 1)
        {
            // end_point = p[0]+2*t[i]*ab+t[i]*t[i]*br;
            msdf_Vec2 end_point, a, b;
            end_point[0] = e->p[0][0] + 2 * t[i] * ab[0] + t[i] * t[i] * br[0];
            end_point[1] = e->p[0][1] + 2 * t[i] * ab[1] + t[i] * t[i] * br[1];

            msdf_v2Sub(a, e->p[2], e->p[0]);
            msdf_v2Sub(b, end_point, origin);
            double distance = msdf_nonZeroSign(msdf_cross(a, b)) * msdf_v2Leng(b);
            if (fabs(distance) <= fabs(min_distance))
            {
                min_distance = distance;
                *param = t[i];
            }
        }
    }

    if (*param >= 0 && *param <= 1)
        return (msdf_signedDistance){min_distance, 0};

    msdf_Vec2 aa, bb;
    msdf_v2Norm(ab, ab);
    msdf_v2Norm(qa, qa);
    msdf_v2Sub(aa, e->p[2], e->p[1]);
    msdf_v2Norm(aa, aa);
    msdf_v2Sub(bb, e->p[2], origin);
    msdf_v2Norm(bb, bb);

    if (*param < .5)
        return (msdf_signedDistance){min_distance, fabs(msdf_v2MulInner(ab, qa))};
    else
        return (msdf_signedDistance){min_distance, fabs(msdf_v2MulInner(aa, bb))};
}

// cubic edge signed distance
msdf_signedDistance msdf_cubicDist(msdf_EdgeSegment *e, msdf_Vec2 origin, double *param)
{
    msdf_Vec2 qa, ab, br, as;
    msdf_v2Sub(qa, e->p[0], origin);
    msdf_v2Sub(ab, e->p[1], e->p[0]);
    br[0] = e->p[2][0] - e->p[1][0] - ab[0];
    br[1] = e->p[2][1] - e->p[1][1] - ab[1];
    as[0] = (e->p[3][0] - e->p[2][0]) - (e->p[2][0] - e->p[1][0]) - br[0];
    as[1] = (e->p[3][1] - e->p[2][1]) - (e->p[2][1] - e->p[1][1]) - br[1];

    msdf_Vec2 ep_dir;
    msdf_direction(ep_dir, e, 0);

    // distance from a
    double min_distance = msdf_nonZeroSign(msdf_cross(ep_dir, qa)) * msdf_v2Leng(qa);
    *param = -msdf_v2MulInner(qa, ep_dir) / msdf_v2MulInner(ep_dir, ep_dir);
    {
        msdf_Vec2 a;
        msdf_v2Sub(a, e->p[3], origin);

        msdf_direction(ep_dir, e, 1);
        // distance from b
        double distance = msdf_nonZeroSign(msdf_cross(ep_dir, a)) * msdf_v2Leng(a);
        if (fabs(distance) < fabs(min_distance))
        {
            min_distance = distance;

            a[0] = origin[0] + ep_dir[0] - e->p[3][0];
            a[1] = origin[1] + ep_dir[1] - e->p[3][1];
            *param = msdf_v2MulInner(a, ep_dir) / msdf_v2MulInner(ep_dir, ep_dir);
        }
    }

    const int search_starts = 4;
    for (int i = 0; i <= search_starts; ++i)
    {
        double t = (double)i / search_starts;
        for (int step = 0;; ++step)
        {
            msdf_Vec2 qpt;
            msdf_point(qpt, e, t);
            msdf_v2Sub(qpt, qpt, origin);
            msdf_Vec2 d;
            msdf_direction(d, e, t);
            double distance = msdf_nonZeroSign(msdf_cross(d, qpt)) * msdf_v2Leng(qpt);
            if (fabs(distance) < fabs(min_distance))
            {
                min_distance = distance;
                *param = t;
            }
            if (step == search_starts)
                break;

            msdf_Vec2 d1, d2;
            d1[0] = 3 * as[0] * t * t + 6 * br[0] * t + 3 * ab[0];
            d1[1] = 3 * as[1] * t * t + 6 * br[1] * t + 3 * ab[1];
            d2[0] = 6 * as[0] * t + 6 * br[0];
            d2[1] = 6 * as[1] * t + 6 * br[1];

            t -= msdf_v2MulInner(qpt, d1) / (msdf_v2MulInner(d1, d1) + msdf_v2MulInner(qpt, d2));
            if (t < 0 || t > 1)
                break;
        }
    }

    if (*param >= 0 && *param <= 1)
        return (msdf_signedDistance){min_distance, 0};

    msdf_Vec2 d0, d1;
    msdf_direction(d0, e, 0);
    msdf_direction(d1, e, 1);
    msdf_v2Norm(d0, d0);
    msdf_v2Norm(d1, d1);
    msdf_v2Norm(qa, qa);
    msdf_Vec2 a;
    msdf_v2Sub(a, e->p[3], origin);
    msdf_v2Norm(a, a);

    if (*param < .5)
        return (msdf_signedDistance){min_distance, fabs(msdf_v2MulInner(d0, qa))};
    else
        return (msdf_signedDistance){min_distance, fabs(msdf_v2MulInner(d1, a))};
}

void msdf_distToPseudo(msdf_signedDistance *distance, msdf_Vec2 origin, double param, msdf_EdgeSegment *e) {
    if (param < 0) {
        msdf_Vec2 dir, p;
        msdf_direction(dir, e, 0);
        msdf_v2Norm(dir, dir);
        msdf_Vec2 aq = {origin[0], origin[1]};
        msdf_point(p, e, 0);
        msdf_v2Sub(aq, origin, p);
        double ts = msdf_v2MulInner(aq, dir);
        if (ts < 0) {
            double pseudo_dist = msdf_cross(aq, dir);
            if (fabs(pseudo_dist) <= fabs(distance->dist)) {
                distance->dist = pseudo_dist;
                distance->d = 0;
            }
        }
    } else if (param > 1) {
        msdf_Vec2 dir, p;
        msdf_direction(dir, e, 1);
        msdf_v2Norm(dir, dir);
        msdf_Vec2 bq = {origin[0], origin[1]};
        msdf_point(p, e, 1);
        msdf_v2Sub(bq, origin, p);
        double ts = msdf_v2MulInner(bq, dir);
        if (ts > 0) {
            double pseudo_dist = msdf_cross(bq, dir);
            if (fabs(pseudo_dist) <= fabs(distance->dist)) {
                distance->dist = pseudo_dist;
                distance->d = 0;
            }
        }
    }
}

int msdf_signedCompare(msdf_signedDistance a, msdf_signedDistance b) {
    return fabs(a.dist) < fabs(b.dist) || (fabs(a.dist) == fabs(b.dist) && a.d < b.d);
}

int msdf_isCorner(msdf_Vec2 a, msdf_Vec2 b, double threshold) {
    return msdf_v2MulInner(a, b) <= 0 || fabs(msdf_cross(a, b)) > threshold;
}

void msdf_switchColor(msdf_edgeColor *color, unsigned long long *seed, msdf_edgeColor banned)
{
    msdf_edgeColor combined = *color & banned;
    if (combined == msdf_edgeColor_red || combined == msdf_edgeColor_green || combined == msdf_edgeColor_blue) {
        *color = (msdf_edgeColor)(combined ^ msdf_edgeColor_white);
        return;
    }

    if (*color == msdf_edgeColor_black || *color == msdf_edgeColor_white) {
        static const msdf_edgeColor start[3] = {msdf_edgeColor_cyan, msdf_edgeColor_magenta, msdf_edgeColor_yellow};
        *color = start[*seed & 3];
        *seed /= 3;
        return;
    }

    int shifted = *color << (1 + (*seed & 1));
    *color = (msdf_edgeColor)((shifted | shifted >> 3) & msdf_edgeColor_white);
    *seed >>= 1;
}

void msdf_linearSplit(msdf_EdgeSegment *e, msdf_EdgeSegment *p1, msdf_EdgeSegment *p2, msdf_EdgeSegment *p3)
{
    msdf_Vec2 p;

    msdf_point(p, e, 1 / 3.0);
    memcpy(&p1->p[0], e->p[0], sizeof(msdf_Vec2));
    memcpy(&p1->p[1], p, sizeof(msdf_Vec2));
    p1->color = e->color;

    msdf_point(p, e, 1 / 3.0);
    memcpy(&p2->p[0], p, sizeof(msdf_Vec2));
    msdf_point(p, e, 2 / 3.0);
    memcpy(&p2->p[1], p, sizeof(msdf_Vec2));
    p2->color = e->color;

    msdf_point(p, e, 2 / 3.0);
    memcpy(&p3->p[0], p, sizeof(msdf_Vec2));
    msdf_point(p, e, 2 / 3.0);
    memcpy(&p3->p[1], e->p[1], sizeof(msdf_Vec2));
    p3->color = e->color;
}

void msdf_quadraticSplit(msdf_EdgeSegment *e, msdf_EdgeSegment *p1, msdf_EdgeSegment *p2, msdf_EdgeSegment *p3)
{
    msdf_Vec2 p, a, b;

    memcpy(&p1->p[0], e->p[0], sizeof(msdf_Vec2));
    msdf_mix(p, e->p[0], e->p[1], 1 / 3.0);
    memcpy(&p1->p[1], p, sizeof(msdf_Vec2));
    msdf_point(p, e, 1 / 3.0);
    memcpy(&p1->p[2], p, sizeof(msdf_Vec2));
    p1->color = e->color;

    msdf_point(p, e, 1 / 3.0);
    memcpy(&p2->p[0], p, sizeof(msdf_Vec2));
    msdf_mix(a, e->p[0], e->p[1], 5 / 9.0);
    msdf_mix(b, e->p[1], e->p[2], 4 / 9.0);
    msdf_mix(p, a, b, 0.5);
    memcpy(&p2->p[1], p, sizeof(msdf_Vec2));
    msdf_point(p, e, 2 / 3.0);
    memcpy(&p2->p[2], p, sizeof(msdf_Vec2));
    p2->color = e->color;

    msdf_point(p, e, 2 / 3.0);
    memcpy(&p3->p[0], p, sizeof(msdf_Vec2));
    msdf_mix(p, e->p[1], e->p[2], 2 / 3.0);
    memcpy(&p3->p[1], p, sizeof(msdf_Vec2));
    memcpy(&p3->p[2], e->p[2], sizeof(msdf_Vec2));
    p3->color = e->color;
}

void msdf_cubicSplit(msdf_EdgeSegment *e, msdf_EdgeSegment *p1, msdf_EdgeSegment *p2, msdf_EdgeSegment *p3)
{
    msdf_Vec2 p, a, b, c, d;

    memcpy(&p1->p[0], e->p[0], sizeof(msdf_Vec2)); // p1 0
    if (e->p[0] == e->p[1]) {
        memcpy(&p1->p[1], e->p[0], sizeof(msdf_Vec2)); // ? p1 1
    } else {
        msdf_mix(p, e->p[0], e->p[1], 1 / 3.0);
        memcpy(&p1->p[1], p, sizeof(msdf_Vec2)); // ? p1 1
    }
    msdf_mix(a, e->p[0], e->p[1], 1 / 3.0);
    msdf_mix(b, e->p[1], e->p[2], 1 / 3.0);
    msdf_mix(p, a, b, 1 / 3.0);
    memcpy(&p1->p[2], p, sizeof(msdf_Vec2)); // p1 2
    msdf_point(p, e, 1 / 3.0);
    memcpy(&p1->p[3], p, sizeof(msdf_Vec2)); // p1 3
    p1->color = e->color;

    msdf_point(p, e, 1 / 3.0);
    memcpy(&p2->p[0], p, sizeof(msdf_Vec2)); // p2 0
    msdf_mix(a, e->p[0], e->p[1], 1 / 3.0);
    msdf_mix(b, e->p[1], e->p[2], 1 / 3.0);
    msdf_mix(c, a, b, 1 / 3.0);
    msdf_mix(a, e->p[1], e->p[2], 1 / 3.0);
    msdf_mix(b, e->p[2], e->p[3], 1 / 3.0);
    msdf_mix(d, a, b, 1 / 3.0);
    msdf_mix(p, c, d, 2 / 3.0);
    memcpy(&p2->p[1], p, sizeof(msdf_Vec2)); // p2 1
    msdf_mix(a, e->p[0], e->p[1], 2 / 3.0);
    msdf_mix(b, e->p[1], e->p[2], 2 / 3.0);
    msdf_mix(c, a, b, 2 / 3.0);
    msdf_mix(a, e->p[1], e->p[2], 2 / 3.0);
    msdf_mix(b, e->p[2], e->p[3], 2 / 3.0);
    msdf_mix(d, a, b, 2 / 3.0);
    msdf_mix(p, c, d, 1 / 3.0);
    memcpy(&p2->p[2], p, sizeof(msdf_Vec2)); // p2 2
    msdf_point(p, e, 2 / 3.0);
    memcpy(&p2->p[3], p, sizeof(msdf_Vec2)); // p2 3
    p2->color = e->color;

    msdf_point(p, e, 2 / 3.0);
    memcpy(&p3->p[0], p, sizeof(msdf_Vec2)); // p3 0

    msdf_mix(a, e->p[1], e->p[2], 2 / 3.0);
    msdf_mix(b, e->p[2], e->p[3], 2 / 3.0);
    msdf_mix(p, a, b, 2 / 3.0);
    memcpy(&p3->p[1], p, sizeof(msdf_Vec2)); // p3 1

    if (e->p[2] == e->p[3]) {
        memcpy(&p3->p[2], e->p[3], sizeof(msdf_Vec2)); // ? p3 2
    } else {
        msdf_mix(p, e->p[2], e->p[3], 2 / 3.0);
        memcpy(&p3->p[2], p, sizeof(msdf_Vec2)); // ? p3 2
    }

    memcpy(&p3->p[3], e->p[3], sizeof(msdf_Vec2)); // p3 3
}

void msdf_edgeSplit(msdf_EdgeSegment *e, msdf_EdgeSegment *p1, msdf_EdgeSegment *p2, msdf_EdgeSegment *p3)
{
    switch (e->type) {
        case STBTT_vline: {
            msdf_linearSplit(e, p1, p2, p3);
            break;
        }
        case STBTT_vcurve: {
            msdf_quadraticSplit(e, p1, p2, p3);
            break;
        }
        case STBTT_vcubic: {
            msdf_cubicSplit(e, p1, p2, p3);
            break;
        }
    }
}

double msdf_shoelace(const msdf_Vec2 a, const msdf_Vec2 b)
{
    return (b[0] - a[0]) * (a[1] + b[1]);
}


void* msdf__alloc(size_t size, void* ctx) {
    return malloc(size);
}
void msdf__free(void* ptr, void* ctx) {
    free(ptr);
}

int msdf_genGlyph(msdf_Result* result, stbtt_fontinfo *font, int stbttGlyphIndex, uint32_t borderWidth, float scale, float range, msdf_AllocCtx* alloc) {
    msdf_AllocCtx allocCtx;

    if (alloc) {
        allocCtx = *alloc;
    } else {
        allocCtx.alloc = msdf__alloc;
        allocCtx.free = msdf__free;
        allocCtx.ctx = NULL;
    }

    //char f = c;
    // Funit to pixel scale
    //float scale = stbtt_ScaleForMappingEmToPixels(font, h);
    int glyphIdx = stbttGlyphIndex;
    // get glyph bounding box (scaled later)
    int ix0, iy0, ix1, iy1;
    float xoff = .0, yoff = .0;
    stbtt_GetGlyphBox(font, glyphIdx, &ix0, &iy0, &ix1, &iy1);

    float glyphWidth = ix1 - ix0;
    float glyphHeight = iy1 - iy0;
    float borderWidthF32 = borderWidth;
    float wF32 = ceilf(glyphWidth  * scale);
    float hF32 = ceilf(glyphHeight * scale);
    wF32 += 2.f * borderWidth;
    hF32 += 2.f * borderWidth;
    int w = wF32;
    int h = hF32;

    float* bitmap = (float*) allocCtx.alloc(w * h * 3 * sizeof(float), allocCtx.ctx);
    memset(bitmap, 0x0, w * h * 3 * sizeof(float));

    // em scale
    //scale = stbtt_ScaleForMappingEmToPixels(font, h);

    //if (autofit)
    //{

        // calculate new height
        //float newh = h + (h - (iy1 - iy0) * scale) - 4;

        // calculate new scale
        // see 'stbtt_ScaleForMappingEmToPixels' in stb_truetype.h
        //uint8_t *p = font->data + font->head + 18;
        //int unitsPerEm = p[0] * 256 + p[1];
        //scale = ((float)h) / ((float)unitsPerEm);

        // make sure we are centered
        //xoff = .0;
        //yoff = .0;
    //}

    // get left offset and advance
    //int left_bearing, advance;
    //stbtt_GetGlyphHMetrics(font, glyphIdx, &advance, &left_bearing);
    //left_bearing *= scale;

    int32_t glyphOrgX = ix0 * scale;
    int32_t glyphOrgY = iy0 * scale;

    int32_t borderWidthX = borderWidth;
    int32_t borderWidthY = borderWidth;

    //   org  8,8
    // - bord 4,4
    // erg:   4,4

    // calculate offset for centering glyph on bitmap
    
    //glyphOrgX >= 2 ? (glyphOrgX) : ();
    int32_t translateX = (glyphOrgX - borderWidth);//borderWidth + ((w / 2) - ((ix1 - ix0) * scale) / 2 - leftBearingScaled);
    int32_t translateY = (glyphOrgY - borderWidth);//borderWidth + ((h / 2) - ((iy1 - iy0) * scale) / 2 - ((float) iy0) * scale);
    //translateY  = 8;
    // set the glyph metrics
    // (pre-scale them)

    #if 0
    if (metrics)
    {
        metrics->left_bearing = left_bearing;
        metrics->advance = advance * scale;
        metrics->ix0 = ix0 * scale;
        metrics->ix1 = ix1 * scale;
        metrics->iy0 = iy0 * scale;
        metrics->iy1 = iy1 * scale;
        metrics->glyphIdx = glyphIdx;
    }
    #endif

    stbtt_vertex *verts;
    int num_verts = stbtt_GetGlyphShape(font, glyphIdx, &verts);

    // figure out how many contours exist
    int contour_count = 0;
    for (int i = 0; i < num_verts; i++) {
        if (verts[i].type == STBTT_vmove) {
            contour_count++;
        }
    }

    if (contour_count == 0) {
        return 0;
    }

    // determin what vertices belong to what contours
    typedef struct {
        size_t start, end;
    } msdf_Indices;
    msdf_Indices *contours = allocCtx.alloc(sizeof(msdf_Indices) * contour_count, allocCtx.ctx);
    int j = 0;
    for (int i = 0; i <= num_verts; i++) {
        if (verts[i].type == STBTT_vmove) {
            if (i > 0) {
                contours[j].end = i;
                j++;
            }

            contours[j].start = i;
        } else if (i >= num_verts) {
            contours[j].end = i;
        }
    }

    typedef struct {
        msdf_signedDistance min_distance;
        msdf_EdgeSegment *near_edge;
        double near_param;
    } msdf_EdgePoint;

    typedef struct {
        msdf_EdgeSegment *edges;
        size_t edge_count;
    } msdf_Contour;

    // process verts into series of contour-specific edge lists
    msdf_Vec2 initial = {0, 0}; // fix this?
    msdf_Contour *contour_data = allocCtx.alloc(sizeof(msdf_Contour) * contour_count, allocCtx.ctx);
    double cscale = 64.0;
    for (int i = 0; i < contour_count; i++) {
        size_t count = contours[i].end - contours[i].start;
        contour_data[i].edges = allocCtx.alloc(sizeof(msdf_EdgeSegment) * count, allocCtx.ctx);
        contour_data[i].edge_count = 0;

        size_t k = 0;
        for (int j = contours[i].start; j < contours[i].end; j++) {
            msdf_EdgeSegment *e = &contour_data[i].edges[k];
            stbtt_vertex *v = &verts[j];
            e->type = v->type;
            e->color = msdf_edgeColor_white;

            switch (v->type) {
                case STBTT_vmove: {
                    msdf_Vec2 p = {v->x / cscale, v->y / cscale};
                    memcpy(&initial, p, sizeof(msdf_Vec2));
                    break;
                }

                case STBTT_vline: {
                    msdf_Vec2 p = {v->x / cscale, v->y / cscale};
                    memcpy(&e->p[0], initial, sizeof(msdf_Vec2));
                    memcpy(&e->p[1], p, sizeof(msdf_Vec2));
                    memcpy(&initial, p, sizeof(msdf_Vec2));
                    contour_data[i].edge_count++;
                    k++;
                    break;
                }

                case STBTT_vcurve: {
                    msdf_Vec2 p = {v->x / cscale, v->y / cscale};
                    msdf_Vec2 c = {v->cx / cscale, v->cy / cscale};
                    memcpy(&e->p[0], initial, sizeof(msdf_Vec2));
                    memcpy(&e->p[1], c, sizeof(msdf_Vec2));
                    memcpy(&e->p[2], p, sizeof(msdf_Vec2));

                    if ((e->p[0][0] == e->p[1][0] && e->p[0][1] == e->p[1][1]) ||
                        (e->p[1][0] == e->p[2][0] && e->p[1][1] == e->p[2][1]))
                    {
                        e->p[1][0] = 0.5 * (e->p[0][0] + e->p[2][0]);
                        e->p[1][1] = 0.5 * (e->p[0][1] + e->p[2][1]);
                    }

                    memcpy(&initial, p, sizeof(msdf_Vec2));
                    contour_data[i].edge_count++;
                    k++;
                    break;
                }

                case STBTT_vcubic: {
                    msdf_Vec2 p = {v->x / cscale, v->y / cscale};
                    msdf_Vec2 c = {v->cx / cscale, v->cy / cscale};
                    msdf_Vec2 c1 = {v->cx1 / cscale, v->cy1 / cscale};
                    memcpy(&e->p[0], initial, sizeof(msdf_Vec2));
                    memcpy(&e->p[1], c, sizeof(msdf_Vec2));
                    memcpy(&e->p[2], c1, sizeof(msdf_Vec2));
                    memcpy(&e->p[3], p, sizeof(msdf_Vec2));
                    memcpy(&initial, p, sizeof(msdf_Vec2));
                    contour_data[i].edge_count++;
                    k++;
                    break;
                }
            }
        }
    }

    // calculate edge-colors
    uint64_t seed = 0;
    double anglethreshold = 3.0;
    double crossthreshold = sin(anglethreshold);
    size_t corner_count = 0;
    for (int i = 0; i < contour_count; ++i) {
        for (int j = 0; j < contour_data[i].edge_count; ++j) {
            corner_count++;
        }
    }

    int *corners = allocCtx.alloc(sizeof(int) * corner_count, allocCtx.ctx);
    int cornerIndex = 0;
    for (int i = 0; i < contour_count; ++i) {

        if (contour_data[i].edge_count > 0) {
            msdf_Vec2 prev_dir, dir;
            msdf_direction(prev_dir, &contour_data[i].edges[contour_data[i].edge_count - 1], 1);

            int index = 0;
            for (int j = 0; j < contour_data[i].edge_count; ++j, ++index) {
                msdf_EdgeSegment *e = &contour_data[i].edges[j];
                msdf_direction(dir, e, 0);
                msdf_v2Norm(dir, dir);
                msdf_v2Norm(prev_dir, prev_dir);
                if (msdf_isCorner(prev_dir, dir, crossthreshold)) {
                    corners[cornerIndex++] = index;
                }
                msdf_direction(prev_dir, e, 1);
            }
        }

        if (cornerIndex == 0) {
            for (int j = 0; j < contour_data[i].edge_count; ++j) {
                contour_data[i].edges[j].color = msdf_edgeColor_white;
            }
        } else if (cornerIndex == 1) {
            msdf_edgeColor colors[3] = {msdf_edgeColor_white, msdf_edgeColor_white};
            msdf_switchColor(&colors[0], &seed, msdf_edgeColor_black);
            colors[2] = colors[0];
            msdf_switchColor(&colors[2], &seed, msdf_edgeColor_black);

            int corner = corners[0];
            if (contour_data[i].edge_count >= 3) {
                int m = contour_data[i].edge_count;
                for (int j = 0; j < m; ++j) {
                    contour_data[i].edges[(corner + j) % m].color = (colors + 1)[(int)(3 + 2.875 * i / (m - 1) - 1.4375 + .5) - 3];
                }
            } else if (contour_data[i].edge_count >= 1) {
                msdf_EdgeSegment *parts[7] = {NULL};
                msdf_edgeSplit(&contour_data[i].edges[0], parts[0 + 3 * corner], parts[1 + 3 * corner], parts[2 + 3 * corner]);
                if (contour_data[i].edge_count >= 2) {
                    msdf_edgeSplit(&contour_data[i].edges[1], parts[3 - 3 * corner], parts[4 - 3 * corner], parts[5 - 3 * corner]);
                    parts[0]->color = parts[1]->color = colors[0];
                    parts[2]->color = parts[3]->color = colors[1];
                    parts[4]->color = parts[5]->color = colors[2];
                } else {
                    parts[0]->color = colors[0];
                    parts[1]->color = colors[1];
                    parts[2]->color = colors[2];
                }
                if (allocCtx.free) {
                    allocCtx.free(contour_data[i].edges, allocCtx.ctx);
                }
                contour_data[i].edges = allocCtx.alloc(sizeof(msdf_EdgeSegment) * 7, allocCtx.ctx);
                contour_data[i].edge_count = 0;
                int index = 0;
                for (int j = 0; parts[j]; ++j) {
                    memcpy(&contour_data[i].edges[index++], &parts[j], sizeof(msdf_EdgeSegment));
                    contour_data[i].edge_count++;
                }
            }
        } else {
            int spline = 0;
            int start = corners[0];
            int m = contour_data[i].edge_count;
            msdf_edgeColor color = msdf_edgeColor_white;
            msdf_switchColor(&color, &seed, msdf_edgeColor_black);
            msdf_edgeColor initial_color = color;
            for (int j = 0; j < m; ++j) {
                int index = (start + j) % m;
                if (spline + 1 < corner_count && corners[spline + 1] == index) {
                    ++spline;

                    msdf_edgeColor s = (msdf_edgeColor)((spline == corner_count - 1) * initial_color);
                    msdf_switchColor(&color, &seed, s);
                }
                contour_data[i].edges[index].color = color;
            }
        }
    }

    if (allocCtx.free) {
        allocCtx.free(corners, allocCtx.ctx);
    }

    // normalize shape
    for (int i = 0; i < contour_count; i++) {
        if (contour_data[i].edge_count == 1) {
            msdf_EdgeSegment *parts[3] = {0};
            msdf_edgeSplit(&contour_data[i].edges[0], parts[0], parts[1], parts[2]);
            if (allocCtx.free) {
                allocCtx.free(contour_data[i].edges, allocCtx.ctx);
            }
            contour_data[i].edges = allocCtx.alloc(sizeof(msdf_EdgeSegment) * 3, allocCtx.ctx);
            contour_data[i].edge_count = 3;
            for (int j = 0; j < 3; j++) {
                memcpy(&contour_data[i].edges[j], &parts[j], sizeof(msdf_EdgeSegment));
            }
        }
    }

    // calculate windings
    int *windings = allocCtx.alloc(sizeof(int) * contour_count, allocCtx.ctx);
    for (int i = 0; i < contour_count; i++) {
        size_t edge_count = contour_data[i].edge_count;
        if (edge_count == 0) {
            windings[i] = 0;
            continue;
        }

        double total = 0;

        if (edge_count == 1) {
            msdf_Vec2 a, b, c;
            msdf_point(a, &contour_data[i].edges[0], 0);
            msdf_point(b, &contour_data[i].edges[0], 1 / 3.0);
            msdf_point(c, &contour_data[i].edges[0], 2 / 3.0);
            total += msdf_shoelace(a, b);
            total += msdf_shoelace(b, c);
            total += msdf_shoelace(c, a);
        } else if (edge_count == 2) {
            msdf_Vec2 a, b, c, d;
            msdf_point(a, &contour_data[i].edges[0], 0);
            msdf_point(b, &contour_data[i].edges[0], 0.5);
            msdf_point(c, &contour_data[i].edges[1], 0);
            msdf_point(d, &contour_data[i].edges[1], 0.5);
            total += msdf_shoelace(a, b);
            total += msdf_shoelace(b, c);
            total += msdf_shoelace(c, d);
            total += msdf_shoelace(d, a);
        } else {
            msdf_Vec2 prev;
            msdf_point(prev, &contour_data[i].edges[edge_count - 1], 0);
            for (int j = 0; j < edge_count; j++) {
                msdf_Vec2 cur;
                msdf_point(cur, &contour_data[i].edges[j], 0);
                total += msdf_shoelace(prev, cur);
                memcpy(prev, cur, sizeof(msdf_Vec2));
            }
        }

        windings[i] = ((0 < total) - (total < 0)); // sign
    }

    typedef struct {
        double r, g, b;
        double med;
    } msdf_MultiDistance;

    msdf_MultiDistance *contour_sd;
    contour_sd = allocCtx.alloc(sizeof(msdf_MultiDistance) * contour_count, allocCtx.ctx);

    float invRange = 1.0 / range;

    for (int y = 0; y < h; ++y) {
        int row = iy0 > iy1 ? y : h - y - 1;
        for (int x = 0; x < w; ++x) {
            float a64 = 64.0;
            msdf_Vec2 p = {(translateX + x + xoff) / (scale * a64), (translateY + y + yoff) / (scale * a64)};
            //p[0] = ;
            //p[1] = ;
            msdf_EdgePoint sr, sg, sb;
            sr.near_edge = sg.near_edge = sb.near_edge = NULL;
            sr.near_param = sg.near_param = sb.near_param = 0;
            sr.min_distance.dist = sg.min_distance.dist = sb.min_distance.dist = MSDF_INF;
            sr.min_distance.d = sg.min_distance.d = sb.min_distance.d = 1;
            double d = fabs(MSDF_INF);
            double neg_dist = -MSDF_INF;
            double pos_dist = MSDF_INF;
            int winding = 0;

            // calculate distance to contours from current point (and if its inside or outside of the shape?)
            for (int j = 0; j < contour_count; ++j) {
                msdf_EdgePoint r, g, b;
                r.near_edge = g.near_edge = b.near_edge = NULL;
                r.near_param = g.near_param = b.near_param = 0;
                r.min_distance.dist = g.min_distance.dist = b.min_distance.dist = MSDF_INF;
                r.min_distance.d = g.min_distance.d = b.min_distance.d = 1;

                for (int k = 0; k < contour_data[j].edge_count; ++k) {
                    msdf_EdgeSegment *e = &contour_data[j].edges[k];
                    double param;
                    msdf_signedDistance distance;
                    distance.dist = MSDF_INF;
                    distance.d = 1;

                    // calculate signed distance
                    switch (e->type) {
                    case STBTT_vline: {
                        distance = msdf_linearDist(e, p, &param);
                        break;
                    }
                    case STBTT_vcurve: {
                        distance = msdf_quadraticDist(e, p, &param);
                        break;
                    }
                    case STBTT_vcubic: {
                        distance = msdf_cubicDist(e, p, &param);
                        break;
                    }
                    }

                    if (e->color & msdf_edgeColor_red && msdf_signedCompare(distance, r.min_distance)) {
                        r.min_distance = distance;
                        r.near_edge = e;
                        r.near_param = param;
                    }
                    if (e->color & msdf_edgeColor_green && msdf_signedCompare(distance, g.min_distance)) {
                        g.min_distance = distance;
                        g.near_edge = e;
                        g.near_param = param;
                    }
                    if (e->color & msdf_edgeColor_blue && msdf_signedCompare(distance, b.min_distance)) {
                        b.min_distance = distance;
                        b.near_edge = e;
                        b.near_param = param;
                    }
                }

                if (msdf_signedCompare(r.min_distance, sr.min_distance)) {
                    sr = r;
                }
                if (msdf_signedCompare(g.min_distance, sg.min_distance)) {
                    sg = g;
                }
                if (msdf_signedCompare(b.min_distance, sb.min_distance)) {
                    sb = b;
                }

                double med_min_dist = fabs(msdf_median(r.min_distance.dist, g.min_distance.dist, b.min_distance.dist));

                if (med_min_dist < d) {
                    d = med_min_dist;
                    winding = -windings[j];
                }

                if (r.near_edge) {
                    msdf_distToPseudo(&r.min_distance, p, r.near_param, r.near_edge);
                }
                if (g.near_edge) {
                    msdf_distToPseudo(&g.min_distance, p, g.near_param, g.near_edge);
                }
                if (b.near_edge) {
                    msdf_distToPseudo(&b.min_distance, p, b.near_param, b.near_edge);
                }

                med_min_dist = msdf_median(r.min_distance.dist, g.min_distance.dist, b.min_distance.dist);
                contour_sd[j].r = r.min_distance.dist;
                contour_sd[j].g = g.min_distance.dist;
                contour_sd[j].b = b.min_distance.dist;
                contour_sd[j].med = med_min_dist;

                if (windings[j] > 0 && med_min_dist >= 0 && fabs(med_min_dist) < fabs(pos_dist)) {
                    pos_dist = med_min_dist;
                }
                if (windings[j] < 0 && med_min_dist <= 0 && fabs(med_min_dist) < fabs(neg_dist)) {
                    neg_dist = med_min_dist;
                }
            }

            if (sr.near_edge) {
                msdf_distToPseudo(&sr.min_distance, p, sr.near_param, sr.near_edge);
            }
            if (sg.near_edge) {
                msdf_distToPseudo(&sg.min_distance, p, sg.near_param, sg.near_edge);
            }
            if (sb.near_edge) {
                msdf_distToPseudo(&sb.min_distance, p, sb.near_param, sb.near_edge);
            }

            msdf_MultiDistance msd;
            msd.r = msd.g = msd.b = msd.med = MSDF_INF;
            if (pos_dist >= 0 && fabs(pos_dist) <= fabs(neg_dist)) {
                msd.med = MSDF_INF;
                winding = 1;
                for (int i = 0; i < contour_count; ++i) {
                    if (windings[i] > 0 && contour_sd[i].med > msd.med && fabs(contour_sd[i].med) < fabs(neg_dist)) {
                        msd = contour_sd[i];
                    }
                }
            } else if (neg_dist <= 0 && fabs(neg_dist) <= fabs(pos_dist)) {
                msd.med = -MSDF_INF;
                winding = -1;
                for (int i = 0; i < contour_count; ++i) {
                    if (windings[i] < 0 && contour_sd[i].med < msd.med && fabs(contour_sd[i].med) < fabs(pos_dist)) {
                        msd = contour_sd[i];
                    }
                }
            }

            for (int i = 0; i < contour_count; ++i) {
                if (windings[i] != winding && fabs(contour_sd[i].med) < fabs(msd.med)) {
                    msd = contour_sd[i];
                }
            }

            if (msdf_median(sr.min_distance.dist, sg.min_distance.dist, sb.min_distance.dist) == msd.med) {
                msd.r = sr.min_distance.dist;
                msd.g = sg.min_distance.dist;
                msd.b = sb.min_distance.dist;
            }

            size_t index = 3 * ((row * w) + x);

            float mr = ((float)msd.r) * invRange + 0.5f;
            float mg = ((float)msd.g) * invRange + 0.5f;
            float mb = ((float)msd.b) * invRange + 0.5f;
            bitmap[index + 0] = mr;
            bitmap[index + 1] = mg;
            bitmap[index + 2] = mb;
            
        }
    }

    if (allocCtx.free) {
        for (int i = 0; i < contour_count; i++) {
            allocCtx.free(contour_data[i].edges, allocCtx.ctx);
        }
        allocCtx.free(contour_data, allocCtx.ctx);
        allocCtx.free(contour_sd, allocCtx.ctx);
        allocCtx.free(contours, allocCtx.ctx);
        allocCtx.free(windings, allocCtx.ctx);
        allocCtx.free(verts, allocCtx.ctx);
    }

    // msdf error correction
    typedef struct {
        int x, y;
    } msdf_Clash;
    msdf_Clash *clashes = allocCtx.alloc(sizeof(msdf_Clash) * w * h, allocCtx.ctx);
    size_t cindex = 0;

    double tx = MSDF_EDGE_THRESHOLD / (scale * range);
    double ty = MSDF_EDGE_THRESHOLD / (scale * range);
    for (int y = 0; y < h; y++) {
        for (int x = 0; x < w; x++) {
            if ((x > 0 && msdf_pixelClash(msdf_pixelAt(x, y, w, bitmap), msdf_pixelAt(msdf_max(x - 1, 0), y, w, bitmap), tx)) || (x < w - 1 && msdf_pixelClash(msdf_pixelAt(x, y, w, bitmap), msdf_pixelAt(msdf_min(x + 1, w - 1), y, w, bitmap), tx)) || (y > 0 && msdf_pixelClash(msdf_pixelAt(x, y, w, bitmap), msdf_pixelAt(x, msdf_max(y - 1, 0), w, bitmap), ty)) || (y < h - 1 && msdf_pixelClash(msdf_pixelAt(x, y, w, bitmap), msdf_pixelAt(x, msdf_min(y + 1, h - 1), w, bitmap), ty))) {
                clashes[cindex].x = x;
                clashes[cindex++].y = y;
            }
        }
    }

    for (int i = 0; i < cindex; i++) {
        size_t index = 3 * ((clashes[i].y * w) + clashes[i].x);
        float med = msdf_median(bitmap[index], bitmap[index + 1], bitmap[index + 2]);
        bitmap[index + 0] = med;
        bitmap[index + 1] = med;
        bitmap[index + 2] = med;
    }

    if (allocCtx.free) {
        allocCtx.free(clashes, allocCtx.ctx);
    }

    result->glyphIdx = glyphIdx;
    result->rgb = bitmap;
    result->width = w;
    result->height = h;
    result->yOffset = translateY;

    return 1;
}
#endif
#endif // MSDF_H
