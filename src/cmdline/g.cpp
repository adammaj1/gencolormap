/*


for all versions
sudo apt-get install qt6-base-dev
cd gencolormap
mkdir build
cd build
cmake ..
make


*/
/*
 * Copyright (C) 2015, 2016, 2017, 2018, 2019, 2020, 2021, 2022
 * Computer Graphics Group, University of Siegen
 * Written by Martin Lambers <martin.lambers@uni-siegen.de>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */


/* Notes about the color spaces used internally:
 *
 * - We use D65 white everywhere
 * - RGB means linear RGB; we also have sRGB
 * - RGB and sRGB values are in [0,1]
 * - XYZ, LUV, and similar values are in the original range (not normalized);
 *   often this is [0,100]
 * - All angles (for hue) are measured in radians
 */
 
 
/* Generate color maps for scientific visualization purposes.
 *
 * Usage:
 * - Decide which type of color map you need and how many colors the map should
 *   contain.
 * - Allocate memory for you color map (3 * unsigned char for each color entry).
 * - Call the function that generates your color map.
 *   The return value is always the number of colors that had to be clipped
 *   to fit into sRGB; you want to keep that number low by adjusting parameters.
 *
 * All colors are represented as unsigned char sRGB triplets, with each value in
 * [0,255].
 */




#include <vector>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <cstring>
#include <cmath>
#include <getopt.h> // https://www.gnu.org/software/libc/manual/html_node/Getopt.html   Parsing Long Options with getopt_long




// To communicate extra information back to the program, a few global extern variables are referenced by the program to fetch information from getopt:
extern char *optarg;
extern int optind;


//#include "colormap.hpp"
//#include "export.hpp"

enum type {
    brewer_seq,
    brewer_div,
    brewer_qual,
    puseq_lightness,
    puseq_saturation,
    puseq_rainbow,
    puseq_blackbody,
    puseq_multihue,
    pudiv_lightness,
    pudiv_saturation,
    puqual_hue,
    cubehelix,
    moreland,
    mcnames
};

enum format {
    csv,
    json,
    ppm
};


/*
 * Brewer-like color maps, as described in
 * M. Wijffelaars, R. Vliegen, J.J. van Wijk, E.-J. van der Linden. Generating
 * color palettes using intuitive parameters. In Computer Graphics Forum, vol. 27,
 * no. 3, pp. 743-750, 2008.
 */

// Create a sequential color map with n colors of the given hue in [0,2*PI].

const float BrewerSequentialDefaultHue = 4.18879020479f; // 240 deg
const float BrewerSequentialDefaultContrast = 0.88f;
//float BrewerSequentialDefaultContrastForSmallN(int n); // only for discrete color maps, i.e. n <= 9
const float BrewerSequentialDefaultSaturation = 0.6f;
const float BrewerSequentialDefaultBrightness = 0.75f;
const float BrewerSequentialDefaultWarmth = 0.15f;



// Create a diverging color map with n colors. Half of them will have the given
// hue (in [0,2*PI]), the other half will have a hue that has the distance given
// by divergence (in [0,2*PI]) to that hue, and they will meet in the middle at
// a neutral color.

const float BrewerDivergingDefaultHue = 4.18879020479f; // 240 deg
const float BrewerDivergingDefaultDivergence = 4.18879020479f; // 240 deg = 2/3 * 2PI
const float BrewerDivergingDefaultContrast = 0.88f;
//float BrewerDivergingDefaultContrastForSmallN(int n); // only for discrete color maps, i.e. n <= 9
const float BrewerDivergingDefaultSaturation = 0.6f;
const float BrewerDivergingDefaultBrightness = 0.75f;
const float BrewerDivergingDefaultWarmth = 0.15f;


// Create a qualitative color map with n colors. The colors will have the same
// saturation; lightness and hue will differ. The parameter hue sets the hue of
// the first color, and the parameter divergence defines the hue range starting
// from that hue that can be used for the colors.

const float BrewerQualitativeDefaultHue = 0.0f;
const float BrewerQualitativeDefaultDivergence = 4.18879020479f; // 2/3 * 2PI
const float BrewerQualitativeDefaultContrast = 0.5f;
const float BrewerQualitativeDefaultSaturation = 0.5f;
const float BrewerQualitativeDefaultBrightness = 0.8f;


/*
 * Perceptually unifrom (PU) color maps.
 *
 * These are computed in CIELUV/LCH to achieve approximate perceptual uniformity.
 * One or two of lightness, saturation, and hue are changed while the other(s)
 * stay constant.
 *
 * For example, color maps of constant lightness are useful for 3D surfaces
 * because they do not interfere with additional shading.
 *
 * Relevant paper:
 * M. Lambers. Interactive Creation of Perceptually Uniform Color Maps.
 * Proc. EuroVis Short Papers, May 2020. https://dx.doi.org/10.2312/evs.20201048
 */

/* Sequential perceptually uniform maps */

// Varying lightness

const float PUSequentialLightnessDefaultLightnessRange = 0.95f;
const float PUSequentialLightnessDefaultSaturationRange = 0.95f;
const float PUSequentialLightnessDefaultSaturation = 0.42f;
const float PUSequentialLightnessDefaultHue = 0.349065850399f; // 20 deg


// Varying saturation

const float PUSequentialSaturationDefaultSaturationRange = PUSequentialLightnessDefaultSaturationRange;
const float PUSequentialSaturationDefaultLightness = 0.5f;
const float PUSequentialSaturationDefaultSaturation = PUSequentialLightnessDefaultSaturation;
const float PUSequentialSaturationDefaultHue = 0.349065850399f; // 20 deg



// Varying hue (through all colors, rainbow-like)

const float PUSequentialRainbowDefaultLightnessRange = PUSequentialLightnessDefaultLightnessRange;
const float PUSequentialRainbowDefaultSaturationRange = PUSequentialLightnessDefaultSaturationRange;
const float PUSequentialRainbowDefaultHue = 0.0f;
const float PUSequentialRainbowDefaultRotations = -1.5f;
const float PUSequentialRainbowDefaultSaturation = 0.8f;



// Varying hue (through physically plausible black body colors at increasing
// temperatures).
// The defaults are chosen so that we start at red and arrive at the D65 white
// point (6500 K), thus excluding the blue colors that occur at higher
// temperatures.
const float PUSequentialBlackBodyDefaultTemperature = 250.0f;
const float PUSequentialBlackBodyDefaultTemperatureRange = 6250.0f;
const float PUSequentialBlackBodyDefaultLightnessRange = PUSequentialLightnessDefaultLightnessRange;
const float PUSequentialBlackBodyDefaultSaturationRange = PUSequentialLightnessDefaultSaturationRange;
const float PUSequentialBlackBodyDefaultSaturation = 1.4f;



// Varying hue (user definable)
const float PUSequentialMultiHueDefaultLightnessRange = PUSequentialLightnessDefaultLightnessRange;
const float PUSequentialMultiHueDefaultSaturationRange = PUSequentialSaturationDefaultSaturationRange;
const float PUSequentialMultiHueDefaultSaturation = 0.38f;
const int PUSequentialMultiHueDefaultHues = 2; // number of hues defined in the following lists
const float PUSequentialMultiHueDefaultHueValues[] = { 0.0f, 1.0471975512f }; // hues values in radians in [0,2pi]
const float PUSequentialMultiHueDefaultHuePositions[] = { 0.25f, 0.75f }; // hue positions in [0,1] sorted in ascending order



/* Diverging perceptually uniform maps */

// Varying lightness

const float PUDivergingLightnessDefaultLightnessRange = PUSequentialLightnessDefaultLightnessRange;
const float PUDivergingLightnessDefaultSaturationRange = PUSequentialLightnessDefaultSaturationRange;
const float PUDivergingLightnessDefaultSaturation = PUSequentialLightnessDefaultSaturation;
const float PUDivergingLightnessDefaultHue = 0.349065850399f; // 20 deg
const float PUDivergingLightnessDefaultDivergence = 4.18879020479f; // 2/3 * 2PI



// Varying saturation

const float PUDivergingSaturationDefaultSaturationRange = PUSequentialSaturationDefaultSaturationRange;
const float PUDivergingSaturationDefaultLightness = 0.5f;
const float PUDivergingSaturationDefaultSaturation = 0.45f;
const float PUDivergingSaturationDefaultHue = 0.349065850399f; // 20 deg
const float PUDivergingSaturationDefaultDivergence = 4.18879020479f; // 2/3 * 2PI


/* Qualitative perceptually uniform maps */

const float PUQualitativeHueDefaultHue = 0.0f;
const float PUQualitativeHueDefaultDivergence = 4.18879020479f; // 2/3 * 2PI
const float PUQualitativeHueDefaultLightness = 0.55f;
const float PUQualitativeHueDefaultSaturation = 0.15f;



/*
 * CubeHelix color maps, as described in
 * Green, D. A., 2011, A colour scheme for the display of astronomical intensity
 * images, Bulletin of the Astronomical Society of India, 39, 289.
 */

// Create a CubeHelix colormap with n colors. The parameter hue (in [0,2*PI])
// sets the hue of the first color. The parameter rot sets the number of
// rotations. It can be negative for backwards rotation. The saturation parameter
// determines the saturation of the colors; higher values may lead to clipping
// of colors in the sRGB space. The gamma parameter sets optional gamma correction.
// The return value is the number of colors that had to be clipped.

const float CubeHelixDefaultHue = 0.523598775598f; // 1/12 * 2PI
const float CubeHelixDefaultRotations = -1.5f;
const float CubeHelixDefaultSaturation = 1.2f;
const float CubeHelixDefaultGamma = 1.0f;


/*
 * Moreland color maps, as described in
 * K. Moreland, Diverging Color Maps for Scientific Visualization, Proc. Int.
 * Symp. Visual Computing, December 2009, DOI 10.1007/978-3-642-10520-3_9.
 */

// Create a Moreland colormap with n colors. Specify the two endpoints
// of the colormap as sRGB colors; all intermediate colors will be generated.

const unsigned char MorelandDefaultR0 = 180;
const unsigned char MorelandDefaultG0 = 4;
const unsigned char MorelandDefaultB0 = 38;
const unsigned char MorelandDefaultR1 = 59;
const unsigned char MorelandDefaultG1 = 76;
const unsigned char MorelandDefaultB1 = 192;


/*
 * McNames color maps, as described in
 * J. McNames, An Effective Color Scale for Simultaneous Color and Gray-Scale Publications,
 * IEEE Signal Processing Magazine 23(1), January 2006, DOI 10.1109/MSP.2006.1593340.
 *
 * Note: Use CubeHelix instead! The McNames color maps are perceptually not linear in luminance!
 */

// Create a McNames colormap with n colors. Specify the number of
// periods.

const float McNamesDefaultPeriods = 2.0f;



/* Generic helpers */

static const float pi = M_PI;
static const float twopi = 2.0 * M_PI;

static float sqr(float x)
{
    return x * x;
}

static float clamp(float x, float lo, float hi)
{
    return std::min(std::max(x, lo), hi);
}

static float hue_diff(float h0, float h1)
{
    float t = std::fabs(h1 - h0);
    return (t < pi ? t : twopi - t);
}

static float uchar_to_float(unsigned char x)
{
    return x / 255.0f;
}

static unsigned char float_to_uchar(float x, bool* clipped = NULL)
{
    if (clipped) {
        *clipped = (x < 0.0f || x > 1.0f);
    }
    int v = std::round(x * 255.0f);
    return v < 0 ? 0 : v > 255 ? 255 : v;
}

/* A color triplet class without assumptions about the color space */

class triplet {
public:
    struct {
        union { float x, l,       m, r; };
        union { float y, a, u, c, s, g; };
        union { float z, b, v, h      ; };
    };

    triplet() {}
    triplet(float t0, float t1, float t2) : x(t0), y(t1), z(t2) {}
    triplet(const triplet& t) : x(t.x), y(t.y), z(t.z) {}
    triplet& operator=(const triplet& t) { x = t.x; y = t.y; z = t.z; return *this; }
};

triplet operator+(triplet t0, triplet t1)
{
    return triplet(t0.x + t1.x, t0.y + t1.y, t0.z + t1.z);
}

triplet operator*(float s, triplet t)
{
    return triplet(s * t.x, s * t.y, s * t.z);
}

/* XYZ and related color spaces helper functions and values */

static float u_prime(triplet xyz)
{
    return 4.0f * xyz.x / (xyz.x + 15.0f * xyz.y + 3.0f * xyz.z);
}

static float v_prime(triplet xyz)
{
    return 9.0f * xyz.y / (xyz.x + 15.0f * xyz.y + 3.0f * xyz.z);
}

static const triplet d65_xyz = triplet(95.047f, 100.000f, 108.883f);
static const float d65_u_prime = u_prime(d65_xyz);
static const float d65_v_prime = v_prime(d65_xyz);

static triplet adjust_y(triplet xyz, float new_y)
{
    float sum = xyz.x + xyz.y + xyz.z;
    // keep old chromaticity in terms of x, y
    float x = xyz.x / sum;
    float y = xyz.y / sum;
    // apply new luminance
    float r = new_y / y;
    return triplet(r * x, new_y, r * (1.0f - x - y));
}

/* Color space conversion: LCH <-> LUV */

static triplet lch_to_luv(triplet lch)
{
    return triplet(lch.l, lch.c * std::cos(lch.h), lch.c * std::sin(lch.h));
}

static triplet luv_to_lch(triplet luv)
{
    triplet lch(luv.l, std::hypot(luv.u, luv.v), std::atan2(luv.v, luv.u));
    if (lch.h < 0.0f)
        lch.h += twopi;
    return lch;
}

static float lch_saturation(float l, float c)
{
    return c / std::max(l, 1e-8f);
}

static float lch_chroma(float l, float s)
{
    return s * l;
}

static float lch_distance(triplet lch0, triplet lch1)
{
    /* We have to compute the euclidean distance in LUV space so that it is
     * perceptually uniform. But using the equations above we can simplify
     * the resulting expression to this: */
    return std::sqrt(sqr(lch0.l - lch1.l) + sqr(lch0.c) + sqr(lch1.c)
            - 2.0f * lch0.c * lch1.c * std::cos(lch0.h - lch1.h));
}

/* Color space conversion: LUV <-> XYZ */

static triplet luv_to_xyz(triplet luv)
{
    if (luv.l <= 0.0f)
        return triplet(0.0f, 0.0f, 0.0f);

    triplet xyz;
    float u_prime = luv.u / (13.0f * luv.l) + d65_u_prime;
    float v_prime = luv.v / (13.0f * luv.l) + d65_v_prime;
    if (luv.l <= 8.0f) {
        xyz.y = d65_xyz.y * luv.l * (3.0f * 3.0f * 3.0f / (29.0f * 29.0f * 29.0f));
    } else {
        float tmp = (luv.l + 16.0f) / 116.0f;
        xyz.y = d65_xyz.y * tmp * tmp * tmp;
    }
    xyz.x = xyz.y * (9.0f * u_prime) / (4.0f * v_prime);
    xyz.z = xyz.y * (12.0f - 3.0f * u_prime - 20.0f * v_prime) / (4.0f * v_prime);
    return xyz;
}

static triplet xyz_to_luv(triplet xyz)
{
    triplet luv;
    float y_ratio = xyz.y / d65_xyz.y;
    if (y_ratio <= (6.0f * 6.0f * 6.0f) / (29.0f * 29.0f * 29.0f)) {
        luv.l = (29.0f * 29.0f * 29.0f) / (3.0f * 3.0f * 3.0f) * y_ratio;
    } else {
        luv.l = 116.0f * std::cbrt(y_ratio) - 16.0f;
    }
    luv.u = 13.0f * luv.l * (u_prime(xyz) - d65_u_prime);
    luv.v = 13.0f * luv.l * (v_prime(xyz) - d65_v_prime);
    return luv;
}

static float luv_saturation(triplet luv)
{
    return lch_saturation(luv.l, std::hypot(luv.u, luv.v));
}

/* Color space conversion: LAB <-> XYZ */

static float lab_invf(float t)
{
    if (t > 6.0f / 29.0f)
        return t * t * t;
    else
        return (3.0f * 6.0f * 6.0f) / (29.0f * 29.0f) * (t - 4.0f / 29.0f);
}

static triplet lab_to_xyz(triplet lab)
{
    float t = (lab.l + 16.0f) / 116.0f;
    return triplet(d65_xyz.x * lab_invf(t + lab.a / 500.0f),
                   d65_xyz.y * lab_invf(t),
                   d65_xyz.z * lab_invf(t - lab.b / 200.0f));
}

static float lab_f(float t)
{
    if (t > (6.0f * 6.0f * 6.0f) / (29.0f * 29.0f * 29.0f))
        return std::cbrt(t);
    else
        return (29.0f * 29.0f) / (3.0f * 6.0f * 6.0f) * t + 4.0f / 29.0f;
}

static triplet xyz_to_lab(triplet xyz)
{
    triplet fxyz = triplet(lab_f(xyz.x / d65_xyz.x), lab_f(xyz.y / d65_xyz.y), lab_f(xyz.z / d65_xyz.z));
    return triplet(116.0f * fxyz.y - 16.0f, 500.0f * (fxyz.x - fxyz.y), 200.0f * (fxyz.y - fxyz.z));
}

/* Color space conversion: RGB <-> XYZ
 * The matrix entries are the same as used by PBRT3 and Mitsuba2.
 * Other values exist, some claiming to be more precise, but we stick
 * to these standard values for comparability. */

static triplet rgb_to_xyz(triplet rgb)
{
    return 100.0f * triplet(
            (0.412453f * rgb.r + 0.357580f * rgb.g + 0.180423f * rgb.b),
            (0.212671f * rgb.r + 0.715160f * rgb.g + 0.072169f * rgb.b),
            (0.019334f * rgb.r + 0.119193f * rgb.g + 0.950227f * rgb.b));
}

static triplet xyz_to_rgb(triplet xyz)
{
    return 0.01f * triplet(
            (+3.240479f * xyz.x - 1.537150f * xyz.y - 0.498535f * xyz.z),
            (-0.969256f * xyz.x + 1.875991f * xyz.y + 0.041556f * xyz.z),
            (+0.055648f * xyz.x - 0.204023f * xyz.y + 1.057311f * xyz.z));
}

/* Color space conversion: RGB <-> sRGB */

static float rgb_to_srgb_helper(float x)
{
    return (x <= 0.0031308f ? (x * 12.92f) : (1.055f * std::pow(x, 1.0f / 2.4f) - 0.055f));
}

static triplet rgb_to_srgb(triplet rgb)
{
    return triplet(rgb_to_srgb_helper(rgb.r), rgb_to_srgb_helper(rgb.g), rgb_to_srgb_helper(rgb.b));
}

static float srgb_to_rgb_helper(float x)
{
    return (x <= 0.04045f ? (x / 12.92f) : std::pow((x + 0.055f) / 1.055f, 2.4f));
}

static triplet srgb_to_rgb(triplet srgb)
{
    return triplet(srgb_to_rgb_helper(srgb.r), srgb_to_rgb_helper(srgb.g), srgb_to_rgb_helper(srgb.b));
}

/* Helpers for the conversion to colormap entries */

static bool xyz_to_colormap(triplet xyz, unsigned char* colormap)
{
    bool clipped[3];
    triplet srgb = rgb_to_srgb(xyz_to_rgb(xyz));
    colormap[0] = float_to_uchar(srgb.r, clipped + 0);
    colormap[1] = float_to_uchar(srgb.g, clipped + 1);
    colormap[2] = float_to_uchar(srgb.b, clipped + 2);
    return clipped[0] || clipped[1] || clipped[2];
}

static bool luv_to_colormap(triplet luv, unsigned char* colormap)
{
    return xyz_to_colormap(luv_to_xyz(luv), colormap);
}

static bool lch_to_colormap(triplet lch, unsigned char* colormap)
{
    return luv_to_colormap(lch_to_luv(lch), colormap);
}

static bool lab_to_colormap(triplet lab, unsigned char* colormap)
{
    return xyz_to_colormap(lab_to_xyz(lab), colormap);
}

/* Various helpers */

static float srgb_to_lch_hue(triplet srgb)
{
    return luv_to_lch(xyz_to_luv(rgb_to_xyz(srgb_to_rgb(srgb)))).h;
}

// Compute most saturated color that fits into the sRGB
// cube for the given LCH hue value. This is the core
// of the Wijffelaars paper.
static triplet most_saturated_in_srgb(float lch_hue)
{
    /* Static values, only computed once */
    static float h[] = {
        srgb_to_lch_hue(triplet(1, 0, 0)),
        srgb_to_lch_hue(triplet(1, 1, 0)),
        srgb_to_lch_hue(triplet(0, 1, 0)),
        srgb_to_lch_hue(triplet(0, 1, 1)),
        srgb_to_lch_hue(triplet(0, 0, 1)),
        srgb_to_lch_hue(triplet(1, 0, 1))
    };

    /* RGB values and variable pointers to them */
    int i, j, k;
    if (lch_hue < h[0]) {
        i = 2;
        j = 1;
        k = 0;
    } else if (lch_hue < h[1]) {
        i = 1;
        j = 2;
        k = 0;
    } else if (lch_hue < h[2]) {
        i = 0;
        j = 2;
        k = 1;
    } else if (lch_hue < h[3]) {
        i = 2;
        j = 0;
        k = 1;
    } else if (lch_hue < h[4]) {
        i = 1;
        j = 0;
        k = 2;
    } else if (lch_hue < h[5]) {
        i = 0;
        j = 1;
        k = 2;
    } else {
        i = 2;
        j = 1;
        k = 0;
    }

    /* Compute the missing component */
    float srgb[3];
    float M[3][3] = {
        { 0.4124f, 0.3576f, 0.1805f },
        { 0.2126f, 0.7152f, 0.0722f },
        { 0.0193f, 0.1192f, 0.9505f }
    };
    float alpha = -std::sin(lch_hue);
    float beta = std::cos(lch_hue);
    float T = alpha * d65_u_prime + beta * d65_v_prime;
    srgb[j] = 0.0f;
    srgb[k] = 1.0f;
    float q0 = T * (M[0][k] + 15.0f * M[1][k] + 3.0f * M[2][k]) - (4.0f * alpha * M[0][k] + 9.0f * beta * M[1][k]);
    float q1 = T * (M[0][i] + 15.0f * M[1][i] + 3.0f * M[2][i]) - (4.0f * alpha * M[0][i] + 9.0f * beta * M[1][i]);
    srgb[i] = rgb_to_srgb_helper(clamp(- q0 / q1, 0.0f, 1.0f));

    /* Convert back to LUV */
    return xyz_to_luv(rgb_to_xyz(srgb_to_rgb(triplet(srgb[0], srgb[1], srgb[2]))));
}

static float s_max(float l, float h)
{
    triplet pmid = most_saturated_in_srgb(h);
    triplet pend = triplet(0.0f, 0.0f, 0.0f);
    if (l > pmid.l)
        pend.l = 100.0f;
    float alpha = (pend.l - l) / (pend.l - pmid.l);
    float pmids = luv_saturation(pmid);
    float pends = luv_saturation(pend);
    return alpha * (pmids - pends) + pends;
}

static triplet get_bright_point()
{
    static triplet pb(-1.0f, -1.0f, -1.0f);
    if (pb.l < 0.0f) {
        pb = xyz_to_luv(rgb_to_xyz(triplet(1.0f, 1.0f, 0.0f)));
    }
    return pb;
}

static float mix_hue(float alpha, float h0, float h1)
{
    float M = std::fmod(pi + h1 - h0, twopi) - pi;
    return std::fmod(h0 + alpha * M, twopi);
}

static void get_color_points(float hue, float saturation, float warmth,
        triplet pb, float pb_hue, float pb_saturation,
        triplet* p0, triplet* p1, triplet* p2,
        triplet* q0, triplet* q1, triplet* q2)
{
    *p0 = lch_to_luv(triplet(0.0f, 0.0f, hue));
    *p1 = most_saturated_in_srgb(hue);
    triplet p2_lch;
    p2_lch.l = (1.0f - warmth) * 100.0f + warmth * pb.l;
    p2_lch.h = mix_hue(warmth, hue, pb_hue);
    p2_lch.c = lch_chroma(p2_lch.l, std::min(s_max(p2_lch.l, p2_lch.h), warmth * saturation * pb_saturation));
    *p2 = lch_to_luv(p2_lch);
    *q0 = (1.0f - saturation) * (*p0) + saturation * (*p1);
    *q2 = (1.0f - saturation) * (*p2) + saturation * (*p1);
    *q1 = 0.5f * ((*q0) + (*q2));
}

static triplet b(triplet b0, triplet b1, triplet b2, float t)
{
    float a = (1.0f - t) * (1.0f - t);
    float b = 2.0f * (1.0f - t) * t;
    float c = t * t;
    return a * b0 + b * b1 + c * b2;
}

static float inv_b(float b0, float b1, float b2, float v)
{
    return (b0 - b1 + std::sqrt(std::max(b1 * b1 - b0 * b2 + (b0 - 2.0f * b1 + b2) * v, 0.0f)))
        / (b0 - 2.0f * b1 + b2);
}

static triplet get_colormap_entry(float t,
        triplet p0, triplet p2,
        triplet q0, triplet q1, triplet q2,
        float contrast, float brightness)
{
    float l = 125 - 125 * std::pow(0.2f, (1.0f - contrast) * brightness + t * contrast);
    float T = (l <= q1.l ? 0.5f * inv_b(p0.l, q0.l, q1.l, l) : 0.5f * inv_b(q1.l, q2.l, p2.l, l) + 0.5f);
    return (T <= 0.5f ? b(p0, q0, q1, 2.0f * T) : b(q1, q2, p2, 2.0f * (T - 0.5f)));
}

/* Brewer-like color maps */

float BrewerSequentialDefaultContrastForSmallN(int n)
{
    return std::min(0.88f, 0.34f + 0.06f * n);
}

int BrewerSequential(int n, unsigned char* colormap, float hue,
        float contrast, float saturation, float brightness, float warmth)
{
    triplet pb, p0, p1, p2, q0, q1, q2;
    pb = get_bright_point();
    triplet pb_lch = luv_to_lch(pb);
    float pbs = lch_saturation(pb_lch.l, pb_lch.c);
    get_color_points(hue, saturation, warmth, pb, pb_lch.h, pbs, &p0, &p1, &p2, &q0, &q1, &q2);

    int clipped = 0;
    for (int i = 0; i < n; i++) {
        float t = (i + 0.5f) / n;
        triplet c = get_colormap_entry(t, p0, p2, q0, q1, q2, contrast, brightness);
        if (luv_to_colormap(c, colormap + 3 * i))
            clipped++;
    }
    return clipped;
}

float BrewerDivergingDefaultContrastForSmallN(int n)
{
    return std::min(0.88f, 0.34f + 0.06f * n);
}

int BrewerDiverging(int n, unsigned char* colormap, float hue, float divergence,
        float contrast, float saturation, float brightness, float warmth)
{
    float hue1 = hue + divergence;
    if (hue1 >= twopi)
        hue1 -= twopi;

    triplet pb;
    triplet p00, p01, p02, q00, q01, q02;
    triplet p10, p11, p12, q10, q11, q12;
    pb = get_bright_point();
    triplet pb_lch = luv_to_lch(pb);
    float pbs = lch_saturation(pb_lch.l, pb_lch.c);
    get_color_points(hue,  saturation, warmth, pb, pb_lch.h, pbs, &p00, &p01, &p02, &q00, &q01, &q02);
    get_color_points(hue1, saturation, warmth, pb, pb_lch.h, pbs, &p10, &p11, &p12, &q10, &q11, &q12);

    int clipped = 0;
    for (int i = 0; i < n; i++) {
        triplet c;
        if (n % 2 == 1 && i == n / 2) {
            // compute neutral color in the middle of the map
            triplet c0 = get_colormap_entry(1.0f, p00, p02, q00, q01, q02, contrast, brightness);
            triplet c1 = get_colormap_entry(1.0f, p10, p12, q10, q11, q12, contrast, brightness);
            if (n <= 9) {
                // for discrete color maps, use an extra neutral color
                float c0s = luv_saturation(c0);
                float c1s = luv_saturation(c1);
                float sn = 0.5f * (c0s + c1s) * warmth;
                c.l = 0.5f * (c0.l + c1.l);
                float cc = lch_chroma(c.l, std::min(s_max(c.l, pb_lch.h), sn));
                c = lch_to_luv(triplet(c.l, cc, pb_lch.h));
            } else {
                // for continuous color maps, use an average, since the extra neutral color looks bad
                c = 0.5f * (c0 + c1);
            }
        } else {
            float t = (i + 0.5f) / n;
            if (i < n / 2) {
                float tt = 2.0f * t;
                c = get_colormap_entry(tt, p00, p02, q00, q01, q02, contrast, brightness);
            } else {
                float tt = 2.0f * (1.0f - t);
                c = get_colormap_entry(tt, p10, p12, q10, q11, q12, contrast, brightness);
            }
        }
        if (luv_to_colormap(c, colormap + 3 * i))
            clipped++;
    }
    return clipped;
}

int BrewerQualitative(int n, unsigned char* colormap, float hue, float divergence,
        float contrast, float saturation, float brightness)
{
    // Get all information about yellow
    static triplet ylch(-1.0f, -1.0f, -1.0f);
    if (ylch.l < 0.0f) {
        ylch = luv_to_lch(xyz_to_luv(rgb_to_xyz(triplet(1.0f, 1.0f, 0.0f))));
    }

    // Get saturation of red (maximum possible saturation)
    static float rs = -1.0f;
    if (rs < 0.0f) {
        triplet rluv = xyz_to_luv(rgb_to_xyz(triplet(1.0f, 0.0f, 0.0f)));
        rs = luv_saturation(rluv);
    }

    // Derive parameters of the method
    float eps = hue / twopi;
    float r = divergence / twopi;
    float l0 = brightness * ylch.l;
    float l1 = (1.0f - contrast) * l0;

    // Generate colors
    int clipped = 0;
    for (int i = 0; i < n; i++) {
        float t = (i + 0.5f) / n;
        float ch = std::fmod(twopi * (eps + t * r), twopi);
        float alpha = hue_diff(ch, ylch.h) / pi;
        float cl = (1.0f - alpha) * l0 + alpha * l1;
        float cs = std::min(s_max(cl, ch), saturation * rs);
        triplet c = lch_to_luv(triplet(cl, lch_chroma(cl, cs), ch));
        if (luv_to_colormap(c, colormap + 3 * i))
            clipped++;
    }
    return clipped;
}

/* Perceptually uniform (PU) */

static triplet lch_compute_uniform_lc(float t,
        float t0, float t1,
        triplet lch0, triplet lch1, float D,
        float hue)
{
    // t is in [0,1]; s is in [0,1] but relative to [t0,t1]
    float s = (t - t0) / (t1 - t0);

    triplet lcht;
    lcht.h = hue;
    lcht.l = (1.0f - s) * lch0.l + s * lch1.l;

    // compute four solutions for lcht.c based on two conditions:
    float lcht_cs[4];
    // 1) the distance between lcht and lch0 is (s * D)
    float tmp00 = lch0.c * std::cos(lcht.h - lch0.h);
    float tmp01 = std::max(0.0f, sqr(tmp00) - sqr(lcht.l - lch0.l) - sqr(lch0.c) + sqr(s * D));
    lcht_cs[0] = tmp00 + std::sqrt(tmp01);
    lcht_cs[1] = tmp00 - std::sqrt(tmp01);
    // 2) the distance between lcht and lch1 is ((1-s) * D)
    float tmp10 = lch1.c * std::cos(lcht.h - lch1.h);
    float tmp11 = std::max(0.0f, sqr(tmp10) - sqr(lcht.l - lch1.l) - sqr(lch1.c) + sqr((1.0f - s) * D));
    lcht_cs[2] = tmp10 + std::sqrt(tmp11);
    lcht_cs[3] = tmp10 - std::sqrt(tmp11);

    // find the best solution
    float min_c = std::min(lch0.c, lch1.c);
    float max_c = std::max(lch0.c, lch1.c);
    float min_solution_error = 9999.9f;
    int min_solution_error_index = -1;
    for (int i = 0; i < 4; i++) {
        if (lcht_cs[i] < min_c || lcht_cs[i] > max_c) {
            // this solution is invalid
        } else {
            // solution error is sum of the deviations from condition 1 and 2
            float dist_lch0_lcht = lch_distance(lch0, triplet(lcht.l, lcht_cs[i], lcht.h));
            float dist_lch1_lcht = lch_distance(lch1, triplet(lcht.l, lcht_cs[i], lcht.h));
            float solution_error = std::abs(dist_lch0_lcht - s * D) + std::abs(dist_lch1_lcht - (1.0f - s) * D);
            //fprintf(stderr, "t=%g: VALID %i with error %g\n", t, i, solution_error);
            if (min_solution_error_index == -1 || solution_error < min_solution_error) {
                min_solution_error = solution_error;
                min_solution_error_index = i;
            }
        }
    }
    if (min_solution_error_index == -1) {
        //fprintf(stderr, "TODO: FALLBACK at t=%g\n", t);
        lcht.c = 0.5f * (lch0.c + lch1.c);
    } else {
        lcht.c = lcht_cs[min_solution_error_index];
    }

    return lcht;
}

int PUSequentialLightness(int n, unsigned char* colormap,
        float lightness_range, float saturation_range, float saturation, float hue)
{
    triplet lch_00, lch_10, lch_05;

    lch_00.l = (1.0f - lightness_range) * 100.0f;
    lch_00.c = lch_chroma(lch_00.l, 1.0f - saturation_range);
    lch_00.h = hue;
    lch_10.l = lightness_range * 100.0f;
    lch_10.c = lch_chroma(lch_10.l, 1.0f - saturation_range);
    lch_10.h = hue;
    lch_05.l = (1.0f - 0.5f) * lch_00.l + 0.5f * lch_10.l;
    lch_05.c = lch_chroma(lch_05.l, 5.0f * saturation_range * saturation);
    lch_05.h = hue;

    // the following distances are actually the same since hue is constant:
    float D_00_05 = lch_distance(lch_00, lch_05);
    float D_05_10 = lch_distance(lch_05, lch_10);

    triplet lch;
    int clipped = 0;
    for (int i = 0; i < n; i++) {
        float t = (i + 0.5f) / n;
        if (t <= 0.5f) {
            lch = lch_compute_uniform_lc(t, 0.0f, 0.5f, lch_00, lch_05, D_00_05, hue);
        } else {
            lch = lch_compute_uniform_lc(t, 0.5f, 1.0f, lch_05, lch_10, D_05_10, hue);
        }
        bool is_clipped = false;
        if (lch.c > 100.0f) {
            is_clipped = true;
            lch.c = 100.0f;
        }
        if (lch_to_colormap(lch, colormap + 3 * i)) {
            is_clipped = true;
        }
        if (is_clipped)
            clipped++;
    }
    return clipped;
}

int PUSequentialSaturation(int n, unsigned char* colormap,
        float saturation_range, float lightness, float saturation, float hue)
{
    lightness = std::max(0.01f, lightness * 100.0f);

    triplet lch_00, lch_10;
    lch_00.l = lightness;
    lch_00.c = lch_chroma(lch_00.l, 1.0f - saturation_range);
    lch_00.h = hue;
    lch_10.l = lightness;
    lch_10.c = lch_chroma(lch_10.l, saturation_range * 5.0f * saturation);
    lch_10.h = hue;

    float D_00_10 = lch_distance(lch_00, lch_10);

    triplet lch;
    int clipped = 0;
    for (int i = 0; i < n; i++) {
        float t = (i + 0.5f) / n;
        lch = lch_compute_uniform_lc(t, 0.0f, 1.0f, lch_00, lch_10, D_00_10, hue);
        bool is_clipped = false;
        if (lch.c > 100.0f) {
            is_clipped = true;
            lch.c = 100.0f;
        }
        if (lch_to_colormap(lch, colormap + 3 * i)) {
            is_clipped = true;
        }
        if (is_clipped)
            clipped++;
    }
    return clipped;
}

int PUSequentialRainbow(int n, unsigned char* colormap,
        float lightness_range, float saturation_range,
        float hue, float rotations, float saturation)
{
    triplet lch_00, lch_10, lch_05;

    lch_00.l = (1.0f - lightness_range) * 100.0f;
    lch_00.c = lch_chroma(lch_00.l, 1.0f - saturation_range);
    lch_00.h = hue + 0.0f * rotations * twopi;
    lch_10.l = lightness_range * 100.0f;
    lch_10.c = lch_chroma(lch_10.l, 1.0f - saturation_range);
    lch_10.h = hue + 1.0f * rotations * twopi;
    lch_05.l = 0.5f * (lch_00.l + lch_10.l);
    lch_05.c = lch_chroma(lch_05.l, saturation_range * saturation);
    lch_05.h = hue + 0.5f * rotations * twopi;

    // the following are not necessarily equal because hue varies:
    float D_00_05 = lch_distance(lch_00, lch_05);
    float D_05_10 = lch_distance(lch_05, lch_10);

    triplet lch;
    int clipped = 0;
    for (int i = 0; i < n; i++) {
        float t = (i + 0.5f) / n;
        float h = hue + t * rotations * twopi;
        if (t <= 0.5f) {
            lch = lch_compute_uniform_lc(t, 0.0f, 0.5f, lch_00, lch_05, D_00_05, h);
        } else {
            lch = lch_compute_uniform_lc(t, 0.5f, 1.0f, lch_05, lch_10, D_05_10, h);
        }
        bool is_clipped = false;
        if (lch.c > 100.0f) {
            is_clipped = true;
            lch.c = 100.0f;
        }
        if (lch_to_colormap(lch, colormap + 3 * i)) {
            is_clipped = true;
        }
        if (is_clipped)
            clipped++;
    }
    return clipped;
}

static float plancks_law(float temperature, float lambda)
{
    const float c = 299792458.0f;     // speed of light in vacuum
    //const float c = 299700000.0f;     // speed of visible light in air
    const float h = 6.626070041e-34f; // Planck's constant
    const float k = 1.38064853e-23f;  // Boltzmann constant
    return 2.0f * h * c * c * std::pow(lambda, -5.0f)
        / (std::exp(h * c / (lambda * k * temperature)) - 1.0f);

}

#if 0
static triplet color_matching_function_approx(float lambda /* in nanometers */)
{
    // Analytic approximation of the CIE 1964 10 deg standard observer CMF.
    // Taken from the paper "Simple Analytic Approximations to the CIE XYZ
    // Color Matching Functions" by Wyman, Sloan, Shirley, JCGT 2(2), 2013.
    triplet xyz;
    xyz.x = 0.398f * std::exp(-1250.0f * sqr(std::log((lambda + 570.1f) / 1014.0f)))
          + 1.132f * std::exp(-234.0f  * sqr(std::log((1338.0f - lambda) / 743.5f)));
    xyz.y = 1.011f * std::exp(-0.5f    * sqr((lambda - 556.1f) / 46.14f));
    xyz.z = 2.060f * std::exp(-32.0f   * sqr(std::log((lambda - 265.8f) / 180.4f)));
    return xyz;
}
#endif

static triplet color_matching_function(int lambda /* in nanometers */)
{
    // Tables in 5nm resolution for the CIE 1931 2 deg Standard Observer CMF,
    // from lambda=380 to lambda=780. Obtained from
    // http://www.cie.co.at/publ/abst/datatables15_2004/CIE_sel_colorimetric_tables.xls
    // on 2016-02-10.
    static const triplet cmf_xyz[] = {
        triplet(0.001368, 0.000039, 0.006450),
        triplet(0.002236, 0.000064, 0.010550),
        triplet(0.004243, 0.000120, 0.020050),
        triplet(0.007650, 0.000217, 0.036210),
        triplet(0.014310, 0.000396, 0.067850),
        triplet(0.023190, 0.000640, 0.110200),
        triplet(0.043510, 0.001210, 0.207400),
        triplet(0.077630, 0.002180, 0.371300),
        triplet(0.134380, 0.004000, 0.645600),
        triplet(0.214770, 0.007300, 1.039050),
        triplet(0.283900, 0.011600, 1.385600),
        triplet(0.328500, 0.016840, 1.622960),
        triplet(0.348280, 0.023000, 1.747060),
        triplet(0.348060, 0.029800, 1.782600),
        triplet(0.336200, 0.038000, 1.772110),
        triplet(0.318700, 0.048000, 1.744100),
        triplet(0.290800, 0.060000, 1.669200),
        triplet(0.251100, 0.073900, 1.528100),
        triplet(0.195360, 0.090980, 1.287640),
        triplet(0.142100, 0.112600, 1.041900),
        triplet(0.095640, 0.139020, 0.812950),
        triplet(0.057950, 0.169300, 0.616200),
        triplet(0.032010, 0.208020, 0.465180),
        triplet(0.014700, 0.258600, 0.353300),
        triplet(0.004900, 0.323000, 0.272000),
        triplet(0.002400, 0.407300, 0.212300),
        triplet(0.009300, 0.503000, 0.158200),
        triplet(0.029100, 0.608200, 0.111700),
        triplet(0.063270, 0.710000, 0.078250),
        triplet(0.109600, 0.793200, 0.057250),
        triplet(0.165500, 0.862000, 0.042160),
        triplet(0.225750, 0.914850, 0.029840),
        triplet(0.290400, 0.954000, 0.020300),
        triplet(0.359700, 0.980300, 0.013400),
        triplet(0.433450, 0.994950, 0.008750),
        triplet(0.512050, 1.000000, 0.005750),
        triplet(0.594500, 0.995000, 0.003900),
        triplet(0.678400, 0.978600, 0.002750),
        triplet(0.762100, 0.952000, 0.002100),
        triplet(0.842500, 0.915400, 0.001800),
        triplet(0.916300, 0.870000, 0.001650),
        triplet(0.978600, 0.816300, 0.001400),
        triplet(1.026300, 0.757000, 0.001100),
        triplet(1.056700, 0.694900, 0.001000),
        triplet(1.062200, 0.631000, 0.000800),
        triplet(1.045600, 0.566800, 0.000600),
        triplet(1.002600, 0.503000, 0.000340),
        triplet(0.938400, 0.441200, 0.000240),
        triplet(0.854450, 0.381000, 0.000190),
        triplet(0.751400, 0.321000, 0.000100),
        triplet(0.642400, 0.265000, 0.000050),
        triplet(0.541900, 0.217000, 0.000030),
        triplet(0.447900, 0.175000, 0.000020),
        triplet(0.360800, 0.138200, 0.000010),
        triplet(0.283500, 0.107000, 0.000000),
        triplet(0.218700, 0.081600, 0.000000),
        triplet(0.164900, 0.061000, 0.000000),
        triplet(0.121200, 0.044580, 0.000000),
        triplet(0.087400, 0.032000, 0.000000),
        triplet(0.063600, 0.023200, 0.000000),
        triplet(0.046770, 0.017000, 0.000000),
        triplet(0.032900, 0.011920, 0.000000),
        triplet(0.022700, 0.008210, 0.000000),
        triplet(0.015840, 0.005723, 0.000000),
        triplet(0.011359, 0.004102, 0.000000),
        triplet(0.008111, 0.002929, 0.000000),
        triplet(0.005790, 0.002091, 0.000000),
        triplet(0.004109, 0.001484, 0.000000),
        triplet(0.002899, 0.001047, 0.000000),
        triplet(0.002049, 0.000740, 0.000000),
        triplet(0.001440, 0.000520, 0.000000),
        triplet(0.001000, 0.000361, 0.000000),
        triplet(0.000690, 0.000249, 0.000000),
        triplet(0.000476, 0.000172, 0.000000),
        triplet(0.000332, 0.000120, 0.000000),
        triplet(0.000235, 0.000085, 0.000000),
        triplet(0.000166, 0.000060, 0.000000),
        triplet(0.000117, 0.000042, 0.000000),
        triplet(0.000083, 0.000030, 0.000000),
        triplet(0.000059, 0.000021, 0.000000),
        triplet(0.000042, 0.000015, 0.000000),
    };
    triplet xyz(0, 0, 0);
    if (lambda >= 380 && lambda <= 780) {
        int i = (lambda - 380) / 5;
        xyz = cmf_xyz[i];
        if (lambda % 5 != 0) {
            int i1 = i + 1;
            triplet xyz1 = cmf_xyz[i1];
            float alpha = (lambda % 5) / 5.0f;
            xyz = (1.0f - alpha) * xyz + alpha * xyz1;
        }
    }
    return xyz;
}

static float black_body_hue_at_temperature(float t)
{
    // Integrate radiance over the visible spectrum; according
    // to literature, sampling at 10nm intervals is enough.
    triplet xyz(0, 0, 0);
    int stepsize = 5;
    float s = stepsize * 1e-9f; // stepsize in meters
    for (int lambda = 360; lambda <= 830; lambda += stepsize) {
        float l = lambda * 1e-9f; // lambda in in meters
        float radiosity = pi * plancks_law(t, l);
        //xyz = xyz + s * radiosity * color_matching_function_approx(l);
        xyz = xyz + s * radiosity * color_matching_function(lambda);
    }
    triplet lch = luv_to_lch(xyz_to_luv(adjust_y(xyz, 50.0f)));
    return lch.h;
}

int PUSequentialBlackBody(int n, unsigned char* colormap,
        float temperature, float temperature_range,
        float lightness_range, float saturation_range, float saturation)
{
    triplet lch_00, lch_10, lch_05;

    lch_00.l = (1.0f - lightness_range) * 100.0f;
    lch_00.c = lch_chroma(lch_00.l, 1.0f - saturation_range);
    lch_00.h = black_body_hue_at_temperature(temperature + 0.0f * temperature_range);
    lch_10.l = lightness_range * 100.0f;
    lch_10.c = lch_chroma(lch_10.l, 1.0f - saturation_range);
    lch_10.h = black_body_hue_at_temperature(temperature + 1.0f * temperature_range);
    lch_05.l = 0.5f * (lch_00.l + lch_10.l);
    lch_05.c = lch_chroma(lch_05.l, saturation_range * saturation);
    lch_05.h = black_body_hue_at_temperature(temperature + 0.5f * temperature_range);

    // the following are not necessarily equal because hue varies:
    float D_00_05 = lch_distance(lch_00, lch_05);
    float D_05_10 = lch_distance(lch_05, lch_10);

    triplet lch;
    int clipped = 0;
    for (int i = 0; i < n; i++) {
        float t = (i + 0.5f) / n;
        float h = black_body_hue_at_temperature(temperature + t * temperature_range);
        if (t <= 0.5f) {
            lch = lch_compute_uniform_lc(t, 0.0f, 0.5f, lch_00, lch_05, D_00_05, h);
        } else {
            lch = lch_compute_uniform_lc(t, 0.5f, 1.0f, lch_05, lch_10, D_05_10, h);
        }
        bool is_clipped = false;
        if (lch.c > 100.0f) {
            is_clipped = true;
            lch.c = 100.0f;
        }
        if (lch_to_colormap(lch, colormap + 3 * i)) {
            is_clipped = true;
        }
        if (is_clipped)
            clipped++;
    }
    return clipped;
}

static float multi_hue_get(float t, int hues, const float* hue_values, const float* hue_positions)
{
    /* Trivial and pathological cases */
    if (hues < 1)
        return 0.0f;
    if (hues == 1)
        return hue_values[0];
    if (t <= hue_positions[0])
        return hue_values[0];
    if (t >= hue_positions[hues - 1])
        return hue_values[hues - 1];
    /* Find index i so that t is in [hue_positions[i], hue_positions[i+1]]  */
    int i;
    for (i = 0; i < hues - 2; i++) {
        if (t >= hue_positions[i] && t < hue_positions[i+1])
            break;
    }
    float h0 = hue_values[i];
    float h1 = hue_values[i + 1];
    float p0 = hue_positions[i];
    float p1 = hue_positions[i + 1];
    /* Check if the distance between h0 and h1 is shorter if we pass the
     * boundary at 2pi */
    if (h0 < h1 && h1 - h0 > h0 + twopi - h1) {
        h0 += twopi;
    } else if (h1 < h0 && h0 - h1 > h1 + twopi - h0) {
        h1 += twopi;
    }
    float alpha = (t - p0) / (p1 - p0);
    float hue = (1.0f - alpha) * h0 + alpha * h1;
    if (hue >= twopi)
        hue -= twopi;
    return hue;
}

int PUSequentialMultiHue(int n, unsigned char* colormap,
        float lightness_range,
        float saturation_range,
        float saturation,
        int hues,
        const float* hue_values,
        const float* hue_positions)
{
    triplet lch_00, lch_10, lch_05;

    lch_00.l = (1.0f - lightness_range) * 100.0f;
    lch_00.c = lch_chroma(lch_00.l, 1.0f - saturation_range);
    lch_00.h = multi_hue_get(0.0f, hues, hue_values, hue_positions);
    lch_10.l = lightness_range * 100.0f;
    lch_10.c = lch_chroma(lch_10.l, 1.0f - saturation_range);
    lch_10.h = multi_hue_get(1.0f, hues, hue_values, hue_positions);
    lch_05.l = (1.0f - 0.5f) * lch_00.l + 0.5f * lch_10.l;
    lch_05.c = lch_chroma(lch_05.l, 5.0f * saturation_range * saturation);
    lch_05.h = multi_hue_get(0.5f, hues, hue_values, hue_positions);

    // the following distances should ideally be the same,
    // but they are usually not since we use different hues.
    // at least they should be close since hue differences
    // are less dominant in the distance measure.
    float D_00_05 = lch_distance(lch_00, lch_05);
    float D_05_10 = lch_distance(lch_05, lch_10);

    triplet lch;
    int clipped = 0;
    for (int i = 0; i < n; i++) {
        float t = (i + 0.5f) / n;
        float h = multi_hue_get(t, hues, hue_values, hue_positions);
        if (t <= 0.5f) {
            lch = lch_compute_uniform_lc(t, 0.0f, 0.5f, lch_00, lch_05, D_00_05, h);
        } else {
            lch = lch_compute_uniform_lc(t, 0.5f, 1.0f, lch_05, lch_10, D_05_10, h);
        }
        bool is_clipped = false;
        if (lch.c > 100.0f) {
            is_clipped = true;
            lch.c = 100.0f;
        }
        if (lch_to_colormap(lch, colormap + 3 * i)) {
            is_clipped = true;
        }
        if (is_clipped)
            clipped++;
    }
    return clipped;
}

int PUDivergingLightness(int n, unsigned char* colormap,
        float lightness_range, float saturation_range, float saturation, float hue, float divergence)
{
    int lowerN = n / 2;
    int higherN = n - lowerN;
    int clipped = 0;

    clipped += PUSequentialLightness(higherN, colormap, lightness_range, saturation_range, saturation, hue + divergence);
    for (int i = 0; i < higherN; i++)
        for (int j = 0; j < 3; j++)
            colormap[3 * lowerN + 3 * i + j] = colormap[3 * (higherN - 1 - i) + j];
    clipped += PUSequentialLightness(lowerN, colormap, lightness_range, saturation_range, saturation, hue);
    return clipped;
}

int PUDivergingSaturation(int n, unsigned char* colormap,
        float saturation_range, float lightness, float saturation, float hue, float divergence)
{
    int lowerN = n / 2;
    int higherN = n - lowerN;
    int clipped = 0;

    clipped += PUSequentialSaturation(lowerN, colormap + 3 * lowerN, saturation_range, lightness, saturation, hue);
    for (int i = 0; i < lowerN; i++)
        for (int j = 0; j < 3; j++)
            colormap[3 * i + j] = colormap[3 * lowerN + 3 * (lowerN - 1 - i) + j];
    clipped += PUSequentialSaturation(higherN, colormap + 3 * lowerN, saturation_range, lightness, saturation, hue + divergence);
    return clipped;
}

int PUQualitativeHue(int n, unsigned char* colormap,
        float hue, float divergence, float lightness, float saturation)
{
    divergence *= (n - 1.0f) / n;
    triplet lch;
    lch.l = std::max(0.01f, lightness * 100.0f);
    lch.c = lch_chroma(lch.l, saturation * 5.0f);
    bool all_clipped = (lch.c > 100.0f);
    if (all_clipped)
        lch.c = 100.0f;
    int clipped = 0;
    for (int i = 0; i < n; i++) {
        float t = (i + 0.5f) / n;
        lch.h = hue + t * divergence;
        bool is_clipped = false;
        if (lch_to_colormap(lch, colormap + 3 * i))
            is_clipped = true;
        if (is_clipped || all_clipped)
            clipped++;
    }
    return clipped;
}

/* CubeHelix */

int CubeHelix(int n, unsigned char* colormap, float hue,
        float rot, float saturation, float gamma)
{
    int clipped = 0;
    for (int i = 0; i < n; i++) {
        float fract = (i + 0.5f) / n;
        float angle = twopi * (hue / 3.0f + 1.0f + rot * fract);
        fract = std::pow(fract, gamma);
        float amp = saturation * fract * (1.0f - fract) / 2.0f;
        float s = std::sin(angle);
        float c = std::cos(angle);
        triplet srgb(
                fract + amp * (-0.14861f * c + 1.78277f * s),
                fract + amp * (-0.29227f * c - 0.90649f * s),
                fract + amp * (1.97294f * c));
        bool clipped_[3];
        colormap[3 * i + 0] = float_to_uchar(srgb.r, clipped_ + 0);
        colormap[3 * i + 1] = float_to_uchar(srgb.g, clipped_ + 1);
        colormap[3 * i + 2] = float_to_uchar(srgb.b, clipped_ + 2);
        if (clipped_[0] || clipped_[1] || clipped_[2])
            clipped++;
    }
    return clipped;
}

/* Moreland */

static triplet lab_to_msh(triplet lab)
{
    triplet msh;
    msh.m = std::sqrt(lab.l * lab.l + lab.a * lab.a + lab.b * lab.b);
    msh.s = (msh.m > 0.001f) ? std::acos(lab.l / msh.m) : 0.0f;
    msh.h = (msh.s > 0.001f) ? std::atan2(lab.b, lab.a) : 0.0f;
    return msh;
}

static triplet msh_to_lab(triplet msh)
{
    return triplet(
            msh.m * std::cos(msh.s),
            msh.m * std::sin(msh.s) * std::cos(msh.h),
            msh.m * std::sin(msh.s) * std::sin(msh.h));
}

static float adjust_hue(triplet msh, float unsaturated_m)
{
    if (msh.m >= unsaturated_m - 0.1f) {
        return msh.h;
    } else {
        float hue_spin = msh.s * std::sqrt(unsaturated_m * unsaturated_m - msh.m * msh.m)
            / (msh.m * std::sin(msh.s));
        if (msh.h > -pi / 3.0f)
            return msh.h + hue_spin;
        else
            return msh.h - hue_spin;
    }
}

int Moreland(int n, unsigned char* colormap,
        unsigned char sr0, unsigned char sg0, unsigned char sb0,
        unsigned char sr1, unsigned char sg1, unsigned char sb1)
{
    triplet omsh0 = lab_to_msh(xyz_to_lab(rgb_to_xyz(srgb_to_rgb(triplet(
                            uchar_to_float(sr0), uchar_to_float(sg0), uchar_to_float(sb0))))));
    triplet omsh1 = lab_to_msh(xyz_to_lab(rgb_to_xyz(srgb_to_rgb(triplet(
                            uchar_to_float(sr1), uchar_to_float(sg1), uchar_to_float(sb1))))));
    bool place_white = (omsh0.s >= 0.05f && omsh1.s >= 0.05f && hue_diff(omsh0.h, omsh1.h) > pi / 3.0f);
    float mmid = std::max(std::max(omsh0.m, omsh1.m), 88.0f);

    int clipped = 0;
    for (int i = 0; i < n; i++) {
        triplet msh0 = omsh0;
        triplet msh1 = omsh1;
        float t = (i + 0.5f) / n;
        if (place_white) {
            if (t < 0.5f) {
                msh1.m = mmid;
                msh1.s = 0.0f;
                msh1.h = 0.0f;
                t *= 2.0f;
            } else {
                msh0.m = mmid;
                msh0.s = 0.0f;
                msh0.h = 0.0f;
                t = 2.0f * t - 1.0f;
            }
        }
        if (msh0.s < 0.05f && msh1.s >= 0.05f) {
            msh0.h = adjust_hue(msh1, msh0.m);
        } else if (msh0.s >= 0.05f && msh1.s < 0.05f) {
            msh1.h = adjust_hue(msh0, msh1.m);
        }
        triplet msh = (1.0f - t) * msh0 + t * msh1;
        if (lab_to_colormap(msh_to_lab(msh), colormap + 3 * i))
            clipped++;
    }
    return clipped;
}

/* McNames */

static void cart2pol(float x, float y, float* theta, float* rho)
{
    *theta = std::atan2(y, x);
    *rho = std::hypot(x, y);
}

static void pol2cart(float theta, float rho, float* x, float* y)
{
    *x = rho * std::cos(theta);
    *y = rho * std::sin(theta);
}

static float windowfunc(float t)
{
    static const float ww = std::sqrt(3.0f / 8.0f);
    /* triangular window function: */
#if 0
    if (t <= 0.5f) {
        return ww * 2.0f * t;
    } else {
        return ww * 2.0f * (1.0f - t);
    }
#endif
    /* window function based on cosh: */
#if 1
    static const float acosh2 = std::acosh(2.0f);
    return 0.95f * ww * (2.0f - std::cosh(acosh2 * (2.0f * t - 1.0f)));
#endif
}

int McNames(int n, unsigned char* colormap, float periods)
{
    static const float sqrt3 = std::sqrt(3.0f);
    static const float a12 = std::asin(1.0f / sqrt3);
    static const float a23 = pi / 4.0f;

    int clipped = 0;
    for (int i = 0; i < n; i++) {
        float t = 1.0f - (i + 0.5f) / n;
        float w = windowfunc(t);
        float tt = (1.0f - t) * sqrt3;
        float ttt = (tt - sqrt3 / 2.0f) * periods * twopi / sqrt3;

        float r0, g0, b0, r1, g1, b1, r2, g2, b2;
        float ag, rd;
        r0 = tt;
        g0 = w * std::cos(ttt);
        b0 = w * std::sin(ttt);
        cart2pol(r0, g0, &ag, &rd);
        pol2cart(ag + a12, rd, &r1, &g1);
        b1 = b0;
        cart2pol(r1, b1, &ag, &rd);
        pol2cart(ag + a23, rd, &r2, &b2);
        g2 = g1;

        bool clipped_[3];
        colormap[3 * i + 0] = float_to_uchar(r2, clipped_ + 0);
        colormap[3 * i + 1] = float_to_uchar(g2, clipped_ + 1);
        colormap[3 * i + 2] = float_to_uchar(b2, clipped_ + 2);
        if (clipped_[0] || clipped_[1] || clipped_[2])
            clipped++;
    }
    return clipped;
}




std::string ToCSV(int n, const unsigned char* srgb_colormap)
{
    std::string s;
    for (int i = 0; i < n; i++) {
        s +=  std::to_string(srgb_colormap[3 * i + 0]) + ", "
            + std::to_string(srgb_colormap[3 * i + 1]) + ", "
            + std::to_string(srgb_colormap[3 * i + 2]) + '\n';
    }
    return s;
}

std::string ToJSON(int n, const unsigned char* srgb_colormap)
{
    std::string s =
        "[\n"
        "{\n"
        "\"ColorSpace\" : \"RGB\",\n"
        "\"Name\" : \"GenColorMapGenerated\",\n"
        "\"NanColor\" : [ -1, -1, -1 ],\n"
        "\"RGBPoints\" : [\n";
    for (int i = 0; i < n; i++) {
        s +=  std::to_string(i / float(n - 1)) + ", "
            + std::to_string(srgb_colormap[3 * i + 0] / 255.0f) + ", "
            + std::to_string(srgb_colormap[3 * i + 1] / 255.0f) + ", "
            + std::to_string(srgb_colormap[3 * i + 2] / 255.0f) + (i == n - 1 ? "\n" : ",\n");
    }
    s += "]\n}\n]\n";
    return s;
}

std::string ToPPM(int n, const unsigned char* srgb_colormap)
{
    std::string s = "P3\n"; // magic number for plain PPM
    s += std::to_string(n) + " 1\n"; // width and height
    s += "255\n"; // max val
    for (int i = 0; i < n; i++) {
        s +=  std::to_string(srgb_colormap[3 * i + 0]) + ' '
            + std::to_string(srgb_colormap[3 * i + 1]) + ' '
            + std::to_string(srgb_colormap[3 * i + 2]) + '\n';
    }
    return s;
}









int main(int argc, char* argv[])
{
    bool print_version = false;
    bool print_help = false;
    int format = csv; // output format
    int type = brewer_seq;
    int n = 256; // number of colors in the map
    float hue = -1.0f;
    float divergence = -1.0f;
    float contrast = -1.0f;
    float saturation = -1.0f;
    float saturation_range = -1.0f;
    float brightness = -1.0f;
    float warmth = -1.0f;
    float lightness = -1.0f;
    float lightness_range = -1.0f;
    float rotations = NAN;
    float temperature = -1.0f;
    float temperature_range = -1.0f;
    std::vector<float> hue_values;
    std::vector<float> hue_positions;
    float gamma = -1.0f;
    bool have_color0 = false;
    unsigned char color0[3];
    bool have_color1 = false;
    unsigned char color1[3];
    float periods = NAN;
    
    
    // -option  
    struct option options[] = {
        { "version",           no_argument,       0, 'v' },
        { "help",              no_argument,       0, 'H' },
        { "format",            required_argument, 0, 'f' },
        { "type",              required_argument, 0, 't' },
        { "n",                 required_argument, 0, 'n' },
        { "hue",               required_argument, 0, 'h' },
        { "divergence",        required_argument, 0, 'd' },
        { "contrast",          required_argument, 0, 'c' },
        { "saturation",        required_argument, 0, 's' },
        { "saturation-range",  required_argument, 0, 'S' },
        { "brightness",        required_argument, 0, 'b' },
        { "warmth",            required_argument, 0, 'w' },
        { "lightness",         required_argument, 0, 'l' },
        { "lightness-range",   required_argument, 0, 'L' },
        { "rotations",         required_argument, 0, 'r' },
        { "temperature",       required_argument, 0, 'T' },
        { "temperature-range", required_argument, 0, 'R' },
        { "hue-values",        required_argument, 0, 'V' },
        { "hue-positions",     required_argument, 0, 'P' },
        { "gamma",             required_argument, 0, 'g' },
        { "color0",            required_argument, 0, 'A' },
        { "color1",            required_argument, 0, 'O' },
        { "periods",           required_argument, 0, 'p' },
        { 0, 0, 0, 0 }
    };

    for (;;) {
        int c = getopt_long(argc, argv, "vHf:t:n:h:d:c:s:S:b:w:l:L:r:T:R:V:P::g:A:O:p:", options, NULL);
        if (c == -1)
            break;
        switch (c) {
        case 'v': print_version = true; break;
        case 'H': print_help = true;  break;
        case 'f': format = (strcmp(optarg, "csv") == 0 ? csv
                          : strcmp(optarg, "json") == 0 ? json
                          : strcmp(optarg, "ppm") == 0 ? ppm
                          : -1); break;
        case 't':
            type = (  strcmp(optarg, "brewer-sequential") == 0 ? brewer_seq
                    : strcmp(optarg, "brewer-diverging") == 0 ? brewer_div
                    : strcmp(optarg, "brewer-qualitative") == 0 ? brewer_qual
                    : strcmp(optarg, "pusequential-lightness") == 0 ? puseq_lightness
                    : strcmp(optarg, "pusequential-saturation") == 0 ? puseq_saturation
                    : strcmp(optarg, "pusequential-rainbow") == 0 ? puseq_rainbow
                    : strcmp(optarg, "pusequential-blackbody") == 0 ? puseq_blackbody
                    : strcmp(optarg, "pusequential-multihue") == 0 ? puseq_multihue
                    : strcmp(optarg, "pudiverging-lightness") == 0 ? pudiv_lightness
                    : strcmp(optarg, "pudiverging-saturation") == 0 ? pudiv_saturation
                    : strcmp(optarg, "puqualitative-hue") == 0 ? puqual_hue
                    : strcmp(optarg, "cubehelix") == 0 ? cubehelix
                    : strcmp(optarg, "moreland") == 0 ? moreland
                    : strcmp(optarg, "mcnames") == 0 ? mcnames
                    : -1);  break;
        case 'n': n = atoi(optarg); break;
        case 'h': hue = atof(optarg) * M_PI / 180.0; break;
        case 'd': divergence = atof(optarg) * M_PI / 180.0; break;
        case 'c': contrast = atof(optarg); break;
        case 's': saturation = atof(optarg); break;
        case 'S': saturation_range = atof(optarg); break;
        case 'b': brightness = atof(optarg); break;
        case 'w': warmth = atof(optarg); break;
        case 'l': lightness = atof(optarg); break;
        case 'L': lightness_range = atof(optarg); break;
        case 'r': rotations = atof(optarg); break;
        case 'T': temperature = atof(optarg); break;
        case 'R': temperature_range = atof(optarg); break;
        case 'V':
            for (;;) {
                hue_values.push_back(atof(optarg) * M_PI / 180.0);
                optarg = strchr(optarg, ',');
                if (!optarg) break;
                optarg++;} break;
                
                
        case 'P':
            for (;;) {
                hue_positions.push_back(atof(optarg));
                optarg = strchr(optarg, ',');
                if (!optarg)
                    break;
                optarg++;
            }
            break;
        case 'g': gamma = atof(optarg); break;
        case 'A': std::sscanf(optarg, "%hhu,%hhu,%hhu", color0 + 0, color0 + 1, color0 + 2); have_color0 = true; break;
        case 'O': std::sscanf(optarg, "%hhu,%hhu,%hhu", color1 + 0, color1 + 1, color1 + 2); have_color1 = true; break;
        case 'p': periods = atof(optarg); break;
        default:  return 1;
        }
    }

    if (print_version) {
        printf("gencolormap version 2.3\n"
                "https://marlam.de/gencolormap\n"
                "Copyright (C) 2022 Computer Graphics Group, University of Siegen.\n"
                "Written by Martin Lambers <martin.lambers@uni-siegen.de>.\n"
                "This is free software under the terms of the MIT/Expat License.\n"
                "There is NO WARRANTY, to the extent permitted by law.\n");
        return 0;
    }

    if (print_help) {
        printf("Usage: %s [option...]\n"
                "Generates a color map and prints it to standard output.\n"
                "Prints the number of colors that had to be clipped to standard error.\n"
                "Common options:\n"
                "  [-f|--format=csv|json|ppm]          Set output format\n"
                "  [-n|--n=N]                          Set number of colors in the map\n"
                "Brewer-like color maps:\n"
                "  [-t|--type=brewer-sequential]       Generate a sequential color map\n"
                "  [-t|--type=brewer-diverging]        Generate a diverging color map\n"
                "  [-t|--type=brewer-qualitative]      Generate a qualitative color map\n"
                "  [-h|--hue=H]                        Set default hue in [0,360] degrees\n"
                "  [-c|--contrast=C]                   Set contrast in [0,1]\n"
                "  [-s|--saturation=S]                 Set saturation in [0,1]\n"
                "  [-b|--brightness=B]                 Set brightness in [0,1]\n"
                "  [-w|--warmth=W]                     Set warmth in [0,1] for seq. and div. maps\n"
                "  [-d|--divergence=D]                 Set diverg. in deg for div. and qual. maps\n"
                "Perceptually uniform color maps:\n"
                "  [-t|--type=pusequential-lightness]  Sequential map, varying lightness\n"
                "  [-t|--type=pusequential-saturation] Sequential map, varying saturation\n"
                "  [-t|--type=pusequential-rainbow]    Sequential map, varying hue (rainbow)\n"
                "  [-t|--type=pusequential-blackbody]  Sequential map, varying hue (black body)\n"
                "  [-t|--type=pusequential-multihue]   Sequential map, varying hue (custom)\n"
                "  [-t|--type=pudiverging-lightness]   Diverging map, varying lightness\n"
                "  [-t|--type=pudiverging-saturation]  Diverging map, varying saturation\n"
                "  [-t|--type=puqualitative-hue]       Qualitative map, evenly distributed hue\n"
                "  [-l|--lightness=L]                  Set lightness in [0,1]\n"
                "  [-L|--lightness-range=LR]           Set lightness range in [0.7,1]\n"
                "  [-s|--saturation=S]                 Set saturation in [0,1]\n"
                "  [-S|--saturation-range=SR]          Set saturation range in [0.7,1]\n"
                "  [-h|--hue=H]                        Set default hue in [0,360] degrees\n"
                "  [-d|--divergence=D]                 Set diverg. in deg for div. and qual. maps\n"
                "  [-r|--rotations=R]                  Set number of rotations for rainbow maps\n"
                "  [-T|--temperature=T]                Set start temp. in K for black body maps\n"
                "  [-R|--temperature-range=TR]         Set range for temperature in K\n"
                "  [-V|--hue-values=H0,H1,...]         Set hue values in [0,360] for multi-hue maps\n"
                "  [-P|--hue-positions=P0,P1,...]      Set hue positions in [0,1] for multi-hue maps\n"
                "CubeHelix color maps:\n"
                "  [-t|--type=cubehelix]               Generate a CubeHelix color map\n"
                "  [-h|--hue=H]                        Set start hue in [0,180] degrees\n"
                "  [-r|--rotations=R]                  Set number of rotations, in (-infty,infty)\n"
                "  [-s|--saturation=S]                 Set saturation, in [0,1]\n"
                "  [-g|--gamma=G]                      Set gamma correction, in (0,infty)\n"
                "Moreland diverging color maps:\n"
                "  [-t|--type=moreland]                Generate a Moreland diverging color map\n"
                "  [-A|--color0=sr,sg,sb]              Set the first color as sRGB in [0,255]\n"
                "  [-O|--color1=sr,sg,sb]              Set the last color as sRGB in [0,255]\n"
                "McNames sequential color maps:\n"
                "  [-t|--type=mcnames]                 Generate a McNames sequential color map\n"
                "  [-p|--periods=P]                    Set the number of periods in (0, infty)\n"
                "Defaults: format=csv, n=256, type=brewer-sequential\n"
                "https://marlam.de/gencolormap\n", argv[0]);
        return 0;
    }

    if (format < 0) {
        fprintf(stderr, "Invalid argument for option -f|--format.\n");
        return 1;
    }
    if (n < 2) {
        fprintf(stderr, "Invalid argument for option -n|--n.\n");
        return 1;
    }
    if (type < 0) {
        fprintf(stderr, "Invalid argument for option -t|--type.\n");
        return 1;
    }
    if (hue < 0.0f) {
        if (type == brewer_seq)
            hue = BrewerSequentialDefaultHue;
        else if (type == brewer_div)
            hue = BrewerDivergingDefaultHue;
        else if (type == brewer_qual)
            hue = BrewerQualitativeDefaultHue;
        else if (type == puseq_lightness)
            hue = PUSequentialLightnessDefaultHue;
        else if (type == puseq_saturation)
            hue = PUSequentialSaturationDefaultHue;
        else if (type == puseq_rainbow)
            hue = PUSequentialRainbowDefaultHue;
        else if (type == pudiv_lightness)
            hue = PUDivergingLightnessDefaultHue;
        else if (type == pudiv_saturation)
            hue = PUDivergingSaturationDefaultHue;
        else if (type == puqual_hue)
            hue = PUQualitativeHueDefaultHue;
        else if (type == cubehelix)
            hue = CubeHelixDefaultHue;
    }
    if (divergence < 0.0f) {
        if (type == brewer_div)
            divergence = BrewerDivergingDefaultDivergence;
        else if (type == brewer_qual)
            divergence = BrewerQualitativeDefaultDivergence;
        else if (type == pudiv_lightness)
            divergence = PUDivergingLightnessDefaultDivergence;
        else if (type == pudiv_saturation)
            divergence = PUDivergingSaturationDefaultDivergence;
        else if (type == puqual_hue)
            divergence = PUQualitativeHueDefaultDivergence;
    }
    if (contrast < 0.0f) {
        if (type == brewer_seq)
            contrast = (n <= 9 ? BrewerSequentialDefaultContrastForSmallN(n)
                    : BrewerSequentialDefaultContrast);
        else if (type == brewer_div)
            contrast = (n <= 9 ? BrewerDivergingDefaultContrastForSmallN(n)
                    : BrewerDivergingDefaultContrast);
        else if (type == brewer_qual)
            contrast = BrewerQualitativeDefaultContrast;
    }
    if (saturation < 0.0f) {
        if (type == brewer_seq)
            saturation = BrewerSequentialDefaultSaturation;
        else if (type == brewer_div)
            saturation = BrewerDivergingDefaultSaturation;
        else if (type == brewer_qual)
            saturation = BrewerQualitativeDefaultSaturation;
        else if (type == puseq_lightness)
            saturation = PUSequentialLightnessDefaultSaturation;
        else if (type == puseq_saturation)
            saturation = PUSequentialSaturationDefaultSaturation;
        else if (type == puseq_rainbow)
            saturation = PUSequentialRainbowDefaultSaturation;
        else if (type == puseq_blackbody)
            saturation = PUSequentialBlackBodyDefaultSaturation;
        else if (type == puseq_multihue)
            saturation = PUSequentialMultiHueDefaultSaturation;
        else if (type == pudiv_lightness)
            saturation = PUDivergingLightnessDefaultSaturation;
        else if (type == pudiv_saturation)
            saturation = PUDivergingSaturationDefaultSaturation;
        else if (type == puqual_hue)
            saturation = PUQualitativeHueDefaultSaturation;
        else if (type == cubehelix)
            saturation = CubeHelixDefaultSaturation;
    }
    if (saturation_range < 0.0f) {
        if (type == puseq_lightness)
            saturation_range = PUSequentialLightnessDefaultSaturationRange;
        else if (type == puseq_saturation)
            saturation_range = PUSequentialSaturationDefaultSaturationRange;
        else if (type == puseq_rainbow)
            saturation_range = PUSequentialRainbowDefaultSaturationRange;
        else if (type == puseq_blackbody)
            saturation_range = PUSequentialBlackBodyDefaultSaturationRange;
        else if (type == puseq_multihue)
            saturation_range = PUSequentialMultiHueDefaultSaturationRange;
        else if (type == pudiv_lightness)
            saturation_range = PUDivergingLightnessDefaultSaturationRange;
        else if (type == pudiv_saturation)
            saturation_range = PUDivergingSaturationDefaultSaturationRange;
    }
    if (brightness < 0.0f) {
        if (type == brewer_seq)
            brightness = BrewerSequentialDefaultBrightness;
        else if (type == brewer_div)
            brightness = BrewerDivergingDefaultBrightness;
        else if (type == brewer_qual)
            brightness = BrewerQualitativeDefaultBrightness;
    }
    if (warmth < 0.0f) {
        if (type == brewer_seq)
            warmth = BrewerSequentialDefaultWarmth;
        else if (type == brewer_div)
            warmth = BrewerDivergingDefaultWarmth;
    }
    if (lightness < 0.0f) {
        if (type == puseq_saturation)
            lightness = PUSequentialSaturationDefaultLightness;
        else if (type == pudiv_saturation)
            lightness = PUDivergingSaturationDefaultLightness;
        else if (type == puqual_hue)
            lightness = PUQualitativeHueDefaultLightness;
    }
    if (lightness_range < 0.0f) {
        if (type == puseq_lightness)
            lightness_range = PUSequentialLightnessDefaultLightnessRange;
        else if (type == puseq_rainbow)
            lightness_range = PUSequentialRainbowDefaultLightnessRange;
        else if (type == puseq_blackbody)
            lightness_range = PUSequentialBlackBodyDefaultLightnessRange;
        else if (type == puseq_multihue)
            lightness_range = PUSequentialMultiHueDefaultLightnessRange;
        else if (type == pudiv_lightness)
            lightness_range = PUDivergingLightnessDefaultLightnessRange;
    }
    if (std::isnan(rotations)) {
        if (type == puseq_rainbow)
            rotations = PUSequentialRainbowDefaultRotations;
        else if (type == cubehelix)
            rotations = CubeHelixDefaultRotations;
    }
    if (temperature < 0.0f) {
        if (type == puseq_blackbody)
            temperature = PUSequentialBlackBodyDefaultTemperature;
    }
    if (temperature_range < 0.0f) {
        if (type == puseq_blackbody)
            temperature_range = PUSequentialBlackBodyDefaultTemperatureRange;
    }
    if (gamma < 0.0f) {
        if (type == cubehelix)
            gamma = CubeHelixDefaultGamma;
    }
    if (hue_values.size() == 0) {
        if (type == puseq_multihue) {
            hue_values.resize(PUSequentialMultiHueDefaultHues);
            for (size_t i = 0; i < hue_values.size(); i++)
                hue_values[i] = PUSequentialMultiHueDefaultHueValues[i];
        }
    }
    if (hue_positions.size() == 0) {
        if (type == puseq_multihue) {
            hue_positions.resize(PUSequentialMultiHueDefaultHues);
            for (size_t i = 0; i < hue_positions.size(); i++)
                hue_positions[i] = PUSequentialMultiHueDefaultHuePositions[i];
        }
    }
    if (hue_values.size() != hue_positions.size()) {
        fprintf(stderr, "Number of hue values and positions do not match.\n");
        return 1;
    }
    if (!have_color0) {
        if (type == moreland) {
            color0[0] = MorelandDefaultR0;
            color0[1] = MorelandDefaultG0;
            color0[2] = MorelandDefaultB0;
        }
    }
    if (!have_color1) {
        if (type == moreland) {
            color1[0] = MorelandDefaultR1;
            color1[1] = MorelandDefaultG1;
            color1[2] = MorelandDefaultB1;
        }
    }
    if (std::isnan(periods)) {
        if (type == mcnames)
            periods = McNamesDefaultPeriods;
    }

    std::vector<unsigned char> colormap(3 * n);
    int clipped;
    switch (type) {
    case brewer_seq:
        clipped = BrewerSequential(n, &(colormap[0]), hue, contrast, saturation, brightness, warmth);
        break;
    case brewer_div:
        clipped = BrewerDiverging(n, &(colormap[0]), hue, divergence, contrast, saturation, brightness, warmth);
        break;
    case brewer_qual:
        clipped = BrewerQualitative(n, &(colormap[0]), hue, divergence, contrast, saturation, brightness);
        break;
    case puseq_lightness:
        clipped = PUSequentialLightness(n, &(colormap[0]), lightness_range, saturation_range, saturation, hue);
        break;
    case puseq_saturation:
        clipped = PUSequentialSaturation(n, &(colormap[0]), saturation_range, lightness, saturation, hue);
        break;
    case puseq_rainbow:
        clipped = PUSequentialRainbow(n, &(colormap[0]), lightness_range, saturation_range, hue, rotations, saturation);
        break;
    case puseq_blackbody:
        clipped = PUSequentialBlackBody(n, &(colormap[0]), temperature, temperature_range, lightness_range, saturation_range, saturation);
        break;
    case puseq_multihue:
        clipped = PUSequentialMultiHue(n, &(colormap[0]), lightness_range, saturation_range, saturation,
                hue_values.size(), hue_values.data(), hue_positions.data());
        break;
    case pudiv_lightness:
        clipped = PUDivergingLightness(n, &(colormap[0]), lightness_range, saturation_range, saturation, hue, divergence);
        break;
    case pudiv_saturation:
        clipped = PUDivergingSaturation(n, &(colormap[0]), saturation_range, lightness, saturation, hue, divergence);
        break;
    case puqual_hue:
        clipped = PUQualitativeHue(n, &(colormap[0]), hue, divergence, lightness, saturation);
        break;
    case cubehelix:
        clipped = CubeHelix(n, &(colormap[0]), hue, rotations, saturation, gamma);
        break;
    case moreland:
        clipped = Moreland(n, &(colormap[0]),
                color0[0], color0[1], color0[2],
                color1[0], color1[1], color1[2]);
        break;
    case mcnames:
        clipped = McNames(n, &(colormap[0]), periods);
        break;
    }


	// output 
    std::string output;
    if (format == csv) {
        output = ToCSV(n, colormap.data());
    } else if (format == json) {
        output = ToJSON(n, colormap.data());
    } else {
        output = ToPPM(n, colormap.data());
    }
    fputs(output.c_str(), stdout); // prints output
    fprintf(stderr, "%d color(s) were clipped\n", clipped);
    
    
	switch (format){
		case 0: fprintf(stdout, "output format = csv = %d\n", format); break;
		case 1: fprintf(stdout, "output format = json = %d\n", format); break;
		case 2: fprintf(stdout, "output format = ppm = %d\n", format); break;}
		
	fprintf(stdout, "All colors are represented as unsigned char sRGB triplets, with each value in [0,255]\n");
	
	fprintf(stdout, "number of colors in the map n = %d\n", n);
	
	switch (type){
			case 0: fprintf(stdout, "map type = brewer_seq = %d\n", type); break;
			case 1: fprintf(stdout, "map type = brewer_div = %d\n", type); break;
			case 2: fprintf(stdout, "map type = brewer_qual = %d\n", type); break;
			case 3: fprintf(stdout, "map type = puseq_lightness = %d\n", type); break;
			case 4: fprintf(stdout, "map type = puseq_saturation = %d\n", type); break;
			case 5: fprintf(stdout, "map type = puseq_rainbow = %d\n", type); break;
			case 6: fprintf(stdout, "map type = puseq_blackbody = %d\n", type); break;
			case 7: fprintf(stdout, "map type = puseq_multihue = %d\n", type); break;
			case 8: fprintf(stdout, "map type = pudiv_lightness = %d\n", type); break;
			case 9: fprintf(stdout, "map type = pudiv_saturation = %d\n", type); break;
			case 10: fprintf(stdout, "map type = puqual_hue = %d\n", type); break;
			case 11: fprintf(stdout, "map type = cubehelix = %d\n", type); break;
			case 12: fprintf(stdout, "map type = moreland = %d\n", type); break;
			case 13: fprintf(stdout, "map type = mcnames = %d\n", type); break;}
    
    
    
    
    

    return 0;
}
