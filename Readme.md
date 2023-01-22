# msdf_c 0.1
A pure C99 multi-channel signed distance field generator.  Handles MSDF bitmap
Generation from given ttf/otf font (outlines).
Single Header STB-style library

---

**Rendering output is not quite what you would get from sdfgen but it should work**

[This library is my take on improving msdf-c with some fixes and improving its API](https://github.com/solenum/msdf-c)

[Based on the C++ implementation by Viktor Chlumský.](https://github.com/Chlumsky/msdfgen)

~~~
// example fragment shader
// scale is render_width/glyph_width
// render_width being the width of each rendered glyph
float median(float r, float g, float b) {
    return max(min(r, g), min(max(r, g), b));
}

void main()
{
  vec3 sample = texture(u_texture, uv).rgb;
  float dist = scale * (median(sample.r, sample.g, sample.b) - 0.5);
  float o = clamp(dist + 0.5, 0.0, 1.0);
  color = vec4(vec3(1.0), o);
}
~~~


Current issues:

* ~~Rendering is not quite what you get from msdfgen~~ [X]
