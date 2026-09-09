#!/usr/bin/env python2.5
"""Random green colour generator, written by dbr, for
http://stackoverflow.com/questions/1586147/how-to-generate-random-greenish-colors
"""

def hsv_to_rgb(h, s, v):
    """Converts HSV value to RGB values
    Hue is in range 0-359 (degrees), value/saturation are in range 0-1 (float)

    Direct implementation of:
    http://en.wikipedia.org/wiki/HSL_and_HSV#Conversion_from_HSV_to_RGB
    """
    h, s, v = [float(x) for x in (h, s, v)]

    hi = (h / 60) % 6
    hi = int(round(hi))

    f = (h / 60) - (h / 60)
    p = v * (1 - s)
    q = v * (1 - f * s)
    t = v * (1 - (1 - f) * s)

    if hi == 0:
        return v, t, p
    elif hi == 1:
        return q, v, p
    elif hi == 2:
        return p, v, t
    elif hi == 3:
        return p, q, v
    elif hi == 4:
        return t, p, v
    elif hi == 5:
        return v, p, q

def test():
    """Check examples on..
    http://en.wikipedia.org/wiki/HSL_and_HSV#Examples
    ..work correctly
    """
    def verify(got, expected):
        if got != expected:
            raise AssertionError("Got %s, expected %s" % (got, expected))

    verify(hsv_to_rgb(0, 1, 1), (1, 0, 0))
    verify(hsv_to_rgb(120, 0.5, 1.0), (0.5, 1, 0.5))
    verify(hsv_to_rgb(240, 1, 0.5), (0, 0, 0.5))

def main(hi=90, hf=140, numbcolors = 5): # hi and hf are random green'ish hue from hue wheel
    """Generate numbcolors (5 default) random RGB colours, and create some simple coloured HTML span tags to verify them.
    """
    test() # Run simple test suite

    from random import randint, uniform

    colors = []

    for i in range(numbcolors):
        # Tweak these values to change colours/variance
        h = randint(hi, hf)
        s = uniform(0.2, 1)
        v = uniform(0.3, 1)

        r, g, b = hsv_to_rgb(h, s, v)

        # Convert to 0-1 range for HTML output
#        r, g, b = [x*255 for x in (r, g, b)]

        colors.append([])
        colors[i].append(r)
        colors[i].append(g)
        colors[i].append(b)

        #print "<span style='background:rgb(%i, %i, %i)'>&nbsp;&nbsp;</span>" % (r, g, b)
    return colors

if __name__ == '__main__':
    main()
