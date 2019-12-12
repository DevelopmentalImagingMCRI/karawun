
# colour conversion stuff ported from:
# https://github.com/boxerab/dcmtk/blob/master/dcmiod/libsrc/cielabutil.cc


# Initialize white points of D65 light point (CIELab standard white point)
D65_WHITEPOINT_X = 0.950456
D65_WHITEPOINT_Y = 1.0
D65_WHITEPOINT_Z = 1.088754

nice_colours_rgb = [(0, 0, 0),
                    (0.99, 0.1, 0.1),
                    (0.1, 0.803921568627451, 0.1),
                    (0.1, 0.1, 0.99),
                    (0.1, 0.99, 0.99),
                    (0.99, 0.1, 0.99),
                    (0.99, 0.498039215686275, 0.1),
                    (0.1, 0.392156862745098, 0.1),
                    (0.541176470588235, 0.168627450980392, 0.886274509803922),
                    (0.545098039215686, 0.137254901960784, 0.137254901960784),
                    (0.1, 0.1, 0.501960784313725),
                    (0.545098039215686, 0.545098039215686, 0.1),
                    (0.99, 0.243137254901961, 0.588235294117647),
                    (0.545098039215686, 0.298039215686275, 0.223529411764706),
                    (0.1, 0.525490196078431, 0.545098039215686),
                    (0.803921568627451, 0.407843137254902, 0.223529411764706),
                    (0.749019607843137, 0.243137254901961, 0.99),
                    (0.1, 0.545098039215686, 0.270588235294118),
                    (0.780392156862745, 0.0823529411764706, 0.52156862745098),
                    (0.803921568627451, 0.215686274509804, 0.1),
                    (0.125490196078431, 0.698039215686274, 0.666666666666667),
                    (0.415686274509804, 0.352941176470588, 0.803921568627451),
                    (0.99, 0.0784313725490196, 0.576470588235294),
                    (0.270588235294118, 0.545098039215686, 0.454901960784314),
                    (0.282352941176471, 0.462745098039216, 0.99),
                    (0.803921568627451, 0.309803921568627, 0.223529411764706),
                    (0.1, 0.1, 0.803921568627451),
                    (0.545098039215686, 0.133333333333333, 0.32156862745098),
                    (0.545098039215686, 0.1, 0.545098039215686),
                    (0.933333333333333, 0.509803921568627, 0.933333333333333),
                    (0.545098039215686, 0.1, 0.1)]


def rgb2DicomLab(RGB):
    Lab = rgb2Lab(RGB)
    return lab2DicomLab(Lab)


def dicomlab2Lab(Labdicom):
    LDicom, aDicom, bDicom = Labdicom
    # results in 0 <= L <= 100
    L = ((LDicom * 100.0) / 65535.0)
    # results in -128 <= a <= 127
    a = ((aDicom * 255.0) / 65535.0) - 128
    # results in -128 <= b <= 127
    b = ((bDicom * 255.0) / 65535.0) - 128
    return (L, a, b)


def lab2DicomLab(Lab):
    L, a, b = Lab
    # results in 0 <= L <= 65535
    LDicom = L * 65535.0 / 100.0
    # results in 0 <= a <= 65535
    aDicom = (a + 128) * 65535.0 / 255.0
    # results in 0 <= b <= 65535
    bDicom = (b + 128) * 65535.0 / 255.0
    return (round(LDicom), round(aDicom), round(bDicom))


def rgb2Lab(RGB):
    XYZ = rgb2Xyz(RGB)
    Lab = xyz2Lab(XYZ)
    return Lab


def gammaCorrection(n):
    if ((n) <= 0.0031306684425005883):
        return 12.92 * (n)
    else:
        return (1.055*pow((n), 0.416666666666666667) - 0.055)


def invGammaCorrection(n):
    if ((n) <= 0.0404482362771076):
        return ((n) / 12.92)
    else:
        return (pow(((n) + 0.055)/1.055, 2.4))


def rgb2Xyz(RGB):
    R, G, B = RGB
    # No good if the inverse gamma correction is done...
    # R = invGammaCorrection(R)
    # G = invGammaCorrection(G)
    # B = invGammaCorrection(B)
    X = (0.4123955889674142161*R +
         0.3575834307637148171*G +
         0.1804926473817015735*B)
    Y = (0.2125862307855955516*R +
         0.7151703037034108499*G +
         0.07220049864333622685*B)
    Z = (0.01929721549174694484*R +
         0.1191838645808485318*G +
         0.9504971251315797660*B)
    return (X, Y, Z)


def xyz2Lab(XYZ):
    X, Y, Z = XYZ
    X /= D65_WHITEPOINT_X
    Y /= D65_WHITEPOINT_Y
    Z /= D65_WHITEPOINT_Z
    X = labf(X)
    Y = labf(Y)
    Z = labf(Z)
    L = 116*Y - 16
    a = 500*(X - Y)
    b = 200*(Y - Z)
    return (L, a, b)


def lab2Rgb(Lab):
    XYZ = lab2Xyz(Lab)
    RGB = xyz2Rgb(XYZ)
    return RGB


def lab2Xyz(Lab):
    L, a, b = Lab
    L = (L + 16)/116
    a = L + a/500
    b = L - b/200
    X = D65_WHITEPOINT_X * labfInv(a)
    Y = D65_WHITEPOINT_Y * labfInv(L)
    Z = D65_WHITEPOINT_Z * labfInv(b)
    return (X, Y, Z)


def xyz2Rgb(XYZ):
    X, Y, Z = XYZ

    R1 = 3.2406*X - 1.5372*Y - 0.4986*Z
    G1 = -0.9689*X + 1.8758*Y + 0.0415*Z
    B1 = 0.0557*X - 0.2040*Y + 1.0570*Z

    Min = min((R1, G1, B1))

    # Force nonnegative values so that gamma correction is well-defined.
    if (Min < 0):
        R1 -= Min
        G1 -= Min
        B1 -= Min

    # Transform from RGB to R'G'B'
    R = gammaCorrection(R1)
    G = gammaCorrection(G1)
    B = gammaCorrection(B1)
    return (R, G, B)


def labf(n):
    if (n >= 8.85645167903563082e-3):
        return pow(n, 0.333333333333333)
    else:
        return (841.0/108.0)*(n) + (4.0/29.0)


def labfInv(n):
    if ((n) >= 0.206896551724137931):
        return (n)*(n)*(n)
    else:
        return (108.0/841.0)*((n) - (4.0/29.0))
