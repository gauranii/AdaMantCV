# -----------------------------------------------------------------------------
# File:    python/adamant/rmath.py
# Paper:   Gauran, Xue, Chen, Ombao, Yu (2022). Ridge Penalization in High
#          Dimensional Testing with Applications to Imaging Genetics.
#          Front. Neurosci. 16:836100. Methods (all figures and tables).
# Purpose: The upper tail of the gamma distribution as R computes it,
#          pgamma(q, shape, rate = rate, lower.tail = FALSE), ported
#          operation by operation so that it returns the same doubles as R,
#          bit for bit. This matters for the moment-matched gamma
#          approximation of a permutation distribution: when the variance of
#          the statistic is rounding noise, shape and rate are huge, and a
#          difference of one unit in the last place inside the algorithm
#          changes the p-value at its first digit.
# Inputs:  None.
# Outputs: None.
# Usage:   from adamant.rmath import pgamma_upper
#          pgamma_upper(3.0, 2.5, 1.2)   # pgamma(3, 2.5, rate = 1.2, lower.tail = FALSE)
# Source:  R 4.6.1: src/nmath/pgamma.c (pgamma, pgamma_raw, pgamma_smallx,
#          pd_upper_series, pd_lower_cf, pd_lower_series, dpois_wrap,
#          ppois_asymp, dpnorm, logcf, log1pmx, lgamma1p), dpois.c
#          (dpois_raw), bd0.c (ebd0), stirlerr.c, pnorm.c (pnorm5,
#          pnorm_both), dnorm.c (dnorm4), lgamma.c (lgammafn), gamma.c
#          (gammafn), lgammacor.c, chebyshev.c (chebyshev_eval), fmax2.c;
#          src/library/stats/R/distn.R (scale = 1/rate).
# Notes:   R's C code is compiled for arm64 with a product that feeds a sum
#          or a difference in the same expression done in one rounding (fused
#          multiply-add). Every such place is written here as muladd(); the
#          other products and sums are rounded one at a time, as in R. Only
#          the upper tail is ported (lower.tail = FALSE), with log.p = TRUE
#          kept inside because pgamma_raw() uses it for very small results.
# -----------------------------------------------------------------------------
import math
import platform

import numpy as np
from llvmlite import ir
from numba import njit
from numba.core import cgutils, types
from numba.extending import intrinsic

# R is compiled with fused multiply-add on arm64 only.
FUSED = platform.machine().lower() in ("arm64", "aarch64")

if FUSED:

    @intrinsic
    def muladd(typingctx, a, b, c):
        """a * b + c with a single rounding (fused multiply-add).

        Parameters: a, b, c (floats). Returns: float.
        """
        sig = types.float64(types.float64, types.float64, types.float64)

        def codegen(context, builder, signature, args):
            double = ir.DoubleType()
            fnty = ir.FunctionType(double, [double, double, double])
            fn = builder.module.declare_intrinsic("llvm.fma", [double], fnty)
            return builder.call(fn, args)

        return sig, codegen

else:

    @njit(cache=True, error_model="numpy")
    def muladd(a, b, c):
        """a * b + c with two roundings, for machines where R does not fuse.

        Parameters: a, b, c (floats). Returns: float.
        """
        return a * b + c


@intrinsic
def _libm_lgamma(typingctx, x):
    """lgamma() of the system C library, the one R's stirlerr() calls.

    Parameters: x (float). Returns: float.
    """
    sig = types.float64(types.float64)

    def codegen(context, builder, signature, args):
        double = ir.DoubleType()
        fnty = ir.FunctionType(double, [double])
        fn = cgutils.get_or_insert_function(builder.module, fnty, "lgamma")
        return builder.call(fn, args)

    return sig, codegen


DBL_EPSILON = 2.220446049250313e-16
DBL_MIN = 2.2250738585072014e-308
DBL_MAX = 1.7976931348623157e308
M_LN2 = 0.693147180559945309417232121458
M_2PI = 6.283185307179586476925286766559
M_SQRT_32 = 5.656854249492380195206754896838
M_1_SQRT_2PI = 0.398942280401432677939946059934
M_LN_2PI = 1.837877066409345483560659472811
M_LN_SQRT_2PI = 0.918938533204672741780329736406
M_SQRT_2PI = 2.50662827463100050241576528481104525301
X_LRG = 2.86111748575702815380240589208115399625e307  # 2^1023 / pi

SCALEFACTOR = 4294967296.0**8  # 2^256
M_CUTOFF = M_LN2 * 1024 / DBL_EPSILON
EULERS_CONST = 0.5772156649015328606065120900824024
LGAMMA1P_C = 0.2273736845824652515226821577978691e-12  # zeta(42) - 1
TOL_LOGCF = 1e-14
MIN_LOG1_VALUE = -0.79149064
DNORM_ZERO = math.sqrt(-2 * M_LN2 * (-1021 + 1 - 53))  # beyond it dnorm is 0
DNORM_HUGE = 2 * math.sqrt(DBL_MAX)
LGAMMA_XMAX = 2.5327372760800758e305
GAMMA_XSML = 2.2474362225598545e-308
LGAMMACOR_XBIG = 94906265.62425156

# Branches of pgamma_raw(), as returned by pgamma_branch().
BRANCHES = (
    "outside (x <= 0, x = Inf, shape = 0, NaN or invalid)",
    "pgamma_smallx (x < 1)",
    "pd_upper_series (x <= shape - 1)",
    "pd_lower_cf (shape < 1 <= x)",
    "pd_lower_series (x > shape - 1, shape >= 1)",
    "ppois_asymp (x near shape, both large)",
    "one (shape < 1, x * eps > 1 - shape)",
)

_LGAMMA1P_COEFFS = np.array(
    [
        0.3224670334241132182362075833230126e-0,
        0.6735230105319809513324605383715000e-1,
        0.2058080842778454787900092413529198e-1,
        0.7385551028673985266273097291406834e-2,
        0.2890510330741523285752988298486755e-2,
        0.1192753911703260977113935692828109e-2,
        0.5096695247430424223356548135815582e-3,
        0.2231547584535793797614188036013401e-3,
        0.9945751278180853371459589003190170e-4,
        0.4492623673813314170020750240635786e-4,
        0.2050721277567069155316650397830591e-4,
        0.9439488275268395903987425104415055e-5,
        0.4374866789907487804181793223952411e-5,
        0.2039215753801366236781900709670839e-5,
        0.9551412130407419832857179772951265e-6,
        0.4492469198764566043294290331193655e-6,
        0.2120718480555466586923135901077628e-6,
        0.1004322482396809960872083050053344e-6,
        0.4769810169363980565760193417246730e-7,
        0.2271109460894316491031998116062124e-7,
        0.1083865921489695409107491757968159e-7,
        0.5183475041970046655121248647057669e-8,
        0.2483674543802478317185008663991718e-8,
        0.1192140140586091207442548202774640e-8,
        0.5731367241678862013330194857961011e-9,
        0.2759522885124233145178149692816341e-9,
        0.1330476437424448948149715720858008e-9,
        0.6422964563838100022082448087644648e-10,
        0.3104424774732227276239215783404066e-10,
        0.1502138408075414217093301048780668e-10,
        0.7275974480239079662504549924814047e-11,
        0.3527742476575915083615072228655483e-11,
        0.1711991790559617908601084114443031e-11,
        0.8315385841420284819798357793954418e-12,
        0.4042200525289440065536008957032895e-12,
        0.1966475631096616490411045679010286e-12,
        0.9573630387838555763782200936508615e-13,
        0.4664076026428374224576492565974577e-13,
        0.2273736960065972320633279596737272e-13,
        0.1109139947083452201658320007192334e-13,
    ]
)
_SFERR_HALVES = np.array(
    [
        0.0,
        0.1534264097200273452913848,
        0.0810614667953272582196702,
        0.0548141210519176538961390,
        0.0413406959554092940938221,
        0.03316287351993628748511048,
        0.02767792568499833914878929,
        0.02374616365629749597132920,
        0.02079067210376509311152277,
        0.01848845053267318523077934,
        0.01664469118982119216319487,
        0.01513497322191737887351255,
        0.01387612882307074799874573,
        0.01281046524292022692424986,
        0.01189670994589177009505572,
        0.01110455975820691732662991,
        0.010411265261972096497478567,
        0.009799416126158803298389475,
        0.009255462182712732917728637,
        0.008768700134139385462952823,
        0.008330563433362871256469318,
        0.007934114564314020547248100,
        0.007573675487951840794972024,
        0.007244554301320383179543912,
        0.006942840107209529865664152,
        0.006665247032707682442354394,
        0.006408994188004207068439631,
        0.006171712263039457647532867,
        0.005951370112758847735624416,
        0.005746216513010115682023589,
        0.005554733551962801371038690,
    ]
)
_S = np.array(
    [
        0.083333333333333333333,
        0.00277777777777777777778,
        0.00079365079365079365079365,
        0.000595238095238095238095238,
        0.0008417508417508417508417508,
        0.0019175269175269175269175262,
        0.0064102564102564102564102561,
        0.029550653594771241830065352,
        0.17964437236883057316493850,
        1.3924322169059011164274315,
        13.402864044168391994478957,
        156.84828462600201730636509,
        2193.1033333333333333333333,
        36108.771253724989357173269,
        691472.26885131306710839498,
        15238221.539407416192283370,
        382900751.39141414141414141,
    ]
)
_PNORM_A = np.array(
    [
        2.2352520354606839287,
        161.02823106855587881,
        1067.6894854603709582,
        18154.981253343561249,
        0.065682337918207449113,
    ]
)
_PNORM_B = np.array(
    [
        47.20258190468824187,
        976.09855173777669322,
        10260.932208618978205,
        45507.789335026729956,
    ]
)
_PNORM_C = np.array(
    [
        0.39894151208813466764,
        8.8831497943883759412,
        93.506656132177855979,
        597.27027639480026226,
        2494.5375852903726711,
        6848.1904505362823326,
        11602.651437647350124,
        9842.7148383839780218,
        1.0765576773720192317e-8,
    ]
)
_PNORM_D = np.array(
    [
        22.266688044328115691,
        235.38790178262499861,
        1519.377599407554805,
        6485.558298266760755,
        18615.571640885098091,
        34900.952721145977266,
        38912.003286093271411,
        19685.429676859990727,
    ]
)
_PNORM_P = np.array(
    [
        0.21589853405795699,
        0.1274011611602473639,
        0.022235277870649807,
        0.001421619193227893466,
        2.9112874951168792e-5,
        0.02307344176494017303,
    ]
)
_PNORM_Q = np.array(
    [
        1.28426009614491121,
        0.468238212480865118,
        0.0659881378689285515,
        0.00378239633202758244,
        7.29751555083966205e-5,
    ]
)
_GAMCS = np.array(
    [
        +.8571195590989331421920062399942e-2,
        +.4415381324841006757191315771652e-2,
        +.5685043681599363378632664588789e-1,
        -.4219835396418560501012500186624e-2,
        +.1326808181212460220584006796352e-2,
        -.1893024529798880432523947023886e-3,
        +.3606925327441245256578082217225e-4,
        -.6056761904460864218485548290365e-5,
        +.1055829546302283344731823509093e-5,
        -.1811967365542384048291855891166e-6,
        +.3117724964715322277790254593169e-7,
        -.5354219639019687140874081024347e-8,
        +.9193275519859588946887786825940e-9,
        -.1577941280288339761767423273953e-9,
        +.2707980622934954543266540433089e-10,
        -.4646818653825730144081661058933e-11,
        +.7973350192007419656460767175359e-12,
        -.1368078209830916025799499172309e-12,
        +.2347319486563800657233471771688e-13,
        -.4027432614949066932766570534699e-14,
        +.6910051747372100912138336975257e-15,
        -.1185584500221992907052387126192e-15,
    ]
)
_ALGMCS = np.array(
    [
        +.1666389480451863247205729650822e+0,
        -.1384948176067563840732986059135e-4,
        +.9810825646924729426157171547487e-8,
        -.1809129475572494194263306266719e-10,
        +.6221098041892605227126015543416e-13,
    ]
)
_BD0_SCALE_HEX = (
    "+0x1.62e430p-1 -0x1.05c610p-29 -0x1.950d88p-54 +0x1.d9cc02p-79",
    "+0x1.5ee02cp-1 -0x1.6dbe98p-25 -0x1.51e540p-50 +0x1.2bfa48p-74",
    "+0x1.5ad404p-1 +0x1.86b3e4p-26 +0x1.9f6534p-50 +0x1.54be04p-74",
    "+0x1.570124p-1 -0x1.9ed750p-25 -0x1.f37dd0p-51 +0x1.10b770p-77",
    "+0x1.5326e4p-1 -0x1.9b9874p-25 -0x1.378194p-49 +0x1.56feb2p-74",
    "+0x1.4f4528p-1 +0x1.aca70cp-28 +0x1.103e74p-53 +0x1.9c410ap-81",
    "+0x1.4b5bd8p-1 -0x1.6a91d8p-25 -0x1.8e43d0p-50 -0x1.afba9ep-77",
    "+0x1.47ae54p-1 -0x1.abb51cp-25 +0x1.19b798p-51 +0x1.45e09cp-76",
    "+0x1.43fa00p-1 -0x1.d06318p-25 -0x1.8858d8p-49 -0x1.1927c4p-75",
    "+0x1.3ffa40p-1 +0x1.1a427cp-25 +0x1.151640p-53 -0x1.4f5606p-77",
    "+0x1.3c7c80p-1 -0x1.19bf48p-34 +0x1.05fc94p-58 -0x1.c096fcp-82",
    "+0x1.38b320p-1 +0x1.6b5778p-25 +0x1.be38d0p-50 -0x1.075e96p-74",
    "+0x1.34e288p-1 +0x1.d9ce1cp-25 +0x1.316eb8p-49 +0x1.2d885cp-73",
    "+0x1.315124p-1 +0x1.c2fc60p-29 -0x1.4396fcp-53 +0x1.acf376p-78",
    "+0x1.2db954p-1 +0x1.720de4p-25 -0x1.d39b04p-49 -0x1.f11176p-76",
    "+0x1.2a1b08p-1 -0x1.562494p-25 +0x1.a7863cp-49 +0x1.85dd64p-73",
    "+0x1.267620p-1 +0x1.3430e0p-29 -0x1.96a958p-56 +0x1.f8e636p-82",
    "+0x1.23130cp-1 +0x1.7bebf4p-25 +0x1.416f1cp-52 -0x1.78dd36p-77",
    "+0x1.1faa34p-1 +0x1.70e128p-26 +0x1.81817cp-50 -0x1.c2179cp-76",
    "+0x1.1bf204p-1 +0x1.3a9620p-28 +0x1.2f94c0p-52 +0x1.9096c0p-76",
    "+0x1.187ce4p-1 -0x1.077870p-27 +0x1.655a80p-51 +0x1.eaafd6p-78",
    "+0x1.1501c0p-1 -0x1.406cacp-25 -0x1.e72290p-49 +0x1.5dd800p-73",
    "+0x1.11cb80p-1 +0x1.787cd0p-25 -0x1.efdc78p-51 -0x1.5380cep-77",
    "+0x1.0e4498p-1 +0x1.747324p-27 -0x1.024548p-51 +0x1.77a5a6p-75",
    "+0x1.0b036cp-1 +0x1.690c74p-25 +0x1.5d0cc4p-50 -0x1.c0e23cp-76",
    "+0x1.077070p-1 -0x1.a769bcp-27 +0x1.452234p-52 +0x1.6ba668p-76",
    "+0x1.04240cp-1 -0x1.a686acp-27 -0x1.ef46b0p-52 -0x1.5ce10cp-76",
    "+0x1.00d22cp-1 +0x1.fc0e10p-25 +0x1.6ee034p-50 -0x1.19a2ccp-74",
    "+0x1.faf588p-2 +0x1.ef1e64p-27 -0x1.26504cp-54 -0x1.b15792p-82",
    "+0x1.f4d87cp-2 +0x1.d7b980p-26 -0x1.a114d8p-50 +0x1.9758c6p-75",
    "+0x1.ee1414p-2 +0x1.2ec060p-26 +0x1.dc00fcp-52 +0x1.f8833cp-76",
    "+0x1.e7e32cp-2 -0x1.ac796cp-27 -0x1.a68818p-54 +0x1.235d02p-78",
    "+0x1.e108a0p-2 -0x1.768ba4p-28 -0x1.f050a8p-52 +0x1.00d632p-82",
    "+0x1.dac354p-2 -0x1.d3a6acp-30 +0x1.18734cp-57 -0x1.f97902p-83",
    "+0x1.d47424p-2 +0x1.7dbbacp-31 -0x1.d5ada4p-56 +0x1.56fcaap-81",
    "+0x1.ce1af0p-2 +0x1.70be7cp-27 +0x1.6f6fa4p-51 +0x1.7955a2p-75",
    "+0x1.c7b798p-2 +0x1.ec36ecp-26 -0x1.07e294p-50 -0x1.ca183cp-75",
    "+0x1.c1ef04p-2 +0x1.c1dfd4p-26 +0x1.888eecp-50 -0x1.fd6b86p-75",
    "+0x1.bb7810p-2 +0x1.478bfcp-26 +0x1.245b8cp-50 +0x1.ea9d52p-74",
    "+0x1.b59da0p-2 -0x1.882b08p-27 +0x1.31573cp-53 -0x1.8c249ap-77",
    "+0x1.af1294p-2 -0x1.b710f4p-27 +0x1.622670p-51 +0x1.128578p-76",
    "+0x1.a925d4p-2 -0x1.0ae750p-27 +0x1.574ed4p-51 +0x1.084996p-75",
    "+0x1.a33040p-2 +0x1.027d30p-29 +0x1.b9a550p-53 -0x1.b2e38ap-78",
    "+0x1.9d31c0p-2 -0x1.5ec12cp-26 -0x1.5245e0p-52 +0x1.2522d0p-79",
    "+0x1.972a34p-2 +0x1.135158p-30 +0x1.a5c09cp-56 +0x1.24b70ep-80",
    "+0x1.911984p-2 +0x1.0995d4p-26 +0x1.3bfb5cp-50 +0x1.2c9dd6p-75",
    "+0x1.8bad98p-2 -0x1.1d6144p-29 +0x1.5b9208p-53 +0x1.1ec158p-77",
    "+0x1.858b58p-2 -0x1.1b4678p-27 +0x1.56cab4p-53 -0x1.2fdc0cp-78",
    "+0x1.7f5fa0p-2 +0x1.3aaf48p-27 +0x1.461964p-51 +0x1.4ae476p-75",
    "+0x1.79db68p-2 -0x1.7e5054p-26 +0x1.673750p-51 -0x1.a11f7ap-76",
    "+0x1.744f88p-2 -0x1.cc0e18p-26 -0x1.1e9d18p-50 -0x1.6c06bcp-78",
    "+0x1.6e08ecp-2 -0x1.5d45e0p-26 -0x1.c73ec8p-50 +0x1.318d72p-74",
    "+0x1.686c80p-2 +0x1.e9b14cp-26 -0x1.13bbd4p-50 -0x1.efeb1cp-78",
    "+0x1.62c830p-2 -0x1.a8c70cp-27 -0x1.5a1214p-51 -0x1.bab3fcp-79",
    "+0x1.5d1bdcp-2 -0x1.4fec6cp-31 +0x1.423638p-56 +0x1.ee3feep-83",
    "+0x1.576770p-2 +0x1.7455a8p-26 -0x1.3ab654p-50 -0x1.26be4cp-75",
    "+0x1.5262e0p-2 -0x1.146778p-26 -0x1.b9f708p-52 -0x1.294018p-77",
    "+0x1.4c9f08p-2 +0x1.e152c4p-26 -0x1.dde710p-53 +0x1.fd2208p-77",
    "+0x1.46d2d8p-2 +0x1.c28058p-26 -0x1.936284p-50 +0x1.9fdd68p-74",
    "+0x1.41b940p-2 +0x1.cce0c0p-26 -0x1.1a4050p-50 +0x1.bc0376p-76",
    "+0x1.3bdd24p-2 +0x1.d6296cp-27 +0x1.425b48p-51 -0x1.cddb2cp-77",
    "+0x1.36b578p-2 -0x1.287ddcp-27 -0x1.2d0f4cp-51 +0x1.38447ep-75",
    "+0x1.31871cp-2 +0x1.2a8830p-27 +0x1.3eae54p-52 -0x1.898136p-77",
    "+0x1.2b9304p-2 -0x1.51d8b8p-28 +0x1.27694cp-52 -0x1.fd852ap-76",
    "+0x1.265620p-2 -0x1.d98f3cp-27 +0x1.a44338p-51 -0x1.56e85ep-78",
    "+0x1.211254p-2 +0x1.986160p-26 +0x1.73c5d0p-51 +0x1.4a861ep-75",
    "+0x1.1bc794p-2 +0x1.fa3918p-27 +0x1.879c5cp-51 +0x1.16107cp-78",
    "+0x1.1675ccp-2 -0x1.4545a0p-26 +0x1.c07398p-51 +0x1.f55c42p-76",
    "+0x1.111ce4p-2 +0x1.f72670p-37 -0x1.b84b5cp-61 +0x1.a4a4dcp-85",
    "+0x1.0c81d4p-2 +0x1.0c150cp-27 +0x1.218600p-51 -0x1.d17312p-76",
    "+0x1.071b84p-2 +0x1.fcd590p-26 +0x1.a3a2e0p-51 +0x1.fe5ef8p-76",
    "+0x1.01ade4p-2 -0x1.bb1844p-28 +0x1.db3cccp-52 +0x1.1f56fcp-77",
    "+0x1.fa01c4p-3 -0x1.12a0d0p-29 -0x1.f71fb0p-54 +0x1.e287a4p-78",
    "+0x1.ef0adcp-3 +0x1.7b8b28p-28 -0x1.35bce4p-52 -0x1.abc8f8p-79",
    "+0x1.e598ecp-3 +0x1.5a87e4p-27 -0x1.134bd0p-51 +0x1.c2cebep-76",
    "+0x1.da85d8p-3 -0x1.df31b0p-27 +0x1.94c16cp-57 +0x1.8fd7eap-82",
    "+0x1.d0fb80p-3 -0x1.bb5434p-28 -0x1.ea5640p-52 -0x1.8ceca4p-77",
    "+0x1.c765b8p-3 +0x1.e4d68cp-27 +0x1.5b59b4p-51 +0x1.76f6c4p-76",
    "+0x1.bdc46cp-3 -0x1.1cbb50p-27 +0x1.2da010p-51 +0x1.eb282cp-75",
    "+0x1.b27980p-3 -0x1.1b9ce0p-27 +0x1.7756f8p-52 +0x1.2ff572p-76",
    "+0x1.a8bed0p-3 -0x1.bbe874p-30 +0x1.85cf20p-56 +0x1.b9cf18p-80",
    "+0x1.9ef83cp-3 +0x1.2769a4p-27 -0x1.85bda0p-52 +0x1.8c8018p-79",
    "+0x1.9525a8p-3 +0x1.cf456cp-27 -0x1.7137d8p-52 -0x1.f158e8p-76",
    "+0x1.8b46f8p-3 +0x1.11b12cp-30 +0x1.9f2104p-54 -0x1.22836ep-78",
    "+0x1.83040cp-3 +0x1.2379e4p-28 +0x1.b71c70p-52 -0x1.990cdep-76",
    "+0x1.790ed4p-3 +0x1.dc4c68p-28 -0x1.910ac8p-52 +0x1.dd1bd6p-76",
    "+0x1.6f0d28p-3 +0x1.5cad68p-28 +0x1.737c94p-52 -0x1.9184bap-77",
    "+0x1.64fee8p-3 +0x1.04bf88p-28 +0x1.6fca28p-52 +0x1.8884a8p-76",
    "+0x1.5c9400p-3 +0x1.d65cb0p-29 -0x1.b2919cp-53 +0x1.b99bcep-77",
    "+0x1.526e60p-3 -0x1.c5e4bcp-27 -0x1.0ba380p-52 +0x1.d6e3ccp-79",
    "+0x1.483bccp-3 +0x1.9cdc7cp-28 -0x1.5ad8dcp-54 -0x1.392d3cp-83",
    "+0x1.3fb25cp-3 -0x1.a6ad74p-27 +0x1.5be6b4p-52 -0x1.4e0114p-77",
    "+0x1.371fc4p-3 -0x1.fe1708p-27 -0x1.78864cp-52 -0x1.27543ap-76",
    "+0x1.2cca10p-3 -0x1.4141b4p-28 -0x1.ef191cp-52 +0x1.00ee08p-76",
    "+0x1.242310p-3 +0x1.3ba510p-27 -0x1.d003c8p-51 +0x1.162640p-76",
    "+0x1.1b72acp-3 +0x1.52f67cp-27 -0x1.fd6fa0p-51 +0x1.1a3966p-77",
    "+0x1.10f8e4p-3 +0x1.129cd8p-30 +0x1.31ef30p-55 +0x1.a73e38p-79",
    "+0x1.08338cp-3 -0x1.005d7cp-27 -0x1.661a9cp-51 +0x1.1f138ap-79",
    "+0x1.fec914p-4 -0x1.c482a8p-29 -0x1.55746cp-54 +0x1.99f932p-80",
    "+0x1.ed1794p-4 +0x1.d06f00p-29 +0x1.75e45cp-53 -0x1.d0483ep-78",
    "+0x1.db5270p-4 +0x1.87d928p-32 -0x1.0f52a4p-57 +0x1.81f4a6p-84",
    "+0x1.c97978p-4 +0x1.af1d24p-29 -0x1.0977d0p-60 -0x1.8839d0p-84",
    "+0x1.b78c84p-4 -0x1.44f124p-28 -0x1.ef7bc4p-52 +0x1.9e0650p-78",
    "+0x1.a58b60p-4 +0x1.856464p-29 +0x1.c651d0p-55 +0x1.b06b0cp-79",
    "+0x1.9375e4p-4 +0x1.5595ecp-28 +0x1.dc3738p-52 +0x1.86c89ap-81",
    "+0x1.814be4p-4 -0x1.c073fcp-28 -0x1.371f88p-53 -0x1.5f4080p-77",
    "+0x1.6f0d28p-4 +0x1.5cad68p-29 +0x1.737c94p-53 -0x1.9184bap-78",
    "+0x1.60658cp-4 -0x1.6c8af4p-28 +0x1.d8ef74p-55 +0x1.c4f792p-80",
    "+0x1.4e0110p-4 +0x1.146b5cp-29 +0x1.73f7ccp-54 -0x1.d28db8p-79",
    "+0x1.3b8758p-4 +0x1.8b1b70p-28 -0x1.20aca4p-52 -0x1.651894p-76",
    "+0x1.28f834p-4 +0x1.43b6a4p-30 -0x1.452af8p-55 +0x1.976892p-80",
    "+0x1.1a0fbcp-4 -0x1.e4075cp-28 +0x1.1fe618p-52 +0x1.9d6dc2p-77",
    "+0x1.075984p-4 -0x1.4ce370p-29 -0x1.d9fc98p-53 +0x1.4ccf12p-77",
    "+0x1.f0a30cp-5 +0x1.162a68p-37 -0x1.e83368p-61 -0x1.d222a6p-86",
    "+0x1.cae730p-5 -0x1.1a8f7cp-31 -0x1.5f9014p-55 +0x1.2720c0p-79",
    "+0x1.ac9724p-5 -0x1.e8ee08p-29 +0x1.a7de04p-54 -0x1.9bba74p-78",
    "+0x1.868a84p-5 -0x1.ef8128p-30 +0x1.dc5eccp-54 -0x1.58d250p-79",
    "+0x1.67f950p-5 -0x1.ed684cp-30 -0x1.f060c0p-55 -0x1.b1294cp-80",
    "+0x1.494accp-5 +0x1.a6c890p-32 -0x1.c3ad48p-56 -0x1.6dc66cp-84",
    "+0x1.22c71cp-5 -0x1.8abe2cp-32 -0x1.7e7078p-56 -0x1.ddc3dcp-86",
    "+0x1.03d5d8p-5 +0x1.79cfbcp-31 -0x1.da7c4cp-58 +0x1.4e7582p-83",
    "+0x1.c98d18p-6 +0x1.a01904p-31 -0x1.854164p-55 +0x1.883c36p-79",
    "+0x1.8b31fcp-6 -0x1.356500p-30 +0x1.c3ab48p-55 +0x1.b69bdap-80",
    "+0x1.3cea44p-6 +0x1.a352bcp-33 -0x1.8865acp-57 -0x1.48159cp-81",
    "+0x1.fc0a8cp-7 -0x1.e07f84p-32 +0x1.e7cf6cp-58 +0x1.3a69c0p-82",
    "+0x1.7dc474p-7 +0x1.f810a8p-31 -0x1.245b5cp-56 -0x1.a1f4f8p-80",
    "+0x1.fe02a8p-8 -0x1.4ef988p-32 +0x1.1f86ecp-57 +0x1.20723cp-81",
    "+0x1.ff00acp-9 -0x1.d4ef44p-33 +0x1.2821acp-63 +0x1.5a6d32p-87",
    "0 0 0 0",
)

# The table of logarithms of ebd0(): 129 rows of four single-precision parts.
_BD0_SCALE = np.array(
    [[float.fromhex(part) for part in row.split()] for row in _BD0_SCALE_HEX]
)

_PPOIS_A = np.array(
    [
        -1e99,
        2 / 3.0,
        -4 / 135.0,
        8 / 2835.0,
        16 / 8505.0,
        -8992 / 12629925.0,
        -334144 / 492567075.0,
        698752 / 1477701225.0,
    ]
)
_PPOIS_B = np.array(
    [
        -1e99,
        1 / 12.0,
        1 / 288.0,
        -139 / 51840.0,
        -571 / 2488320.0,
        163879 / 209018880.0,
        5246819 / 75246796800.0,
        -534703531 / 902961561600.0,
    ]
)


@njit(cache=True, error_model="numpy")
def _log1_exp(x):
    """log(1 - exp(x)) for x < 0: the macro R_Log1_Exp.

    Parameters: x (float). Returns: float.
    """
    if x > -M_LN2:
        return math.log(-math.expm1(x))
    return math.log1p(-math.exp(x))


@njit(cache=True, error_model="numpy")
def _fmax2(x, y):
    """The larger of two numbers, NaN if either is: fmax2().

    Parameters: x, y (floats). Returns: float.
    """
    if math.isnan(x) or math.isnan(y):
        return x + y
    return y if x < y else x


@njit(cache=True, error_model="numpy")
def _logcf(x, i, d, eps):
    """Continued fraction for sum_k x^k / (i + k d): logcf().

    Parameters: x, i (> 0), d (>= 0), eps (relative tolerance).
    Returns: float.
    """
    c1 = 2 * d
    c2 = i + d
    c4 = c2 + d
    a1 = c2
    b1 = i * muladd(-i, x, c2)
    b2 = d * d * x
    a2 = muladd(c4, c2, -b2)
    b2 = muladd(c4, b1, -(i * b2))
    while abs(muladd(a2, b1, -(a1 * b2))) > abs(eps * b1 * b2):
        c3 = c2 * c2 * x
        c2 += d
        c4 += d
        a1 = muladd(c4, a2, -(c3 * a1))
        b1 = muladd(c4, b2, -(c3 * b1))

        c3 = c1 * c1 * x
        c1 += d
        c4 += d
        a2 = muladd(c4, a1, -(c3 * a2))
        b2 = muladd(c4, b1, -(c3 * b2))

        if abs(b2) > SCALEFACTOR:
            a1 /= SCALEFACTOR
            b1 /= SCALEFACTOR
            a2 /= SCALEFACTOR
            b2 /= SCALEFACTOR
        elif abs(b2) < 1 / SCALEFACTOR:
            a1 *= SCALEFACTOR
            b1 *= SCALEFACTOR
            a2 *= SCALEFACTOR
            b2 *= SCALEFACTOR
    return a2 / b2


@njit(cache=True, error_model="numpy")
def _log1pmx(x):
    """log(1 + x) - x, accurate for small x: log1pmx().

    Parameters: x (float, > -1). Returns: float.
    """
    if x > 1 or x < MIN_LOG1_VALUE:
        return math.log1p(x) - x
    r = x / (2 + x)
    y = r * r
    if abs(x) < 1e-2:
        t = muladd(2.0 / 9, y, 2.0 / 7)
        t = muladd(t, y, 2.0 / 5)
        t = muladd(t, y, 2.0 / 3)
        return r * muladd(t, y, -x)
    return r * muladd(2 * y, _logcf(y, 3.0, 2.0, TOL_LOGCF), -x)


@njit(cache=True, error_model="numpy")
def _chebyshev_eval(x, a, n):
    """Chebyshev series with the first n coefficients of a: chebyshev_eval().

    Parameters: x (in [-1.1, 1.1]), a (array), n (number of terms).
    Returns: float.
    """
    if x < -1.1 or x > 1.1:
        return np.nan
    twox = x * 2
    b2 = 0.0
    b1 = 0.0
    b0 = 0.0
    for i in range(1, n + 1):
        b2 = b1
        b1 = b0
        b0 = muladd(twox, b1, -b2) + a[n - i]
    return (b0 - b2) * 0.5


@njit(cache=True, error_model="numpy")
def _gammafn_small(x):
    """Gamma(x) for 0 < x <= 10: that range of gammafn().

    Parameters: x (float in (0, 10]). Returns: float.
    """
    n = int(x)
    y = x - n  # in [0, 1)
    n -= 1
    value = _chebyshev_eval(muladd(y, 2.0, -1.0), _GAMCS, 22) + 0.9375
    if n == 0:
        return value
    if n < 0:  # 0 < x < 1
        if y < GAMMA_XSML:
            return np.inf
        return value / x
    for i in range(1, n + 1):
        value *= y + i
    return value


@njit(cache=True, error_model="numpy")
def _lgammacor(x):
    """The correction of Stirling's formula for x >= 10: lgammacor().

    Parameters: x (float, >= 10). Returns: float.
    """
    if x < 10:
        return np.nan
    if x < LGAMMACOR_XBIG:
        tmp = 10 / x
        return _chebyshev_eval(muladd(tmp * tmp, 2.0, -1.0), _ALGMCS, 5) / x
    return 1 / (x * 12)


@njit(cache=True, error_model="numpy")
def _lgammafn(x):
    """log(Gamma(x)) for x >= 0 as R computes it: lgammafn().

    Parameters: x (float, >= 0; this port is not for negative x).
    Returns: float.
    """
    if math.isnan(x):
        return x
    if x <= 0:
        return np.inf
    if x < 1e-306:
        return -math.log(x)
    if x <= 10:
        return math.log(abs(_gammafn_small(x)))
    if x > LGAMMA_XMAX:
        return np.inf
    if x > 1e17:
        return x * (math.log(x) - 1.0)
    if x > 4934720.0:
        return muladd(x - 0.5, math.log(x), M_LN_SQRT_2PI) - x
    return muladd(x - 0.5, math.log(x), M_LN_SQRT_2PI) - x + _lgammacor(x)


@njit(cache=True, error_model="numpy")
def _lgamma1p(a):
    """log(Gamma(a + 1)), accurate for small a: lgamma1p().

    Parameters: a (float). Returns: float.
    """
    if abs(a) >= 0.5:
        return _lgammafn(a + 1)
    lgam = LGAMMA1P_C * _logcf(-a / 2, 42.0, 1.0, TOL_LOGCF)
    for i in range(39, -1, -1):
        lgam = muladd(-a, lgam, _LGAMMA1P_COEFFS[i])
    return muladd(muladd(a, lgam, -EULERS_CONST), a, -_log1pmx(a))


@njit(cache=True, error_model="numpy")
def _stirlerr(n):
    """log(n!) minus Stirling's approximation of it: stirlerr().

    Parameters: n (float, > 0). Returns: float.
    """
    if n <= 23.5:
        nn = n + n
        if n <= 15.0 and nn == int(nn):
            return _SFERR_HALVES[int(nn)]
        if n <= 5.25:
            if n >= 1.0:
                l_n = math.log(n)
                return muladd(n, 1 - l_n, _libm_lgamma(n)) + (l_n - M_LN_2PI) * 0.5
            return muladd(-(n + 0.5), math.log(n), _lgamma1p(n)) + n - M_LN_SQRT_2PI
        if n > 12.8:
            k = 6
        elif n > 12.3:
            k = 7
        elif n > 8.9:
            k = 8
        elif n > 7.3:
            k = 10
        elif n > 6.6:
            k = 12
        elif n > 6.1:
            k = 14
        else:
            k = 16
    else:
        if n > 15.7e6:
            return _S[0] / n
        if n > 6180:
            k = 1
        elif n > 205:
            k = 2
        elif n > 86:
            k = 3
        elif n > 27:
            k = 4
        else:
            k = 5
    # (S0 - (S1 - ... (S[k-1] - S[k] / nn) / nn ...) / nn) / n
    nn = n * n
    v = _S[k] / nn
    for j in range(k - 1, 0, -1):
        v = (_S[j] - v) / nn
    return (_S[0] - v) / n


@njit(cache=True, error_model="numpy")
def _add1(d, yh, yl):
    """Add d to the pair (yh, yl), integer part to yh: the macro ADD1.

    Parameters: d, yh, yl (floats). Returns: (yh, yl).
    """
    d1 = np.floor(d + 0.5)
    d2 = d - d1
    return yh + d1, yl + d2


@njit(cache=True, error_model="numpy")
def _ebd0(x, M):
    """x log(x / M) + (M - x) in two parts: ebd0().

    Parameters: x, M (floats, >= 0). Returns: (yh, yl), the sum of which is
    the value.
    """
    yh = 0.0
    yl = 0.0
    if x == M:
        return yh, yl
    if x == 0:
        return M, yl
    if M == 0:
        return np.inf, yl
    if M / x == np.inf:
        return M, yl
    r, e = math.frexp(M / x)
    if M_LN2 * float(-e) > 1.0 + DBL_MAX / x:
        return np.inf, yl
    i = int(np.floor(muladd(r - 0.5, 256.0, 0.5)))
    f = np.floor(1024.0 / (0.5 + i / 256.0) + 0.5)
    fg = math.ldexp(f, -(e + 10))
    if fg == np.inf:
        return fg, yl
    yh, yl = _add1(-x * _log1pmx(muladd(M, fg, -x) / x), yh, yl)
    if fg == 1:
        return yh, yl
    for j in range(4):
        yh, yl = _add1(x * _BD0_SCALE[i, j], yh, yl)
        yh, yl = _add1(-x * _BD0_SCALE[0, j] * e, yh, yl)
        if not math.isfinite(yh):
            return np.inf, 0.0
    yh, yl = _add1(M, yh, yl)
    yh, yl = _add1(-M * fg, yh, yl)
    return yh, yl


@njit(cache=True, error_model="numpy")
def _dpois_raw(x, lam, give_log):
    """Poisson density at a real x >= 0: dpois_raw().

    Parameters: x (float), lam (mean, >= 0), give_log (bool).
    Returns: float.
    """
    zero = -np.inf if give_log else 0.0
    if lam == 0:
        if x == 0:
            return 0.0 if give_log else 1.0
        return zero
    if not math.isfinite(lam):
        return zero
    if x < 0:
        return zero
    if x <= lam * DBL_MIN:
        return -lam if give_log else math.exp(-lam)
    if lam < x * DBL_MIN:
        if not math.isfinite(x):
            return zero
        v = muladd(x, math.log(lam), -lam) - _lgammafn(x + 1)
        return v if give_log else math.exp(v)
    yh, yl = _ebd0(x, lam)
    yl += _stirlerr(x)
    large = x >= X_LRG
    r = M_SQRT_2PI * math.sqrt(x) if large else M_2PI * x
    if give_log:
        return -yl - yh - (math.log(r) if large else 0.5 * math.log(r))
    return math.exp(-yl) * math.exp(-yh) / (r if large else math.sqrt(r))


@njit(cache=True, error_model="numpy")
def _dpois_wrap(x_plus_1, lam, give_log):
    """dpois(x_plus_1 - 1, lam) for real arguments: dpois_wrap().

    Parameters: x_plus_1 (float), lam (float), give_log (bool).
    Returns: float.
    """
    if not math.isfinite(lam):
        return -np.inf if give_log else 0.0
    if x_plus_1 > 1:
        return _dpois_raw(x_plus_1 - 1, lam, give_log)
    if lam > abs(x_plus_1 - 1) * M_CUTOFF:
        v = -lam - _lgammafn(x_plus_1)
        return v if give_log else math.exp(v)
    d = _dpois_raw(x_plus_1, lam, give_log)
    if give_log:
        return d + math.log(x_plus_1 / lam)
    return d * (x_plus_1 / lam)


@njit(cache=True, error_model="numpy")
def _pgamma_smallx(x, alph, log_p):
    """Upper tail for x < 1 by a series: pgamma_smallx().

    Parameters: x (float in (0, 1)), alph (shape), log_p (bool).
    Returns: float.
    """
    total = 0.0
    c = alph
    n = 0.0
    while True:
        n += 1
        c *= -x / n
        term = c / (alph + n)
        total += term
        if not (abs(term) > DBL_EPSILON * abs(total)):
            break
    lf2 = muladd(alph, math.log(x), -_lgamma1p(alph))
    if log_p:
        return _log1_exp(math.log1p(total) + lf2)
    f2m1 = math.expm1(lf2)
    return -muladd(total, f2m1, total + f2m1)


@njit(cache=True, error_model="numpy")
def _pd_upper_series(x, y, log_p):
    """sum_n x^n / (y (y + 1) ... (y + n - 1)): pd_upper_series().

    Parameters: x, y (floats), log_p (bool). Returns: float.
    """
    term = x / y
    total = term
    while True:
        y += 1
        term *= x / y
        total += term
        if not (term > total * DBL_EPSILON):
            break
    return math.log(total) if log_p else total


@njit(cache=True, error_model="numpy")
def _pd_lower_cf(y, d):
    """Continued fraction for the scaled upper tail: pd_lower_cf().

    Parameters: y, d (floats). Returns: float.
    """
    if y == 0:
        return 0.0
    f0 = y / d
    if abs(y - 1) < abs(d) * DBL_EPSILON:
        return f0
    if f0 > 1.0:
        f0 = 1.0
    c2 = y
    c4 = d
    a1 = 0.0
    b1 = 1.0
    a2 = y
    b2 = d
    while b2 > SCALEFACTOR:
        a1 /= SCALEFACTOR
        b1 /= SCALEFACTOR
        a2 /= SCALEFACTOR
        b2 /= SCALEFACTOR
    i = 0.0
    of = -1.0
    f = 0.0
    while i < 200000:
        i += 1
        c2 -= 1
        c3 = i * c2
        c4 += 2
        a1 = muladd(c4, a2, c3 * a1)
        b1 = muladd(c4, b2, c3 * b1)

        i += 1
        c2 -= 1
        c3 = i * c2
        c4 += 2
        a2 = muladd(c4, a1, c3 * a2)
        b2 = muladd(c4, b1, c3 * b2)

        if b2 > SCALEFACTOR:
            a1 /= SCALEFACTOR
            b1 /= SCALEFACTOR
            a2 /= SCALEFACTOR
            b2 /= SCALEFACTOR

        if b2 != 0:
            f = a2 / b2
            if abs(f - of) <= DBL_EPSILON * _fmax2(f0, abs(f)):
                return f
            of = f
    return f  # not converged; R warns here


@njit(cache=True, error_model="numpy")
def _pd_lower_series(lam, y):
    """sum_n y (y - 1) ... (y - n) / lam^(n + 1): pd_lower_series().

    Parameters: lam, y (floats). Returns: float.
    """
    term = 1.0
    total = 0.0
    while y >= 1 and term > total * DBL_EPSILON:
        term *= y / lam
        total += term
        y -= 1
    if y != np.floor(y):
        f = _pd_lower_cf(y, lam + 1 - y)
        total = muladd(term, f, total)
    return total


@njit(cache=True, error_model="numpy")
def _dnorm(x):
    """Standard normal density: dnorm4(x, 0, 1, FALSE).

    Parameters: x (float). Returns: float.
    """
    if math.isnan(x):
        return x
    if not math.isfinite(x):
        return 0.0
    x = abs(x)
    if x >= DNORM_HUGE:
        return 0.0
    if x < 5:
        return M_1_SQRT_2PI * math.exp(-0.5 * x * x)
    if x > DNORM_ZERO:
        return 0.0
    x1 = np.rint(x * 65536.0) / 65536.0
    x2 = x - x1
    return M_1_SQRT_2PI * (
        math.exp(-0.5 * x1 * x1) * math.exp(muladd(-0.5, x2, -x1) * x2)
    )


@njit(cache=True, error_model="numpy")
def _pnorm_upper(x, log_p):
    """Upper tail of the standard normal: pnorm5(x, 0, 1, FALSE, log_p).

    Parameters: x (float), log_p (bool). Returns: float.
    """
    if math.isnan(x):
        return x
    if not math.isfinite(x):
        if x < 0:
            return 0.0 if log_p else 1.0
        return -np.inf if log_p else 0.0
    y = abs(x)
    if y <= 0.67448975:
        if y > DBL_EPSILON * 0.5:
            xsq = x * x
            xnum = _PNORM_A[4] * xsq
            xden = xsq
            for i in range(3):
                xnum = (xnum + _PNORM_A[i]) * xsq
                xden = (xden + _PNORM_B[i]) * xsq
        else:
            xnum = 0.0
            xden = 0.0
        temp = x * (xnum + _PNORM_A[3]) / (xden + _PNORM_B[3])
        ccum = 0.5 - temp
        return math.log(ccum) if log_p else ccum
    if y <= M_SQRT_32:
        xnum = _PNORM_C[8] * y
        xden = y
        for i in range(7):
            xnum = (xnum + _PNORM_C[i]) * y
            xden = (xden + _PNORM_D[i]) * y
        temp = (xnum + _PNORM_C[7]) / (xden + _PNORM_D[7])
        z = y
    elif (log_p and y < 1e170) or (-8.2924 < x and x < 38.4674):
        xsq = 1.0 / (x * x)
        xnum = _PNORM_P[5] * xsq
        xden = xsq
        for i in range(4):
            xnum = (xnum + _PNORM_P[i]) * xsq
            xden = (xden + _PNORM_Q[i]) * xsq
        temp = xsq * (xnum + _PNORM_P[4]) / (xden + _PNORM_Q[4])
        temp = (M_1_SQRT_2PI - temp) / y
        z = x
    else:
        if x > 0:
            return -np.inf if log_p else 0.0
        return 0.0 if log_p else 1.0
    # the macro do_del(z), then the tail wanted
    xsq = np.trunc(z * 16.0) / 16.0
    dl = (z - xsq) * (z + xsq)
    if log_p:
        if x > 0:
            return muladd(-xsq, xsq * 0.5, -(dl * 0.5)) + math.log(temp)
        return math.log1p(-math.exp(-xsq * (xsq * 0.5)) * math.exp(-(dl * 0.5)) * temp)
    cum = math.exp(-xsq * (xsq * 0.5)) * math.exp(-(dl * 0.5)) * temp
    if x > 0:
        return cum
    return 1.0 - cum


@njit(cache=True, error_model="numpy")
def _dpnorm(x, lower_tail, lp):
    """dnorm(x) / pnorm(x), given lp = log(pnorm(x)): dpnorm().

    Parameters: x (float), lower_tail (bool), lp (float). Returns: float.
    """
    if x < 0:
        x = -x
        lower_tail = not lower_tail
    if x > 10 and not lower_tail:
        term = 1 / x
        total = term
        x2 = x * x
        i = 1.0
        while True:
            term *= -i / x2
            total += term
            i += 2
            if not (abs(term) > DBL_EPSILON * total):
                break
        return 1 / total
    return _dnorm(x) / math.exp(lp)


@njit(cache=True, error_model="numpy")
def _ppois_asymp(x, lam, log_p):
    """P(Poisson(lam) <= x) by an asymptotic expansion: ppois_asymp() with
    lower_tail = TRUE.

    Parameters: x, lam (floats), log_p (bool). Returns: float.
    """
    dfm = lam - x
    pt_ = -_log1pmx(dfm / x)
    s2pt = math.sqrt(2 * x * pt_)
    if dfm < 0:
        s2pt = -s2pt

    res12 = 0.0
    res1_term = math.sqrt(x)
    res1_ig = res1_term
    res2_term = s2pt
    res2_ig = s2pt
    for i in range(1, 8):
        res12 = muladd(res1_ig, _PPOIS_A[i], res12)
        res12 = muladd(res2_ig, _PPOIS_B[i], res12)
        res1_term *= pt_ / i
        res2_term *= 2 * pt_ / (2 * i + 1)
        res1_ig = res1_ig / x + res1_term
        res2_ig = res2_ig / x + res2_term

    elfb = x
    elfb_term = 1.0
    for i in range(1, 8):
        elfb = muladd(elfb_term, _PPOIS_B[i], elfb)
        elfb_term /= x

    f = res12 / elfb
    np_ = _pnorm_upper(s2pt, log_p)
    if log_p:
        n_d_over_p = _dpnorm(s2pt, False, np_)
        return np_ + math.log1p(f * n_d_over_p)
    nd = _dnorm(s2pt)
    return muladd(f, nd, np_)


@njit(cache=True, error_model="numpy")
def _pgamma_raw_once(x, alph, log_p):
    """Upper tail of the gamma distribution with scale 1, one pass:
    pgamma_raw() without its final recomputation.

    Parameters: x (float), alph (shape, > 0), log_p (bool). Returns: float.
    """
    if x <= 0.0:
        return 0.0 if log_p else 1.0
    if x >= np.inf:
        return -np.inf if log_p else 0.0
    if x < 1:
        return _pgamma_smallx(x, alph, log_p)
    if x <= alph - 1 and x < 0.8 * (alph + 50):
        total = _pd_upper_series(x, alph, log_p)
        d = _dpois_wrap(alph, x, log_p)
        if log_p:
            return _log1_exp(d + total)
        return muladd(-d, total, 1.0)
    if alph - 1 < x and alph < 0.8 * (x + 50):
        d = _dpois_wrap(alph, x, log_p)
        if alph < 1:
            if x * DBL_EPSILON > 1 - alph:
                total = 0.0 if log_p else 1.0
            else:
                f = _pd_lower_cf(alph, x - (alph - 1)) * x / alph
                total = math.log(f) if log_p else f
        else:
            total = _pd_lower_series(x, alph - 1)
            total = math.log1p(total) if log_p else 1 + total
        return total + d if log_p else total * d
    return _ppois_asymp(alph - 1, x, log_p)


@njit(cache=True, error_model="numpy")
def _pgamma_upper(x, alph, rate):
    """pgamma(x, alph, rate = rate, lower.tail = FALSE) for scalars.

    Parameters: x (quantile), alph (shape), rate (floats). Returns: float.
    """
    scale = 1.0 / rate  # done by R's pgamma() before the C code
    if math.isnan(x) or math.isnan(alph) or math.isnan(scale):
        return x + alph + scale
    if alph < 0.0 or scale <= 0.0:
        return np.nan
    x /= scale
    if math.isnan(x):
        return x
    if alph == 0.0:
        return 1.0 if x <= 0 else 0.0
    res = _pgamma_raw_once(x, alph, False)
    if res < DBL_MIN / DBL_EPSILON:
        # too close to underflow: redone through the logarithm, as R does
        return math.exp(_pgamma_raw_once(x, alph, True))
    return res


@njit(cache=True, error_model="numpy")
def _branch(x, alph, rate):
    """Which branch of pgamma_raw() the scalar call takes.

    Parameters: x, alph, rate (floats). Returns: index into BRANCHES.
    """
    scale = 1.0 / rate
    if math.isnan(x) or math.isnan(alph) or math.isnan(scale):
        return 0
    if alph <= 0.0 or scale <= 0.0:
        return 0
    x /= scale
    if math.isnan(x) or x <= 0.0 or x >= np.inf:
        return 0
    if x < 1:
        return 1
    if x <= alph - 1 and x < 0.8 * (alph + 50):
        return 2
    if alph - 1 < x and alph < 0.8 * (x + 50):
        if alph < 1:
            return 6 if x * DBL_EPSILON > 1 - alph else 3
        return 4
    return 5


@njit(cache=True, error_model="numpy")
def _pgamma_upper_loop(x, alph, rate, out):
    """Fill out with _pgamma_upper() of each triple.

    Parameters: x, alph, rate, out (1-d float64 arrays of one length).
    Returns: None.
    """
    for i in range(x.size):
        out[i] = _pgamma_upper(x[i], alph[i], rate[i])


@njit(cache=True, error_model="numpy")
def _branch_loop(x, alph, rate, out):
    """Fill out with _branch() of each triple.

    Parameters: x, alph, rate (1-d float64 arrays), out (int64 array).
    Returns: None.
    """
    for i in range(x.size):
        out[i] = _branch(x[i], alph[i], rate[i])


def _broadcast(x, shape, rate):
    """Broadcast the three arguments to flat float64 arrays.

    Parameters: x, shape, rate (array-likes).
    Returns: (x, shape, rate, common shape).
    """
    x, shape, rate = np.broadcast_arrays(
        np.asarray(x, dtype=np.float64),
        np.asarray(shape, dtype=np.float64),
        np.asarray(rate, dtype=np.float64),
    )
    flat = [np.ascontiguousarray(a).ravel() for a in (x, shape, rate)]
    return flat[0], flat[1], flat[2], x.shape


def pgamma_upper(x, shape, rate):
    """Upper tail of the gamma distribution, the doubles R returns for
    pgamma(x, shape, rate = rate, lower.tail = FALSE).

    Parameters: x (quantiles), shape, rate (array-likes, broadcast against
    each other). Returns: float64 array of the broadcast shape; NaN where R
    gives NaN or NA (a missing argument, shape < 0, rate < 0 or infinite).
    """
    x, alph, rate, dims = _broadcast(x, shape, rate)
    out = np.empty(x.size, dtype=np.float64)
    _pgamma_upper_loop(x, alph, rate, out)
    return out.reshape(dims)


def pgamma_branch(x, shape, rate):
    """The branch of R's pgamma_raw() each call of pgamma_upper() takes.

    Parameters: x, shape, rate (array-likes, broadcast against each other).
    Returns: int64 array of indices into BRANCHES.
    """
    x, alph, rate, dims = _broadcast(x, shape, rate)
    out = np.empty(x.size, dtype=np.int64)
    _branch_loop(x, alph, rate, out)
    return out.reshape(dims)
