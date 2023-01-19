# Authors: Sean Bowe, Alessandro Chiesa, Matthew Green, Ian Miers, Pratyush Mishra, Howard Wu
#
# This script rigidly generates Montgomery or (twisted) Edwards curves over a given base prime field.
# Note that we ensure twist security
# Throughout, we write "Edwards" to mean "twisted Edwards".

# TODO: add function to check embedding degree, etc

# References:
# [BBJLP]: Bernstein, Birkner, Joye, Lange, Peters --- "Twisterd Edwards curves"
# [CLN]: Costello, Longa, Naehrig --- "A brief discussion on selecting new elliptic curves"
# [CS]: Costello, Smith --- "Montgomery curves and their arithmetic"

from sage.all import *
from multiprocessing import Pool
import time

def EnsureValidEdwards(F,a,d):
    # see [BBJLP, Thm 3.2]
    assert(F.characteristic() != 2)
    assert(a != F(0))
    assert(d != F(0))
    assert(a != d)
    assert(not d.is_square())

def EnsureOnEdwards(F,a,d,x,y):
    # see [BBJLP, Definition 2.1]
    EnsureValidEdwards(F,a,d)
    assert(a*x^2 + y^2 == 1 + d * x^2 * y^2)
    # TODO: no special cases right?

def EnsureOnWeierstrass(F,a,b,u,v):
    # see [BBJLP, Definition 2.1]
    assert(v^2 == u^3 + a * u + b)


def EnsureValidMontgomery(F,A,B):
    # see [BBJLP, Thm 3.2]
    assert(F.characteristic() != 2)
    assert(A != F(2))
    assert(A != F(-2))
    assert(B != F(0))

def EnsureOnMontgomery(F,A,B,u,v):
    # see [BBJLP, Definition 3.1]
    EnsureValidMontgomery(F,A,B)
    assert(B * v^2 == u^3 + A * u^2 + u)
    # TODO: how to handle special cases?
    #

# The caller is responsible for ensuring that ec is a curve that can be
# expressed in Montgomery form.
def Curve_WeierstrassToMontgomery(F, a, b, A_given=None, B_given=None):
    Falpha.<x> = F['x'];
    f = x^3 + a * x + b;
    for (alpha, _) in f.roots():
        s2 = F(3 * alpha^2 + a)^(-1)
        if s2.is_square():
            s = s2.sqrt()
            A = 3 * alpha * s
            B = s
            if A_given is None and B_given is None:
                return (A, B)
            elif (A_given == A or A_given == -A) and (B_given == B or B_given == -B):
                return (A_given, B_given)
    return None
    
# The caller is responsible for ensuring that ec is a curve that can be
# expressed in Montgomery form.
def Point_WeierstrassToMontgomery(F, a, b, u, v, A_given=None, B_given=None):
    # see [CS, Sec 2.4]
    EnsureOnWeierstrass(F, a, b, u, v)
    (A, B) = Curve_WeierstrassToMontgomery(F, a, b, A_given, B_given)
    x = u * B - A/F(3)
    y = v*B
    EnsureOnMontgomery(F, A, B, x, y)
    return (x, y)


def Curve_MontgomeryToWeierstrass(F, A, B):
    EnsureValidMontgomery(F,A,B)
    # see [CS, Sec 2.4]
    # From [CS], a Montgomery curve Bu^2 = v^3 + Av^2 + v is isomorphic to the
    # Weierstrass model curve y^2 = x^3 + (B^2 * (1 - A^2/3))x + (B^3 * A/3) * (2A^2/9 -1)
    # i.e., a short Weierstrass curve with 
    #    a = (B^2 * (1 - A^2/3)) and
    #    b = (B^3 * A/3) * (2A^2/9 -1)
    s = B
    alpha = A/(F(3) * B)
    a = s^(-2) - F(3) * alpha^2
    b = -alpha * (alpha^2 + a)
    assert(Curve_WeierstrassToMontgomery(F, a, b, A, B) == (A, B))
    return (a, b)

def Point_MontgomeryToWeierstrass(F, A, B, x, y):
    # see [CS, Sec 2.4]
    EnsureOnMontgomery(F, A, B, x, y)
    u = B * (x + A/F(3))
    v = B^2 * y 
    return (u, v)

def Curve_WeierstrassToEdwards(F, a, b, a_given=None, d_given=None):
    if a_given is not None and d_given is not None:
        return (a_given, d_given)
    else:
        (A, B) = Curve_WeierstrassToMontgomery(F, a, b)
        return Curve_MontgomeryToEdwards(F, A, B)

    
def Point_WeierstrassToEdwards(F, a, b, u, v, a_given = None, d_given = None):
    if a_given is not None and d_given is not None:
        (A, B) = Curve_EdwardsToMontgomery(F, a_given, d_given);
        (x, y) = Point_WeierstrassToMontgomery(F, a, b, u, v, A_given=A, B_given=B);
        return Point_MontgomeryToEdwards(F, A, B, x, y)
    else:
        (A, B) = Curve_WeierstrassToMontgomery(F, a, b)
        (x, y) = Point_WeierstrassToMontgomery(F, a, b, u, v, A, B)
        return Point_MontgomeryToEdwards(F, A, B, x, y)


def Curve_EdwardsToWeierstrass(F, a, d):
    (A, B) = Curve_EdwardsToMontgomery(F, a, d);
    return Curve_MontgomeryToWeierstrass(F, A, B)


def Point_EdwardsToWeierstrass(F, a, d, x, y):
    (x, y) = Point_EdwardsToMontgomery(F, a, d, x, y);
    (A, B) = Curve_EdwardsToMontgomery(F, a, d);
    (x, y) = Point_MontgomeryToWeierstrass(F, A, B, x, y)
    EnsureOnEdwards(F, a, d, x, y)
    return (x, y)

def Curve_EdwardsToMontgomery(F,a,d):
    # see [BBJLP, Thm 3.2]
    EnsureValidEdwards(F,a,d)
    A = F(2)*(a+d)/(a-d)
    B = F(4)/(a-d)
    EnsureValidMontgomery(F,A,B)
    return (A,B)

def Point_EdwardsToMontgomery(F,a,d,x,y):
    # see [BBJLP, Thm 3.2]
    EnsureOnEdwards(F,a,d,x,y)
    if (x == F(0) and y == F(-1)): return (F(0),F(0)) # see [BBJLP, p5]
    if (x == F(0) and y == F(1)): return (F(0),F(1)) # point at infinity, see [BBJLP, p5]
    u = (1+y)/(1-y)
    v = (1+y)/((1-y)*x)
    return (u,v) # TODO fix exceptional points

def Curve_MontgomeryToEdwards(F,A,B):
    # see [BBJLP, Thm 3.2]
    EnsureValidMontgomery(F,A,B)
    a = (A+F(2))/B
    d = (A-F(2))/B
    EnsureValidEdwards(F,a,d)
    return (a,d)

def Point_MontgomeryToEdwards(F,A,B,u,v):
    # see [BBJLP, Thm 3.2]
    EnsureOnMontgomery(F,A,B,u,v)
    if (u == F(0) and v == F(0)): return (F(0),F(-1)) # see [BBJLP, p5]
    if (u == F(0) and v == F(1)): return (F(0),F(1)) # point at infinity, see [BBJLP, p5]
    x = u/v
    y = (u-F(1))/(u+F(1))
    (a, d) = Curve_MontgomeryToEdwards(F, A, B)
    EnsureOnEdwards(F, a, d, x, y)
    return (x,y) # TODO fix exceptional points

def MontgomeryCurve(F,A,B):
    EnsureValidMontgomery(F,A,B)
    # Recall that:
    #
    #    EllipticCurve([a1,a2,a3,a4,a6])
    #     -> y^2 + a1 * x * y + a3 * y = x^3 + a2 * x^2 + a4 * x + a6
    (a, b) = Curve_MontgomeryToWeierstrass(F, A, B)
    ec = EllipticCurve(F,[0,0,0,a,b])
    return ec

def EdwardsCurve(F,a,d):
    EnsureValidEdwards(F,a,d)
    (A,B) = Curve_EdwardsToMontgomery(F,a,d)
    ec = MontgomeryCurve(F,A,B)
    return ec

def ExamineMontgomeryParameters(F,A,B):
    # avoids two exceptional points of order 2 (see [BBJLP, p5])
    if ((A+F(2))*(A-F(2))).is_square(): return False

    # avoids two exceptional points of order 4 (see [BBJLP, p5])
    if ((A-F(2))/B).is_square(): return False

    ec = MontgomeryCurve(F,A,B)
    o1 = ec.order() # order of one curve ...
    p = F.order() # order of the field
    t = p + 1 - o1 #    [and its trace of Frobenius]
    o2 = p + 1 + t # ... and its twist
    o1_8 = o1 % 8
    o2_8 = o2 % 8
    b1 = ((o1_8 == 4) and (o2_8 == 0))
    b2 = ((o1_8 == 0) and (o2_8 == 4))
    if not (b1 or b2):
        return False
    elif b1:
        r1 = o1 // 4
        r2 = o2 // 8
    elif b2:
        r1 = o1 // 8
        r2 = o2 // 4
    # print("r2: %s" % r2)
    if not (is_prime(r1) and is_prime(r2)): return False # only consider curves with co-factors (4,8) or (8,4)
    if is_prime(r1) and not is_prime(r2): print("r1 is prime but not r2")
    if is_prime(r2) and not is_prime(r1): print("r2 is prime but not r1")
    print("----")
    print("Found Montgomery curve with A=%s and B=%s" % (A, B))
    print("Additional information: =")
    print("Is -1 a square: %s" % (F(-1).is_square()))
    order = ec.order();
    if order % 8 == 0:
        r = order // 8
        cofactor = 8
        print("r = %s " % r)
        print("cofactor = %s" % cofactor)
    elif order % 4 == 0:
        r = order // 4
        cofactor = 4
        print("r = %s " % r)
        print("cofactor = %s" % cofactor)
    
    # Time to sample a generator for the prime order subgroup.
    gen = ec.gens()[0] # We are in an abelian group, so only one generator exists.
    print(gen)
    gen = cofactor * gen # Multiply by cofactor to get a point in the prime_group()
    print(gen)
    assert(gen.order() == r)
    # The generator's coeffcients are in Weierstrass coordinates; we need to
    # convert them to Montgomery coordinates before printing them.
    (sw_x, sw_y) = gen.xy();
    (a, b) = Curve_MontgomeryToWeierstrass(F, A, B);
    (mont_x, mont_y) = Point_WeierstrassToMontgomery(F, a, b, sw_x, sw_y, A, B);
    print("Generator: (%s, %s)" % (mont_x, mont_y));
    EnsureOnMontgomery(F, A, B, mont_x, mont_y)
    return True

def EdwardsSquareTest(F, a, d):
    # avoids two exceptional points of order 2 (see [BBJLP, p5])
    if (a*d).is_square(): return False

    # avoids two exceptional points of order 4 (see [BBJLP, p5])
    if d.is_square(): return False
    return True

def ExamineEdwardsParameters(F,a,d, check = True):
    if check:
        if not EdwardsSquareTest(F, a, d): return False
    

    ec = EdwardsCurve(F,a,d)
    o1 = ec.order() # order of one curve ...
    p = F.order() # order of the field
    t = p + 1 - o1 #    [and its trace of Frobenius]
    o2 = p + 1 + t # ... and its twist
    o1_8 = o1 % 8
    o2_8 = o2 % 8
    b1 = ((o1_8 == 4) and (o2_8 == 0))
    b2 = ((o1_8 == 0) and (o2_8 == 4))
    if (b1 and b2):
        return False
    elif not (b1 or b2):
        return False
    elif b1:
        r1 = o1 // 4
        r2 = o2 // 8
    elif b2:
        r1 = o1 // 8
        r2 = o2 // 4
    if not (is_prime(r1) and is_prime(r2)): return False # only consider curves with co-factors (4,8) or (8,4)
    print("----")
    print("Found Edwards curve with a=%s and d=%s" % (a, d))
    print("Additional information:")
    print("Is -1 a square: %s" % (F(-1).is_square()))
    print("q = %s mod 4" % (F.cardinality() % 4))
    order = ec.order();
    if order % 8 == 0:
        r = order // 8
        cofactor = 8
        print("r = %s " % r)
        print("cofactor = %s" % cofactor)
    elif order % 4 == 0:
        r = order // 4
        cofactor = 4
        print("r = %s " % r)
        print("cofactor = %s" % cofactor)
    # Time to sample a generator for the prime order subgroup.
    gen = ec.gens()[0] # We are in an abelian group, so only one generator exists.
    gen = cofactor * gen # Multiply by cofactor to get a point in the prime_group()
    # The generator's coeffcients are in Weierstrass coordinates; we need to
    # convert them to Edwards coordinates before printing them.
    (sw_x, sw_y) = gen.xy();
    (sw_a, sw_b) = Curve_EdwardsToWeierstrass(F, a, d);
    (te_x, te_y) = Point_WeierstrassToEdwards(F, sw_a, sw_b, sw_x, sw_y, a_given=a, d_given=d);
    EnsureOnEdwards(F, a, d, te_x, te_y)
    print("Generator: (%s, %s)" % (te_x, te_y));
    assert(gen.order() == r)

    return True

def ExamineMontgomeryParametersHelper(args):
    return ExamineMontgomeryParameters(*args)

def SearchMontgomeryParameters(q, num_threads):
    # only prime fields of odd characteristic
    assert (is_prime(q) and q != 2)
    Fq = GF(q)

    print("Desired base prime field is:")
    print("q = %s" % q)

    As = []
    for i in range (1, 10000):
        if i % 100 == 0:
            print(i)
        
        d = (Fq(1) + Fq(i)) * 4
        A = d + Fq(2)
        
        As.append((Fq, A, Fq(1)))
        As.append((Fq, -A, Fq(1)))

    p = Pool(num_threads)
    result = p.map_async(ExamineMontgomeryParametersHelper, As)
    while not result.ready():
        print("num left: {}".format(result._number_left * result._chunksize))
        time.sleep(100)
    real_result = result.get()
    p.close()
    p.join()
    
def ExamineEdwardsParametersHelper(args):
    return ExamineEdwardsParameters(*args)

def SearchEdwardsParameters(q, num_threads):
    # only prime fields of odd characteristic
    assert (is_prime(q) and q != 2)
    Fq = GF(q)

    q_4 = mod(q,4)

    print("Desired base prime field is:")
    print("q = %s" % q)
    print("(congruent to %s mod 4)" % q_4)

    if q_4 == 1: a = Fq(-1)
    if q_4 == 3: a = Fq(1)

    params = []
    for i in range (1, 50000):
        if EdwardsSquareTest(Fq, a, Fq(i)):
            params.append((Fq, a, Fq(i)))
        elif EdwardsSquareTest(Fq, a, -Fq(i)):
            params.append((Fq, a, -Fq(i)))
    print("Number of parameters being considered: %s" % len(params))

    p = Pool(num_threads)
    result = p.map_async(ExamineEdwardsParametersHelper, params)
    while not result.ready():
        print("num left: {}".format(result._number_left * result._chunksize))
        time.sleep(60)
    real_result = result.get()
    p.close()
    p.join()

def scaling(a, d, prime):
    Fp = GF(prime)
    if Fp(-a).is_square():
        f = sqrt ( Fp ( - a ))
        a_ = Fp ( a / ( f * f ))
        d_ = Fp ( d / ( - a ))
        if a_ == Fp (-1):
            a_ = -1
        else:
            a_ = a
            d_ = a
    return a_ , d_ , f

def StarkJubOldCheck():
    r = 3618502788666131213697322783095070105623107215331596699973092056135872020481
    Fr = GF(r)
    A = Fr(110430)
    B = Fr(1)

    ExamineMontgomeryParameters(Fr, A, B)
    (a, d) = Curve_MontgomeryToEdwards(Fr, A, B)
    print(a, d)

    ExamineEdwardsParameters(Fr, Fr(a), Fr(d))
    (A, B) = Curve_EdwardsToMontgomery(Fr, a, d)
    ExamineMontgomeryParameters(Fr, A, B)

def StarkJubCheckFoundNewEdwards():
    r = 3618502788666131213697322783095070105623107215331596699973092056135872020481
    Fr = GF(r)
    a = 3618502788666131213697322783095070105623107215331596699973092056135872020480
    d = 36659
    x0 = Fr(2192997259653830321980110858627110452394882798306613832002354879583275243830)
    y0 = Fr(1870298999876184600770095838537746999597803359289427334927271486896920468444)

    EnsureValidEdwards(Fr, a, d)
    EnsureOnEdwards(Fr, a, d, x0, y0)

    (a_, d_, f) = scaling(a, d, r)
    print(a_)
    print(d_)
    print(f)


def StarkJubNewCheck():
    r = 3618502788666131213697322783095070105623107215331596699973092056135872020481
    Fr = GF(r)
    A = Fr(146638)
    B = Fr(1)

    ExamineMontgomeryParameters(Fr, A, B)
    (a, d) = Curve_MontgomeryToEdwards(Fr, A, B)
    print(a, d)

    ExamineEdwardsParameters(Fr, Fr(a), Fr(d))
    (A, B) = Curve_EdwardsToMontgomery(Fr, a, d)
    ExamineMontgomeryParameters(Fr, A, B)

def StarkJubCheck():
    #r = 3618502788666131213697322783095070105623107215331596699973092056135872020481
    Fr = GF(r)
    A = Fr(110430)
    B = Fr(1)

    #ExamineMontgomeryParameters(Fr, A, B)
    #(a, d) = Curve_MontgomeryToEdwards(Fr, A, B)
    #print(a, d)
    #ExamineEdwardsParameters(Fr, Fr(a), Fr(d))
    a = 110432
    d = 110428
    ExamineEdwardsParameters(Fr, Fr(a), Fr(d))

    ExamineEdwardsParameters(Fr, Fr(a), Fr(d))
    (A, B) = Curve_EdwardsToMontgomery(Fr, a, d)
    ExamineMontgomeryParameters(Fr, A, B)


    print("---Examined Montgomery Parameters---")

def SearchTedwardsStarkParameters():
    q = 3618502788666131213697322783095070105623107215331596699973092056135872020481

    SearchEdwardsParameters(q, 3)

def EnsureValidEdwardsStarkJub():
    r = 3618502788666131213697322783095070105623107215331596699973092056135872020481
    Fr = GF(r)
    A = Fr(110430)
    B = Fr(1)

    EnsureValidEdwards(Fr, 110432, 110428)

def JubJubCheck_Less():
    r = 52435875175126190479447740508185965837690552500527637822603658699938581184513
    Fr = GF(r)
    A = Fr(40962)
    B = Fr(1)

    ExamineMontgomeryParameters(Fr, A, B)
    (a, d) = Curve_MontgomeryToEdwards(Fr, A, B)
    ExamineEdwardsParameters(Fr, Fr(a), Fr(d))

# JubJub sanity check
def JubJubCheck():
    #r is the order of the elliptic curve, but we suppose that it is the prime
    r = 52435875175126190479447740508185965837690552500527637822603658699938581184513
    Fr = GF(r)
    A = Fr(40962)
    B = Fr(1)

    ExamineMontgomeryParameters(Fr, A, B)
    (a, d) = Curve_MontgomeryToEdwards(Fr, A, B)
    ExamineEdwardsParameters(Fr, Fr(a), Fr(d))

    # Zcash JubJub test parameters
    # https://github.com/daira/jubjub
    a = Fr(-1)
    d = -Fr(10240)/Fr(10241)
    x0 = Fr(11076627216317271660298050606127911965867021807910416450833192264015104452986)
    y0 = Fr(44412834903739585386157632289020980010620626017712148233229312325549216099227)
    x1 = Fr(8076246640662884909881801758704306714034609987455869804520522091855516602923)
    y1 = Fr(13262374693698910701929044844600465831413122818447359594527400194675274060458)

    #ExamineEdwardsParameters(Fr, a, d)
    (A, B) = Curve_EdwardsToMontgomery(Fr, a, d)
    print(A)
    print(B)
    #ExamineMontgomeryParameters(Fr, A, B)

    (u0, v0) = Point_EdwardsToMontgomery(Fr, a, d, x0, y0)
    (nx0, ny0) = Point_MontgomeryToEdwards(Fr, A, B, u0, v0)
    assert(x0 == nx0)
    assert(y0 == ny0)

    (u1, v1) = Point_EdwardsToMontgomery(Fr, a, d, x1, y1)
    (nx1, ny1) = Point_MontgomeryToEdwards(Fr, A, B, u1, v1)
    assert(x1 == nx1)
    assert(y1 == ny1)

#StarkJubOldCheck()
StarkJubNewCheck()
#StarkJubCheckFoundNewEdwards() 
#EnsureValidEdwardsStarkJub()
#JubJubCheck_Less()
#StarkJubCheck()
#SearchTedwardsStarkParameters()

#JubJubCheck()
# q = 52435875175126190479447740508185965837690552500527637822603658699938581184513
# SearchEdwardsParameters(q, 1)

# q = 8444461749428370424248824938781546531375899335154063827935233455917409239041
# SearchEdwardsParameters(q, 3)

#q = 258664426012969094010652733694893533536393512754914660539884262666720468348340822774968888139573360124440321458177
#Fq = GF(q)
#ExamineEdwardsParameters(Fq, Fq(-1), Fq(79743))

q = 22369874298875696930346742206501054934775599465297184582183496627646774052458024540232479018147881220178054575403841904557897715222633333372134756426301062487682326574958588001132586331462553235407484089304633076250782629492557320825577
Fq = GF(q)
sw_a = Fq(5)
sw_b = Fq(17764315118651679038286329069295091506801468118146712649886336045535808055361274148466772191243305528312843236347777260247138934336850548243151534538734724191505953341403463040067571652261229308333392040104884438208594329793895206056414)

(a, d) = Curve_WeierstrassToEdwards(Fq, sw_a, sw_b)
print("Found curve with a, d:")
print((a, d))
print("-a is a square? %s" % (-a).is_square())
print("d is a square? %s" % d.is_square())
print("a*d is a square? %s" % (a*d).is_square())
print("q mod 4: %d" % mod(q, 4))

result = ExamineEdwardsParameters(Fq, a, d, check=False)
print(result)
