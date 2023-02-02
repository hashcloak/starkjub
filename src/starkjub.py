from math import log2
import secrets

def inv(x, p):
    return pow(x, -1, p)

class StarkJub:
    def __init__(self):
        self.p = pow(2, 251) + 17 * pow(2, 192) + 1
        self.a = 146640
        self.d = 146636
        self.identity = Point(0, 1, self)
        self.O = 452312848583266401712165347886883763197416885958242462530951491185349408851
        

    def point(self, x, y):
        return Point(x, y, self)

class Point:
    def __init__(self, x, y, EC):
        self.x = x
        self.y = y
        self.EC = EC

        if not Point.validate(self):
            raise Exception(f"Provided coordinates {self} don't form a point on that curve")

    def __neg__(self):
        x, y, p = self.x, self.y, self.p
        if x == 0 and y == 1:
            return self
        return Point((-x) % p, y, self.EC)


    def __add__(self, other):
        x1, y1, a, d, p = self.x, self.y, self.EC.a, self.EC.d, self.EC.p
        x2, y2 = other.x, other.y
        if self.EC != other.EC: raise Exception('You cannot add points on different curves')
        s = (d * x1 * x2 * y1 * y2) % p
        x = ((x1 * y2 + y1 * x2) * inv(1 + s, p)) % p
        y = ((y1 * y2 - a * x1 * x2) * inv(1 - s, p)) % p
        return Point(x, y, self.EC)

    def __rmul__(self, other):
        if not isinstance(other, int) and not isinstance(self, Point):
            raise Exception('You can multiply only point by integer')
        else:
            res = self.EC.identity
            addend = self
            
            for i in bin(other)[2:][::-1]:
                if i == "1":
                    res = res + addend
                addend = addend + addend
            return res
                

    def __mul__(self, other):
        return self.__rmul__(other)

    def __sub__(self, other):
        return self + (-other)

    def __eq__(self, other):
        if (self.x == other.x) and (self.y == other.y):
            return True
        return False

    def __repr__(self):
        return f'({self.x}, {self.y})'

    def double(self):
        return self + self

    def validate(self):
        x, y, a, d, p = self.x, self.y, self.EC.a, self.EC.d, self.EC.p
        return ((a*x**2+y**2-1-d*x**2*y**2) % p == 0 and
                (0 <= x < p) and (0 <= y < p))

    def print_validation(self):
        if self.validate():
            print("point")
        else:
            print("not a point")

    def compress(self):
        y_bin = bin(self.y)[2:].zfill(self.EC.b)
        if self.x & 1 == 1:
            new_y = '1' + y_bin[1:]
        else:
            new_y = '0' + y_bin[1:]
        return int(new_y, 2)  # .to_bytes(int(log2(b)), byteorder='little')

    @staticmethod
    def decompress(compressed_point, EC):
        # y = int.from_bytes(compressed_point, byteorder='little')
        p, d, b = EC.p, EC.d, EC.b
        y_msb = int(bin(compressed_point)[2:].zfill(b)[0])
        y = int('0' + bin(compressed_point)[2:].zfill(b)[1:], 2)
        u = (y**2-1) % p
        v = ((d*y**2)+1) % p
        z = ((u*v**3)*pow((u*v**7), ((p-5)*inv(8,p))%p, p)) % p
        vz2 = (v*z**2) % p
        if vz2 == u:
            x = z
        elif vz2 == -u % p:
            x = (z*pow(2, ((p-1)*inv(4,p))%p, p)) % p
        else:
            raise Exception('Error in decompression formulas')
        if y_msb != x & 1:
            x = -x % p
        return EC.point(x, y)

def main():
    ECC = StarkJub()
    P = Point(2065699795511519733436237708177164622668357918131020778486714673024550645584, 1666035597895264107928948444893966434436309134596180408598119672656400359905, ECC)    
    print(P)
    P.print_validation()

    G = P + P
    print(G)
    G.print_validation()

    Double = P.double()
    print(Double)
    Double.print_validation()

    assert(G == Double)

    C = P * 11
    print(C)
    C.print_validation()

    D11 = P * 4 + 7 * P
    print(D11)
    D11.print_validation()

    assert(C == D11)

    MAX = P * (P.EC.O + 1)
    print(MAX)
    MAX.validate()
    MAX.print_validation()

if __name__ == "main":
    main()
