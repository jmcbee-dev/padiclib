import math
from fractions import Fraction
from mathutil import extended_gcd


# Returns the decimal value of a given periodic string
# By default strings are read right-left
def find_padic_number(periodic_string, base, read_right_to_left=True):
    l = len(periodic_string)
    p_str = str(periodic_string)

    sum = 0
    if read_right_to_left:
        for idx, c, in enumerate(periodic_string[::-1]):
            sum += int(c) * (base ** idx)
    else:
        for idx, c, in enumerate(periodic_string):
            sum += int(c) * (base ** idx)

    return - sum, (base ** l - 1)


class Padic:
    def __init__(self, frac_numerator=None, frac_denominator=None, dec_number=None, prime_base=None, periodic_string=None, read_left_right=False):
        if isinstance(dec_number, float):
            self.dec_number = dec_number.as_integer_ratio()
        else:
            self.dec_number = None
        self.base = prime_base
        self.periodic_string = periodic_string
        self.frac_numerator = frac_numerator
        self.frac_denominator = frac_denominator
        self.periodic_start = None
        self.read_left_right = read_left_right
        self.finite_string = None
        self.periodic_n = None
        self.periodic_d = None
        self.point_string = None
        self.reduce_fraction()

    def is_padic_integer(self):
        return not (self.frac_denominator % self.base == 0)
    '''
    Finding a digit d for a p-adic expansion to a rational number m/n is equivalent to the following
    m / n === d (mod p)
    m === n * d (mod p)
    
    Let m = m * 1,
    1 = n * d + p * y
    Solve the diophantine equation for d, y using extended euclidian
    
    Finally
    m === m * n * d + m * p * y (mod p)
    '''
    def find_digit(self, numerator, denominator):
        _, x, _ = extended_gcd(denominator % self.base, self.base)
        return int(x * numerator % self.base)

    '''
    After finding a digit d, we need to repeat the process by first
    finding the next divisor q. That is
    q_n = d + p ( q_{n+1} ) 
    q_n = m / n
    ->
    ( m - n * d ) / p = z
    z = new numerator for next divisor
    '''
    def find_next_dividend(self, q_0_num, q_0_den, d):
        # print(q_0_num, q_0_den, d)
        return (q_0_num - q_0_den * d) / self.base

    def find_string(self):
        self.__find_finite_string()
        self.__find_periodic_string()
        return

    def reduce_fraction(self):
        n, d = self.frac_numerator, self.frac_denominator
        while n % d == 0 and d != 1:
            n = n // d
            d = d // d
        self.frac_numerator = n
        self.frac_denominator = d
        return

    '''
    For a given number, if its fractional representation n/d > 0 or n/d < -1, 
    then the p-adic representation will have a finite non-repeating string 
    '''
    def __find_finite_string(self):
        n, d = self.frac_numerator, self.frac_denominator
        self.finite_string = ""

        if -1 <= n / d < 0:
            return

        while n / d > 0 or n / d < -1:
            digit = self.find_digit(n, d)
            self.finite_string = str(digit) + self.finite_string
            n = self.find_next_dividend(n, d, digit)
        self.periodic_n = n
        return

    '''
    To find a periodic string for a given number, we need to find the digits
    of the string. After finding the digits using a divisor q_n,
    we update to q_{n+1}. If we have already seen q_{n+1}, then we are
    guaranteed to be repeating the string and have found the periodic string.
    
    The first divisor q_0 is simply the input number to the class constructor. 
    '''
    def __find_periodic_string(self):
        if (self.dec_number is None and (self.frac_numerator is None or self.frac_denominator is None)) or self.base is None:
            print("Set the base and decimal number first before finding the periodic string")
            return

        if not self.is_padic_integer():
            print("Fraction denominator is divisible by p, therefore there is no periodic string")
            return

        if self.periodic_string is None:
            self.periodic_string = ""
        else:
            return self.periodic_string

        if self.periodic_n is None:
            self.periodic_n = self.frac_numerator

        # Mark the beginning of the periodic string
        n, d = self.periodic_n, self.frac_denominator
        q_mark = n
        digit = self.find_digit(n, d)
        self.periodic_string = str(digit) + self.periodic_string
        n = self.find_next_dividend(n, d, digit)

        while n != q_mark:
            digit = self.find_digit(n, d)
            self.periodic_string = str(digit) + self.periodic_string
            n = self.find_next_dividend(n, d, digit)
        return

    def valuation(self, integer):
        if integer == 0:
            return math.inf
        else:
            val = 0
            while integer % self.base == 0:
                val += 1
                integer = integer // self.base
            return val

    def find_valuation(self):
        numerator_val = self.valuation(self.frac_numerator)
        denominator_val = self.valuation(self.frac_denominator)
        return numerator_val - denominator_val

    def abs_value(self):
        val = self.find_valuation()
        return f'{self.base} ^ ({-val})'

    def distance(self, other: "Padic"):
        tmp_adic = self - other
        return tmp_adic.abs_value()

    def print_string(self):
        if not self.read_left_right:
            if self.periodic_string is not None:
                print(f'...({self.periodic_string}){self.finite_string}')
            else:
                print(f'{self.finite_string}')
        else:
            if self.periodic_string is not None:
                print(f'{self.finite_string[::-1]}({self.periodic_string[::-1]})...')
            else:
                print(f'{self.finite_string[::-1]}')

        # if not self.read_left_right:
        #     if self.periodic_start != 0:
        #         print(f". . .({self.periodic_string[:-self.periodic_start]}){self.periodic_string[-self.periodic_start:]}")
        #     else:
        #         print(f". . .({self.periodic_string})")
        # else:
        #     if self.periodic_start != 0:
        #         print(f"{self.periodic_string[-1:-(1+self.periodic_start):-1]}({self.periodic_string[-(1+self.periodic_start)::-1]}). . .")
        #     else:
        #         print(f"({self.periodic_string[::-1]}). . .")

    def print_frac(self):
        print(f'{self.frac_numerator} / {self.frac_denominator}')

    def __sub__(self, other: "Padic"):
        if self.base != other.base:
            print("Different p base between numbers")
            return
        numerator = (self.frac_numerator * other.frac_denominator) - (self.frac_denominator * other.frac_numerator)
        denominator = self.frac_denominator * other.frac_denominator
        return Padic(frac_numerator=numerator, frac_denominator=denominator, prime_base=self.base)

    def __add__(self, other: "Padic"):
        if self.base != other.base:
            print("Different p base between numbers")
            return
        numerator = (self.frac_numerator * other.frac_denominator) + (self.frac_denominator * other.frac_numerator)
        denominator = self.frac_denominator * other.frac_denominator
        return Padic(frac_numerator=numerator, frac_denominator=denominator, prime_base=self.base)

