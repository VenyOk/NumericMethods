from math import sqrt

eps = 1e-3
a, b = 0, 3
def f(x):
    return x * sqrt(x + 1)

def trapezoid(f, n):
    h = (b - a) / n
    return h * ((f(a) + f(b)) / 2 + sum(f(a + i * h) for i in range(1, n)))

def simpson(f, n):
    h = (b - a) / n
    return h / 3 * (f(a) + f(b) + 4 * sum(f(a + h * i) for i in range(1, n, 2)) + 2 * sum(f(a + h * i) for i in range(2, n, 2)))
def richardson_error(I1, I2, k):
    return abs(I1 - I2) / (2 ** k - 1)

def found_n_for_method(func, k):
    n, R = 1, 1
    Ih = 0
    while eps <= abs(R):
        n += 2
        Ih_prev = Ih
        Ih = func(f, n)
        R = richardson_error(Ih, Ih_prev, k)
    return n - 1
# f(x) = x * sqrt(x + 1)
def main():
    I = 116 / 15
    n1 = found_n_for_method(simpson, 4)
    n2 = found_n_for_method(trapezoid, 2)

    S1 = simpson(f, n1)
    S2 = simpson(f, 2 * n1)
    R = richardson_error(S1, S2, 4)
    print("Метод Симпсона")
    print(f"I = {I}")
    print(f"n = {2 * n1}")
    print(f"I* = {S2}")
    print(f"R = {R}")
    print(f"I* + R = {S2 + R}")

    print()
    T1 = trapezoid(f, n2)
    T2 = trapezoid(f, 2 * n2)
    R2 = richardson_error(T1, T2, 2)
    print("Метод трапеций")
    print(f"I = {I}")
    print(f"n = {n2}")
    print(f"I* = {T2}")
    print(f"R = {R2}")
    print(f"I* + R = {T2 + R2}")
if __name__ == "__main__":
    main()