from math import sqrt, e

eps = 1e-3
a, b = 0, 3
I_exact = 116 / 15
def f(x):
    return x * sqrt(x + 1)
def rectangle(f, n):
    h = (b - a) / n
    return h * sum(f(a + (i - 0.5) * h) for i in range(1, n + 1))

def trapezoid(f, n):
    h = (b - a) / n
    return h / 2 * (f(a) + f(b) + 2 * sum(f(a + i * h) for i in range(1, n)))

def simpson(f, n):
    h = (b - a) / n
    return h / 3 * (f(a) + f(b) +
                    4 * sum(f(a + h * i) for i in range(1, n, 2)) +
                    2 * sum(f(a + h * i) for i in range(2, n, 2)))

def richardson_error(I1, I2, k):
    return (I2 - I1) / (2 ** k - 1)

def found_n_for_method(method, k):
    n = 2
    Ih = 0
    R = 1
    while abs(R) > eps:
        n *= 2
        Ih_2 = Ih
        Ih = method(f, n)
        R = richardson_error(Ih, Ih_2, k)
    return n

def main():

    n_simpson   = found_n_for_method(simpson, 4)
    n_trapezoid = found_n_for_method(trapezoid, 2)
    n_rectangle = found_n_for_method(rectangle, 2)
    
    S1 = simpson(f, n_simpson)
    S2 = simpson(f, 2 * n_simpson)
    R_simpson = richardson_error(S1, S2, 4)
    
    T1 = trapezoid(f, n_trapezoid)
    T2 = trapezoid(f, 2 * n_trapezoid)
    R_trapezoid = richardson_error(T1, T2, 2)
    
    R1 = rectangle(f, n_rectangle)
    R2 = rectangle(f, 2 * n_rectangle)
    R_rectangle = richardson_error(R1, R2, 2)
    
    print(f"eps = {eps}, I = {I_exact}")

    headers = ["", "Симпсон\t", "Трапеции\t", "Прямоугольники\t"]
    
    col_width = 20
    header_row = "".join(f"{h:<{col_width}}" for h in headers)
    print(header_row)
    print("-" * col_width * len(headers))
    
    rows = []

    rows.append((
        "n",
        f"{n_simpson:<{col_width}}\t",
        f"{n_trapezoid:<{col_width}}\t",
        f"{n_rectangle:<{col_width}}\t"
    ))

    rows.append((
        "I*",
        f"{S2:<{col_width}}\t",
        f"{T2:<{col_width}}\t",
        f"{R2:<{col_width}}\t"
    ))

    rows.append((
        "R",
        f"{R_simpson:<{col_width}}\t",
        f"{R_trapezoid:<{col_width}}\t",
        f"{R_rectangle:<{col_width}}\t"
    ))

    rows.append((
        "I* + R",
        f"{S2 + R_simpson:<{col_width}}\t",
        f"{T2 + R_trapezoid:<{col_width}}\t",
        f"{R2 + R_rectangle:<{col_width}}\t"
    ))
    
    for row in rows:
        print(f"{row[0]:<{col_width}}" + "".join(f"{cell:<{col_width}}" for cell in row[1:]))

if __name__ == "__main__":
    main()
