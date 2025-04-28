from math import sqrt
import matplotlib.pyplot as plt
import numpy as np
from time import time

eps = 1e-3
def f(x):
    return x ** 3 + 3 * x ** 2 - 2

def df(x):
    return 3 * x ** 2 + 6 * x

def ddf(x):
    return 6 * x + 6

def bisect_main(a, b, n):
    ans = dict()
    p = (b - a) / n
    a2 = a
    b2 = a + p
    for _ in range(n):
        ans1, ans2 = bisect(a2, b2)
        if ans1 is not None and ans2 is not None:
            ans[ans1] = ans2
        a2 += p
        b2 += p
    return ans

def bisect(a, b):
    cnt = 0
    for _ in range(1000):
        cnt += 1
        c = (a + b) / 2
        fc = f(c)
        if abs(fc) < eps and abs(df(c)) > 1e-2 and abs(ddf(c)) > 1e-2:
            return c, cnt
        if f(a) * f(c) < 0:
            b = c
        else:
            a = c
    return None, None

def newton(x0):
    x = x0
    cnt = 0
    for _ in range(1000):
        fx = f(x)
        dfx = df(x)
        cnt += 1
        if abs(dfx) < 1e-12:
            return None, cnt
        x_new = x - fx / dfx
        if abs(x_new - x) < 1e-6 and df(x_new) != 0 and ddf(x_new) != 0:
            return x_new, cnt
        x = x_new
    return None, cnt

def newton_main(a, b, n):
    ans = dict()
    p = (b - a) / n
    x = a
    found_roots = []

    for _ in range(n):
        root, cnt = newton(x)
        if root is not None:
            if not any(abs(root - r) < eps for r in found_roots):
                found_roots.append(root)
                ans[root] = cnt
        x += p
    return ans

def plot_function_with_roots(roots_bisect, roots_newton):
    x = np.linspace(-5, 5, 1000)
    y = f(x)

    plt.figure(figsize=(10, 6))
    plt.plot(x, y, label='f(x)', color='blue')
    plt.axhline(0, color='gray', linewidth=0.8, linestyle='--')

    for root in roots_bisect:
        plt.plot(root, f(root), 'ro', label='Bisection root' if 'Bisection root' not in plt.gca().get_legend_handles_labels()[1] else "")

    for root in roots_newton:
        plt.plot(root, f(root), 'go', label='Newton root' if 'Newton root' not in plt.gca().get_legend_handles_labels()[1] else "")

    plt.title('График функции f(x) с корнями')
    plt.xlabel('x')
    plt.ylabel('f(x)')
    plt.grid(True)
    plt.legend()
    plt.show()

def main():
    t1 = time()
    result_bisect = bisect_main(-5, 5, 10)
    t2 = time()
    print("Метод деления отрезка пополам")
    for k, v in result_bisect.items():
        error = abs(f(k))
        print(f"корень = {k:.6f}, итерации = {v}, погрешность = {error:.6e}, время = {(t2 - t1):.6f} сек")

    print()
    t1 = time()
    result_newton = newton_main(-5, 5, 20)
    t2 = time()
    print("Метод Ньютона")
    for k, v in result_newton.items():
        error = abs(f(k))
        print(f"корень = {k:.6f}, итерации = {v}, погрешность = {error:.6e}, время = {(t2 - t1):.6f} сек")

    plot_function_with_roots(result_bisect.keys(), result_newton.keys())


if __name__ == "__main__":
    main()
