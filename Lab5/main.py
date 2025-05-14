import math

eps = 1e-3
max_iter = 10000

def f(x, y):
    return (x - 1) ** 2 + (y - 2) ** 4 + math.exp(1 / (x * x + y * y))

def f_x(x, y):
    r2 = x * x + y * y
    return 2 * (x - 1) - (2 * x / r2 ** 2) * math.exp(1 / r2)

def f_y(x, y):
    r2 = x * x + y * y
    return 4 * (y - 2) ** 3 - (2 * y / r2 ** 2) * math.exp(1 / r2)

def f_x_x(x, y):
    r2 = x * x + y * y
    exp1 = math.exp(1 / r2)
    q = 2 * x / r2 ** 2
    dq_dx = 2 / r2 ** 2 - 8 * x * x / r2 ** 3
    return 2 - exp1 * (dq_dx - q * q)

def f_y_y(x, y):
    r2 = x * x + y * y
    exp1 = math.exp(1 / r2)
    s = 2 * y / r2 ** 2
    ds_dy = 2 / r2 ** 2 - 8 * y * y / r2 ** 3
    return 12 * (y - 2) ** 2 - exp1 * (ds_dy - s * s)

def f_x_y(x, y):
    r2 = x * x + y * y
    exp1 = math.exp(1 / r2)
    return exp1 * (8 * x * y / r2 ** 3 + 4 * x * y / r2 ** 4)

def analytical_min():
    return 1.029180, 2.319519

def solve(xk, yk):
    k = 0
    while max(abs(f_x(xk, yk)), abs(f_y(xk, yk))) >= eps and k < max_iter:
        gx, gy = f_x(xk, yk), f_y(xk, yk)
        phi1 = - (gx * gx + gy * gy)
        gHg  = (f_x_x(xk, yk) * gx * gx +
                2 * f_x_y(xk, yk) * gx * gy +
                f_y_y(xk, yk) * gy * gy)
        t_star = - phi1 / gHg
        xk -= t_star * gx
        yk -= t_star * gy
        k += 1
    return k, xk, yk

def main():
    xk, yk = 1.0, 2.0
    print(f"Стартовая точка: ({xk}, {yk})")
    it, xk, yk = solve(xk, yk)
    f1 = f(xk, yk)
    f2 = f(*analytical_min())
    print(f'Количество итераций: {it}')
    print(f"Метод: ({xk:.6f}, {yk:.6f}), f = {f1:.6f}")
    print(f"Аналитически: ({analytical_min()[0]:.6f}, {analytical_min()[1]:.6f}), f = {f2:.6f}")
    print(f'Разность координат: ({abs(xk - analytical_min()[0]):.6e}, {abs(yk - analytical_min()[1]):.6e})')
    print(f'Разность значений функции: {abs(f1 - f2):.6e}')

if __name__ == "__main__":
    main()