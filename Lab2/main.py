import copy
import matplotlib.pyplot as plt
import numpy as np

def count_f(x, func):
    return eval(func)

def count_spline(x, A, B, C, D, X, n):
    for i in range(n - 1):
        if X[i] <= x <= X[i + 1]:
            return A[i] + B[i] * (x - X[i]) + C[i] * (x - X[i]) ** 2 + D[i] * (x - X[i]) ** 3
    return None

def count_c(y, n, h):
    A = [[0] * n for _ in range(n)]
    B = [0] * n

    for i in range(1, n - 1):
        A[i][i - 1] = 1
        A[i][i] = 4
        A[i][i + 1] = 1
        B[i] = (y[i + 1] - 2 * y[i] + y[i - 1]) / (h ** 2)

    A[0][0] = 1
    B[0] = 0
    A[n - 1][n - 1] = 1
    B[n - 1] = 0

    return A, B

def solve(A, B):
    n = len(A)

    for i in range(n):
        max_row = i
        for k in range(i + 1, n):
            if abs(A[k][i]) > abs(A[max_row][i]):
                max_row = k

        A[i], A[max_row] = A[max_row], A[i]
        B[i], B[max_row] = B[max_row], B[i]

        for k in range(i + 1, n):
            c = -A[k][i] / A[i][i]
            for j in range(i, n):
                if i == j:
                    A[k][j] = 0
                else:
                    A[k][j] += c * A[i][j]
            B[k] += c * B[i]

    X = [0 for _ in range(n)]
    for i in range(n - 1, -1, -1):
        X[i] = B[i] / A[i][i]
        for k in range(i - 1, -1, -1):
            B[k] -= A[k][i] * X[i]

    return X

def main():
    a, b = map(int, input("Введите границы интервала: ").split())
    func = input("Введите функцию f(x): ")
    n = 32
    h = (b - a) / n
    X = [a + i * h for i in range(n)]
    print("X =", *X)
    Y = [count_f(X[i], func) for i in range(n)]
    M1, M2 = count_c(Y, n, h)
    C = solve(M1, M2)
    print("C =", *C)
    A = [Y[i] for i in range(n)]
    B = [(Y[i + 1] - Y[i] - h * (C[i + 1] + 2 * C[i])) / (3 * h) for i in range(n - 1)]
    D = [(C[i + 1] - C[i]) / (2 * h) for i in range(n - 1)]
    print("A =", *A)
    print("=" * 20)

    print("B =", *B)
    print("=" * 20)

    print("C =", *C)
    print("=" * 20)

    print("D =", *D)
    print("=" * 20)

    x = float(input("Значение точки: "))
    print("Значение функции:", count_f(x, func))
    print("Интерполированное значение:", count_spline(x, A, B, C, D, X, n))

    x_vals = X
    y_vals = [count_f(x, func) for x in x_vals]
    spline_vals = [count_spline(x, A, B, C, D, X, n) for x in x_vals]

    plt.figure(figsize=(10, 6)) 
    plt.plot(x_vals, y_vals, label='Исходная функция')
    plt.plot(x_vals, spline_vals, label='Сплайн', linestyle='--')
    plt.scatter(X, Y, color='red', label='Узлы')
    plt.title('График')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend()
    plt.grid(True)
    plt.show()

if __name__ == "__main__":
    main()
