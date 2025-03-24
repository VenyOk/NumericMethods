import csv
from math import sqrt

def count_f(x):
    return x * sqrt(x + 1)

def count_spline(x, A, B, C, D, X, i):
    return A[i] + B[i] * (x - X[i]) + C[i] * (x - X[i]) ** 2 + D[i] * (x - X[i]) ** 3

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

def count_c(A, B, C, D):
    n = len(D)
    matrix = [[0.0] * n for _ in range(n)]
    for i in range(n):
        if i > 0:
            matrix[i][i - 1] = A[i - 1]
        matrix[i][i] = B[i]
        if i < n - 1:
            matrix[i][i + 1] = C[i]
    return solve(matrix, D)

def write_to_csv(filename, headers, data):
    with open(filename, 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(headers)
        writer.writerows(data)

def main():
    a, b = 0, 3
    n = 32
    h = (b - a) / n
    X = [a + i * h for i in range(n)]
    Y = list(map(count_f, X))
    coef_a = [0] * (n-1)
    coef_b = [0] * n
    coef_c = [0] * (n-1)
    coef_d = [0] * n
    coef_b[0] = 1.0
    coef_b[n-1] = 1.0
    for i in range(1, n-1):
        coef_a[i - 1] = h
        coef_b[i] = 4 * h
        coef_c[i] = h
        coef_d[i] = 3 * ((Y[i+1] - Y[i]) / h - (Y[i] - Y[i-1]) / h)
    C = count_c(coef_a, coef_b, coef_c, coef_d)

    A = [Y[i] for i in range(n)]
    B = [((Y[i + 1] - Y[i]) / h) - h / 3 * (C[i + 1] + 2 * C[i])
         for i in range(n - 1)]
    B.append(((Y[n-1] - Y[n-2]) / h) - (2 * h * C[n-2] / 3))
    D = [(C[i + 1] - C[i]) / (3 * h) for i in range(n - 1)]
    D.append(C[n-1] / (3 * h))
    headers = ["i", "X", "A", "B", "C", "D"]
    data = []
    for i in range(n):
        row = [i, X[i], A[i], B[i], C[i], D[i]]
        data.append(row)

    write_to_csv('output.csv', headers, data)

    x = float(input("Координаты точки: "))
    print("Значение функции:", count_f(x))
    if x == 3:
        k = n - 1
    else:
        k = 0
        while True:
            if k < n - 1:
                if X[k] <= x <= X[k+1]:
                    break
            else:
                break
            k += 1
    c_s = count_spline(x, A, B, C, D, X, k)
    print(f"Интерполированное значение: {c_s}")
    print(f"Погрешность: {abs(c_s - count_f(x))}")

    x_vals = [a + (i - 0.5) * h for i in range(n)]
    y_vals = [count_f(x) for x in x_vals]
    spline_vals = [count_spline(x_vals[i], A, B, C, D, X, i)
                   for i in range(len(x_vals))]
    p = [abs(spline_vals[i] - y_vals[i]) for i in range(len(x_vals))]

    additional_headers = ["i", "x_i", "f(x_i)", "S(x_i)", "inaccuracy"]
    additional_data = []
    for i in range(len(x_vals)):
        row = [i, x_vals[i], y_vals[i], spline_vals[i], p[i]]
        additional_data.append(row)

    write_to_csv('additional_output.csv', additional_headers, additional_data)

if __name__ == "__main__":
    main()
