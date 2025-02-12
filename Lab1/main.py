import copy


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


def compute_residual(A, B, X):
    n = len(A)
    residual = [abs(B[i] - sum(A[i][j] * X[j] for j in range(n))) for i in range(n)]
    max_residual = max(residual)
    max_x = max(abs(x) for x in X)
    relative_residual = max_residual / max_x if max_x != 0 else float('inf')
    return max_residual, relative_residual


def main():
    k = 0
    alpha = beta = 0.01 * k
    A = [[2.0 + alpha, 1.0, -0.1, 1.0],
         [0.4, 0.5 - alpha, 4.0, -8.5],
         [0.3, -1.0, 1.0 + 2 * alpha, 5.2],
         [1.0, 0.2, 2.5, -1.0 + alpha]]
    B = [2.7 - beta, 21.9, -3.9 + beta, 9.9]
    A_copy = copy.deepcopy(A)
    B_copy = copy.deepcopy(B)
    X = solve(A_copy, B_copy)
    print("Решение системы:")
    for i in range(4):
        print(f"X[{i}] = {X[i]}")
    max_residual, relative_residual = compute_residual(A, B, X)
    print(f"Невязка: {max_residual}")
    print(f"Относительная невязка: {relative_residual}")


if __name__ == "__main__":
    main()
