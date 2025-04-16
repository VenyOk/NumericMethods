import math

def runge_kutta_4_system(f1, f2, x0, z1_0, z2_0, x_end, n_steps):
    h = (x_end - x0) / n_steps
    x_values = [x0]
    z1_values = [z1_0]
    z2_values = [z2_0]
    z1 = z1_0
    z2 = z2_0
    x = x0
    for _ in range(n_steps):
        k1_z1 = f1(x, z1, z2)
        k1_z2 = f2(x, z1, z2)
        k2_z1 = f1(x + 0.5 * h, z1 + 0.5 * h * k1_z1, z2 + 0.5 * h * k1_z2)
        k2_z2 = f2(x + 0.5 * h, z1 + 0.5 * h * k1_z1, z2 + 0.5 * h * k1_z2)
        k3_z1 = f1(x + 0.5 * h, z1 + 0.5 * h * k2_z1, z2 + 0.5 * h * k2_z2)
        k3_z2 = f2(x + 0.5 * h, z1 + 0.5 * h * k2_z1, z2 + 0.5 * h * k2_z2)
        k4_z1 = f1(x + h, z1 + h * k3_z1, z2 + h * k3_z2)
        k4_z2 = f2(x + h, z1 + h * k3_z1, z2 + h * k3_z2)
        z1_new = z1 + (h / 6.0) * (k1_z1 + 2 * k2_z1 + 2 * k3_z1 + k4_z1)
        z2_new = z2 + (h / 6.0) * (k1_z2 + 2 * k2_z2 + 2 * k3_z2 + k4_z2)
        x += h
        z1, z2 = z1_new, z2_new
        x_values.append(x)
        z1_values.append(z1)
        z2_values.append(z2)
    return x_values, z1_values, z2_values

def exact_solution(x):
    sqrt6 = math.sqrt(6)
    C1 = 1 / (4 * sqrt6) + 1 / 48
    C2 = 1 / 48 - 1 / (4 * sqrt6)
    return C1 * math.exp((4 + 2 * sqrt6) * x) + C2 * math.exp((4 - 2 * sqrt6) * x) - (1 / 24) * math.exp(4 * x)

def f1(x, z1, z2):
    return z2

def f2(x, z1, z2):
    return math.exp(4 * x) + 8 * z2 + 8 * z1

def solve_with_required_accuracy(eps=0.001):
    x0, y0, yp0 = 0.0, 0.0, 1.0
    X_END = 1.0
    n_steps = 1
    while True:
        x_n, y_n, _ = runge_kutta_4_system(f1, f2, x0, y0, yp0, X_END, n_steps)
        x_2n, y_2n, _ = runge_kutta_4_system(f1, f2, x0, y0, yp0, X_END, 2 * n_steps)
        max_diff = 0.0
        for i in range(n_steps + 1):
            diff = abs(y_n[i] - y_2n[2 * i])
            if diff > max_diff:
                max_diff = diff
        if max_diff < eps:
            return runge_kutta_4_system(f1, f2, x0, y0, yp0, X_END, 2 * n_steps)
        else:
            n_steps += 1

def main():
    x_values, y_values, _ = solve_with_required_accuracy(eps=0.001)
    print("   i     x_i         y_num         y_exact      |diff|")
    for i, (x_i, y_num) in enumerate(zip(x_values, y_values)):
        y_ex = exact_solution(x_i)
        diff = abs(y_num - y_ex)
        print(f"{i:4d}  {x_i:8.5f}  {y_num:12.8f}  {y_ex:12.8f}  {diff:12.8f}")

if __name__ == "__main__":
    main()
