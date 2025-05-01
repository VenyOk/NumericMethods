import math
import numpy as np

def f(x, y):
    return np.array([
        y[1],                 
        math.exp(4 * x)  + 8 * y[0] + 8*y[1]
    ])

def f2(x):
    return -2 * math.sin(2*x) + (2.0/3.0) * math.cos(2*x) + (1.0/3.0) * math.cos(x)

def f1(x):
    sqrt6 = math.sqrt(6)
    y = ((2*sqrt6 + 1)/48 * math.exp((4 + 2*sqrt6)*x) +
    (1 - 2*sqrt6)/48 * math.exp((4 - 2*sqrt6)*x) +
    (-1/24 * math.exp(4*x)))
    return y


def runge_kutta_step(f, x, y, h):
    k1 = f(x, y)
    
    k2 = f(x + h/2, np.array([
        y[0] + h * k1[0]/2,
        y[1] + h * k1[1]/2
    ]))
    
    k3 = f(x + h/2, np.array([
        y[0] + h * k2[0]/2,
        y[1] + h * k2[1]/2
    ]))
    
    k4 = f(x + h, np.array([
        y[0] + h * k3[0],
        y[1] + h * k3[1]
    ]))
    
    return np.array([
        y[0] + (k1[0] + 2*k2[0] + 2*k3[0] + k4[0]) * h/6,
        y[1] + (k1[1] + 2*k2[1] + 2*k3[1] + k4[1]) * h/6
    ])

def runge_kutta(f, y0, x0, xend, h):
    result = []
    x = x0
    y = y0.copy()
    
    while x <= xend:
        result.append(y.copy())
        y = runge_kutta_step(f, x, y, h)
        x += h
    
    return result

def norm_full(a, b):
    max_val = 0.0
    for i in range(len(a)):
        if 2*i >= len(b):
            continue
        for j in range(1, len(a[i])):
            cur_val = abs(a[i][j] - b[2*i][j])
            if cur_val > max_val:
                max_val = cur_val
    return max_val

def main():
    x0 = 0.0
    xend = 1.0
    h = 1.0
    
    eps = 0.001
    
    y0 = np.array([0.0, 1.0])
    cnt = 0
    while True:
        sol_h = runge_kutta(f, y0, x0, xend, h)
        sol_h2 = runge_kutta(f, y0, x0, xend, h/2)
        
        diff = norm_full(sol_h, sol_h2)
        
        if diff < eps:
            print(diff / (2 **4 - 1))
            break
        cnt += 1
        h /= 2
    print(cnt)
    print(f"h = {h}")
    
    solution = runge_kutta(f, y0, x0, xend, h)
    
    print("x\t\t y(x)\t\ty(x) (точное)\t\tпогрешность")
    x = x0
    for i in range(len(solution)):
        print(f"{x:.2f}\t\t{solution[i][0]:.6f}\t{f1(x):.6f}\t{abs(f1(x)-solution[i][0]):.6f}")
        x += h
    

if __name__ == "__main__":
    main()