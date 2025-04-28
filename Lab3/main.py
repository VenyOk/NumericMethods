import math
from typing import Callable, List, Tuple


def differential_equation(x: float, Y: List[float]) -> List[float]:
    return [Y[1], math.exp(4 * x) + 8 * Y[1] + 8 * Y[0]]


def runge_kutta_step(
    func: Callable[[float, List[float]], List[float]], 
    x: float, 
    Y: List[float], 
    step_size: float
) -> List[float]:
    k1 = [step_size * val for val in func(x, Y)]
    Y1 = [Y[i] + 0.5 * k1[i] for i in range(len(Y))]

    k2 = [step_size * val for val in func(x + 0.5 * step_size, Y1)]
    Y2 = [Y[i] + 0.5 * k2[i] for i in range(len(Y))]

    k3 = [step_size * val for val in func(x + 0.5 * step_size, Y2)]
    Y3 = [Y[i] + k3[i] for i in range(len(Y))]

    k4 = [step_size * val for val in func(x + step_size, Y3)]

    Y_next = [
        Y[i] + (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) / 6 
        for i in range(len(Y))
    ]
    return Y_next


def exact_solution(x: float) -> Tuple[float, float]:
    sqrt6 = math.sqrt(6)
    y = ((2*sqrt6 + 1)/48 * math.exp((4 + 2*sqrt6)*x) + 
         (1 - 2*sqrt6)/48 * math.exp((4 - 2*sqrt6)*x) + 
         (-1/24 * math.exp(4*x)))
    
    coef1 = (5*sqrt6 + 14)/24
    coef2 = (14 - 5*sqrt6)/24
    dy = (coef1 * math.exp((4 + 2*sqrt6)*x) + 
          coef2 * math.exp((4 - 2*sqrt6)*x) + 
          -1/6 * math.exp(4*x))
    return y, dy


def runge_kutta_solver(
    func: Callable[[float, List[float]], List[float]], 
    x0: float, 
    Y0: List[float], 
    x_end: float, 
    step_size: float
) -> List[Tuple[float, float, float, float, float, float]]:
    results = []
    x = x0
    Y = Y0.copy()

    exact_y, exact_y_prime = exact_solution(x)
    results.append((x, Y[0], Y[1], exact_y, exact_y_prime, 0.0))

    while x < x_end:
        if x + 2*step_size > x_end:
            return results

        Y_h_prev = runge_kutta_step(func, x, Y, step_size)
        Y_h = runge_kutta_step(func, x + step_size, Y_h_prev, step_size)

        Y_2h = runge_kutta_step(func, x, Y, 2*step_size)

        error = abs(Y_h[0] - Y_2h[0]) / (2 ** 4 - 1)
        x += 2*step_size
        Y = Y_h

        exact_y, exact_y_prime = exact_solution(x)
        results.append((x, Y[0], Y[1], exact_y, exact_y_prime, error))

    return results


def print_results(results: List[Tuple[float, float, float, float, float, float]]) -> None:
    headers = ("x", "Approx y", "Exact y", "Approx y'", "Exact y'", "Error")
    print("{:>8} {:>14} {:>14} {:>14} {:>14} {:>14}".format(*headers))
    
    for x_val, approx_y, approx_yprime, exact_y, exact_dy, err in results:
        print("{:8.4f} {:14.6f} {:14.6f} {:14.6f} {:14.6f} {:14.6f}".format(
            x_val, approx_y, exact_y, approx_yprime, exact_dy, err))


def main():
    initial_x = 0.0
    final_x = 1.0
    initial_conditions = [0.0, 1.0]
    step_size = 0.01

    solution = runge_kutta_solver(
        differential_equation, 
        initial_x, 
        initial_conditions, 
        final_x, 
        step_size
    )

    print_results(solution)


if __name__ == "__main__":
    main()