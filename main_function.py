import numpy as np
import matplotlib.pyplot as plt

# Stałe
G = 6.6743e-11


def e(i, x, h):
    if (i - 1) * h < x < i * h:
        return x / h - i + 1
    elif i * h < x < (i + 1) * h:
        return i + 1 - x / h
    else:
        return 0


def de(i, x, h):
    if (i - 1) * h < x < i * h:
        return 1 / h
    elif i * h < x < (i + 1) * h:
        return -1 / h
    else:
        return 0


def B(i, j, h):
    if abs(j - i) > 1:
        return 0
    elif j == i:
        a = (i - 1) * h
        b = (i + 1) * h
    else:
        a = min(i, j) * h
        b = max(i, j) * h

    midpoint1 = (b - a) / (2 * np.sqrt(3)) + (a + b) / 2
    midpoint2 = (a - b) / (2 * np.sqrt(3)) + (a + b) / 2

    return ((a - b) / 2) * (de(i, midpoint1, h)*de(j, midpoint1, h) + de(i, midpoint2, h)*de(j, midpoint2, h))


def Lt(j,h):
    a = (j - 1) * h
    b = (j + 1) * h
    Bt = 0
    L = 0

    midpoint1 = (b - a) / (2 * np.sqrt(3)) + (a + b) / 2
    midpoint2 = (a - b) / (2 * np.sqrt(3)) + (a + b) / 2

    if j * h < 1 or j * h > 2:
        Bt = G*((a - b) / 2) * (de(j, midpoint1, h) + de(j, midpoint2, h))/3
        return L - Bt
    L = ((a - b) / 2) * (e(j, midpoint1, h) + e(j, midpoint2, h)) *4*np.pi*1e11*G
    return L-Bt


def make_matrix(N):
    M = np.zeros((N - 1, N - 1))
    Y = np.zeros(N - 1)
    for i in range(1, N):
        Y[i - 1] = Lt(i, 3 / N)
        for j in range(1, N):
            M[i - 1, j - 1] = B(i, j, 3/N)
    return M,Y

def Phi(x, W, N):
    suma = 0.0
    for i in range(1, N):
        suma += W[i - 1] * e(i, x, 3/N)
    return -5*G + G*x / 3 + suma


if __name__ == '__main__':
    N = int(input("Wpisz liczbę przedziałów N: "))
    M,Y = make_matrix(N)
    W = np.linalg.solve(M, Y)

    x_vals = np.linspace(0, 3, 1000)
    y_vals = [Phi(x, W, N) for x in x_vals]

    plt.plot(x_vals, y_vals, linewidth=3)
    plt.title("Wykres przybliżonej funkcji Φ dla n = " + str(N))
    plt.xlabel("x")
    plt.ylabel("Φ(x)")
    plt.grid()
    plt.show()
