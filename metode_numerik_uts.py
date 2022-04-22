import os
from numpy import array, zeros, fabs, linalg, diag, diagflat, dot
import numpy as np

a = array([[1, 2, 3, -1],
           [2, 5, 4, 8],
           [4, 2, 2, 1],
           [6, 4, -1, -2]], int)
# suku-suku konstanta matriks b dari persamaan
b = array([10, 8, -2, 4], int)

A = array([[1, 2, 3, -1], [2, 5, 4, 8], [4, 2, 2, 1], [6, 4, -1, -2]])
b = array([10, 8, -2, 4])
guess = array([1, 1, 1, 1])


def f(x):
    return x**3 + 2*x**2 + 10*x


def g(x):
    return 3*x**2 + 4*x + 10


def invers():
    # Python program to inverse
    # a matrix using numpy

    # Import required package

    # Taking a 4 * 4 matrix
    A = np.array([[1, 2, 3, -1],
                  [2, 5, 4, 8],
                  [4, 2, 2, 1],
                  [6, 4, -1, -2]])

    # Calculating the inverse of the matrix
    print(np.linalg.inv(A))
    print(" ")


def jacobi(A, b, N=25, x=None):
    """Solves the equation Ax=b via the Jacobi iterative method."""
    # Create an initial guess if needed
    if x is None:
        x = zeros(len(A[0]))

    # Create a vector of the diagonal elements of A
    # and subtract them from A
    D = diag(A)
    R = A - diagflat(D)

    # Iterate for N times
    for i in range(N):
        x = (b - dot(R, x)) / D
    return x


def gauss():
    print("Solusi dari Numpy :")
    print(linalg.solve(a, b))

    n = len(b)
    x = zeros(n, float)

    # loop pertama menentukan baris tetap
    for k in range(n-1):
        if fabs(a[k, k]) < 1.0e-12:

            for i in range(k+1, n):
                if fabs(a[i, k]) > fabs(a[k, k]):
                    a[[k, i]] = a[[i, k]]
                    b[[k, i]] = b[[i, k]]
                    break

    # menerapkan eliminasi di bawah baris tetap

        for i in range(k+1, n):
            if a[i, k] == 0:
                continue

            factor = a[k, k]/a[i, k]
            for j in range(k, n):
                a[i, j] = a[k, j] - a[i, j]*factor
                # kami juga menghitung vektor b dari setiap baris
            b[i] = b[k] - b[i]*factor
    print(a)
    print(b)

    x[n-1] = b[n-1] / a[n-1, n-1]
    for i in range(n-2, -1, -1):
        sum_ax = 0

        for j in range(i+1, n):
            sum_ax += a[i, j] * x[j]

        x[i] = (b[i] - sum_ax) / a[i, i]

    print("\nSolusi dari sistem tersebut adalah :")
    print(x)
    print(" ")


def secant(x0, x1, e, N):
    print('\n\n*** IMPLEMENTASI METODE SECANT ***')
    step = 1
    condition = True
    while condition:
        if f(x0) == f(x1):
            print("ERROR ! PEMBAGIAN DENGAN NOL")
            break

        x2 = x0 - (x1-x0)*f(x0)/(f(x1) - f(x0))
        print('Iterasi ke-%d, x2 = %0.6f dan f(x2) = %0.6f' % (step, x2, f(x2)))
        x0 = x1
        x1 = x2
        step = step + 1

        if step > N:
            print('Tidak Konvergen')
            break

        condition = abs(f(x2)) > e
    print('\nAkar yang diperlukan adalah : %0.8f' % x2)


def newtonRaphson(x0, e, N):
    print("\n\n*** IMPLEMENTASI METODE NEWTON RAPHSON ***")
    step = 1
    flag = 1
    Condition = True
    while Condition:
        if g(x0) == 0.0:
            print("ERROR ! PEMBAGIAN DENGAN NOL")
            break
        x1 = x0 - f(x0) / g(x0)
        print("Iterasi ke-%d, x1 = %0.6f and f(x1) = %0.6f" % (step, x1, f(x1)))
        x0 = x1
        step = step + 1
        if step > N:
            flag = 0
            break

        Condition = abs(f(x1)) > e

    if flag == 1:
        print("\nAkar yang diperlukan adalah : %0.8f" % x1)
    else:
        print("\nTidak Konvergen")


def soal1():
    print("Diketahui : f(x) = x^3 + 2x^2 + 10x = 20")
    print("Gunakan Metode : ")
    print(" ")
    print("(i)   Newton-Raphson")
    print("(ii)  Secant")
    print("(iii) Kembali")
    pilih1 = input("Pilih [i/ii/iii] : ")
    if pilih1 == 'i':
        print(" ")
        x0 = int(input("Masukkan Tebakan (x0) : "))
        e = float(input("Kesalahan di toleransi (e): "))
        N = int(input("Langkah Maksimum (n): "))
        newtonRaphson(x0, e, N)
        print(" ")
        os.system("pause")
        os.system("cls")
        soal1()
    elif pilih1 == 'ii':
        print(" ")
        x0 = int(input('Masukkan tebakan pertama : '))
        x1 = int(input('Masukkan tebakan kedua   : '))
        e = float(input('kesalahan di toleransi   : '))
        N = int(input('Maksimal step            : '))

        # convert type data
        x0 = float(x0)
        x1 = float(x1)
        e = float(e)
        N = int(N)
        secant(x0, x1, e, N)
        print(" ")
        os.system("pause")
        os.system("cls")
        soal1()
    elif pilih1 == 'iii':
        print(" ")
        pass
        os.system("pause")
        os.system("cls")
        utama()
    else:
        print("Anda Memilih Input yang salah !")
        print("      Silahkan Coba Lagi       ")
        print(" ")
        os.system("pause")
        os.system("cls")


def soal2():
    print("Diberikan sistem persamaan linear Ax = b dengan A dan b sebagai berikut : ")
    print(" ")
    print("      | 1   2   3  -1 |              | 10 |")
    print("  A = | 2   5   4   8 |     dan b =  |  8 |")
    print("      | 4   2   2   1 |              | -2 |")
    print("      | 6   4  -1  -2 |              |  4 |")
    print(" ")
    print("a. Tentukan solusi dengan metode eliminasi Gauss")
    print("b. Tentukan solusi dengan metode matriks balikan")
    print("c. Tentukan solusi dengan metode Iterasi Jacobi")
    print("d. Kembali")
    pilih2 = input("Pilih [a/b/c/d] : ")
    if pilih2 == 'a':
        print(" ")
        gauss()
        os.system("pause")
        os.system("cls")
        soal2()
    elif pilih2 == 'b':
        print(" ")
        invers()
        os.system("pause")
        os.system("cls")
        soal2()
    elif pilih2 == 'c':
        print(" ")
        sol = jacobi(A, b, N=25, x=guess)
        print("A:")
        print(A)

        print("b:")
        print(b)

        print("x:")
        print(sol)
        os.system("pause")
        os.system("cls")
        soal2()
    elif pilih2 == 'd':
        print(" ")
        pass
        os.system("pause")
        os.system("cls")
        utama()
    else:
        print("Anda Memilih Input yang salah !")
        print("      Silahkan Coba Lagi       ")
        print(" ")
        os.system("pause")
        os.system("cls")


def utama():
    print("           TUGAS TENGAH SEMESTER          ")
    print("MATAKULIAH ANALISIS STATISTIK DAN  NUMERIK")
    print("      PROGRAM STUDI INFORMATIKA UNMUL     ")
    print("                TAHUN 2002                ")
    print("==========================================")
    print(" ")
    print("a. Soal Nomer 1")
    print("b. Soal Nomer 2")
    print("c. Keluar")
    print(" ")
    pilih = input("pilih [a/b/c] : ")
    if pilih == 'a':
        os.system("cls")
        soal1()
    elif pilih == 'b':
        os.system("cls")
        soal2()
    elif pilih == 'c':
        os.system("cls")
        exit(0)


utama()
