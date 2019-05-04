from sympy import *
from numpy import size
from matplotlib import pyplot as plp
import math
import pandas as pd


init_printing(use_unicode=True)


class Tabela_Butchera:
    def __init__(self, Wektor_b, Wektor_c, Macierz_a):
        self.Wektor_b = Wektor_b
        self.Wektor_c = Wektor_c
        self.Macierz_a = Macierz_a

    def getWektor_b(self):
        return self.Wektor_b

    def getWektor_c(self):
        return self.Wektor_c

    def getMacierz_a(self):
        return self.Macierz_a

    def condition1(self):
        c = self.Wektor_c
        a = self.Macierz_a
        for i in range(1, len(c)):
            if (sum([a[i - 1, j - 1] for j in range(1, len(c))]) != c[i]):
                raise ValueError("Błędna Tabela Butchera")

    def condition2(self):
        c = self.Wektor_c
        b = self.Wektor_b
        if (sum([b[i] * ((c[i]) ** (len(c) - 1)) for i in range(0, len(c))]) != 1 / len(c)):
            raise ValueError("Błędna Tabela Butchera")

    def condition3(self):
        c = self.Wektor_c
        a = self.Macierz_a
        for i in range(1, len(c)):
            print(sum([a[i - 1, j] * ((c[j + 1]) ** (len(c) - 1)) for j in range(0, size(a, 0))]),
                  (c[i] ** (len(c))) / len(c))
            # if (sum([a[i, j] * ((c[j]) ** (len(c) - 1)) for j in range(0, size(a, 0))]) != (c[i] ** (len(c))) / len(c)):
            # raise ValueError("Błędna Tabela Butchera")
            print('oo')

    def condition4(self):
        c = self.Wektor_c
        a = self.Macierz_a
        b = self.Wektor_b
        for j in range(0, len(c)):
            if (sum([b[i] * ((c[i]) ** (len(c) - 1)) * a[i - 1, j - 1] for i in range(1, len(c))]) != b[j] * (
                    1 - (c[j] ** (len(c))))):
                raise ValueError("Błędna Tabela Butchera")


class Runge_Kutta:
    def __init__(self, Rownanie_1, Rownanie_2, Symbole_1, Symbole_2, Tabela_Butchera):
        self.Rownanie_1 = Rownanie_1
        self.Rownanie_2 = Rownanie_2
        self.Symbole_1 = Symbole_1
        self.Symbole_2 = Symbole_2
        self.IloscRownan = len(Rownanie_1)
        self.Tabela_Butchera = Tabela_Butchera

    def calculate(self, krok, t0, Wartosc_Poczatkowa_1, Wartosc_Poczatkowa_2, tn):
        n = len(self.Tabela_Butchera.getWektor_c())
        c = self.Tabela_Butchera.getWektor_c()
        a = self.Tabela_Butchera.getMacierz_a()
        b = self.Tabela_Butchera.getWektor_b()
        k = [[i for i in range(len(c))] for i in range(self.IloscRownan)]
        l = [[i for i in range(len(c))] for i in range(self.IloscRownan)]
        i1 = [(self.Symbole_1[i], Wartosc_Poczatkowa_1[i]) for i in range(0, len(self.Symbole_1))]
        i2 = [(self.Symbole_2[j], Wartosc_Poczatkowa_2[j]) for j in range(0, len(self.Symbole_2))]
        print(i1, i2)
        yield t0, tuple(Wartosc_Poczatkowa_1), tuple(Wartosc_Poczatkowa_2)
        while (tn > t0):

            for e1 in range(0, len(i1)):
                k[e1][0] = self.Rownanie_1[e1].subs(i1 + i2)
            for e2 in range(0, len(i2)):
                l[e2][0] = self.Rownanie_2[e2].subs(i1 + i2)

            for i in range(1, n):

                for p in range(0, self.IloscRownan):
                    i1 = [(self.Symbole_1[v], Wartosc_Poczatkowa_1[v] + sum([krok * a[i - 1, j] * k[v][j] for j in range(0, i)])) for v in
                          range(0, len(self.Symbole_1))]
                    i2 = [(self.Symbole_2[v],
                           Wartosc_Poczatkowa_2[v] + sum([krok * a[i - 1, j] * l[v][j] for j in range(0, i)])) for v in
                          range(0, len(self.Symbole_2))]

                    k[p][i] = (self.Rownanie_1[p].subs(i1 + i2))

                for p in range(0, self.IloscRownan):
                    i1 = [(self.Symbole_1[v], Wartosc_Poczatkowa_1[v] + sum([krok * a[i - 1, j] * k[v][j] for j in range(0, i)])) for v in
                          range(0, len(self.Symbole_1))]
                    i2 = [(self.Symbole_2[v],
                           Wartosc_Poczatkowa_2[v] + sum([krok * a[i - 1, j] * l[v][j] for j in range(0, i)])) for v in
                          range(0, len(self.Symbole_2))]

                    l[p][i] = (self.Rownanie_2[p].subs(i1 + i2))

            for p in range(0, len(Wartosc_Poczatkowa_1)):
                Wartosc_Poczatkowa_1[p] = Wartosc_Poczatkowa_1[p] + krok * sum([b[i] * k[p][i] for i in range(0, n)])
            for p in range(0, len(Wartosc_Poczatkowa_2)):
                Wartosc_Poczatkowa_2[p] = Wartosc_Poczatkowa_2[p] + krok * sum([b[i] * l[p][i] for i in range(0, n)])
            t0 = t0 + krok
            yield t0, tuple(Wartosc_Poczatkowa_1),tuple(Wartosc_Poczatkowa_2)


if __name__ == "__main__":
    t = Symbol('t')
    x1 = Symbol('x1')
    y1 = Symbol('y1')
    x2 = Symbol('x2')
    y2 = Symbol('y2')

    vx1 = Symbol('vx1')
    vy1 = Symbol('vy1')
    vx2 = Symbol('vx2')
    vy2 = Symbol('vy2')

    #########################################################################################
    promien_Ziemi = 6371_000
    Promien_Ksiezyca = 384_400_000
    G = 6.673 * (10 ** (-11))
    Masa_Ziemi = 5.972 * (10 ** 24)
    Masa_Ksiezyca = 7348 * 10**22
    Masa_Apollo = 5000
    Promien_Ziemi = promien_Ziemi + 300_000
    GM = 398600.435436 * 1e9 #Standardowy wspolczynnik grawitacyjny Ziemi
    GM1 = 4902.800066 * 1e9 #Standardowy wspolczynnik grawitacyjny Ksiezyca
    alfa = -pi / 4.8
    Sinalfa = math.sin(alfa)
    Cosalfa = math.cos(alfa)
    v = 10849.505
    #########################################################################################

    Symbole1 = [x1, y1, x2, y2]
    Symbole2 = [vx1, vy1, vx2, vy2]

    Rownanie1 = [vx1, vy1, vx2, vy2]
    Rownanie2 = [((1 / (((x1 ** 2 + y1 ** 2)) ** (3 / 2))) * (- x1)) * (GM) ,
                 (1 / (((x1 ** 2 + y1 ** 2)) ** (3 / 2))) * (- y1) *(GM) ,
                 ((1 / (((x2 ** 2 + y2 ** 2)) ** (3 / 2))) * (- x2)) * (GM) + ((1 / ((((x2 - x1) ** 2 + (y2-y1) ** 2)) ** (3 / 2))) * (x1-x2)) * (GM1),
                 (1 / (((x2 ** 2 + y2 ** 2)) ** (3 / 2))) * (- y2) * (GM) + ((1 / ((((x2 - x1) ** 2 + (y2-y1) ** 2)) ** (3 / 2))) * (y1-y2)) * (GM1)]
    #########################################################################################
    b = [(1 / 6), (1 / 3), (1 / 3), (1 / 6)]
    c = [0, (1 / 2), (1 / 2), 1]
    a = Matrix([[(1 / 2), 0, 0], [0, (1 / 2), 0], [0, 0, 1]])
    Tabelka_Butchera_RKIV = Tabela_Butchera(b, c, a)
    Tabelka_Butchera_RKIV.condition1()
    Tabelka_Butchera_RKIV.condition2()
    Tabelka_Butchera_RKIV.condition3()

    Runge_KuttaIV = Runge_Kutta(Rownanie1, Rownanie2, Symbole1, Symbole2,
                                  Tabelka_Butchera_RKIV)

    Wartosc_Poczatkowa_1 = [0,Promien_Ksiezyca,Promien_Ziemi*Cosalfa,Promien_Ziemi*Sinalfa]
    Wartosc_Poczatkowa_2 = [-1019.4,0, -v*Sinalfa,v*Cosalfa]
    lst = list(Runge_KuttaIV.calculate(1, 0, Wartosc_Poczatkowa_1, Wartosc_Poczatkowa_2, 36_000))

    t0,Wartosc_Poczatkowa_1,Wartosc_Poczatkowa_2 = lst[-1]
    lst2 = list(Runge_KuttaIV.calculate(10, t0, list(Wartosc_Poczatkowa_1), list(Wartosc_Poczatkowa_2), 200_000))

    t0,zeros1,zeros2 = lst2[-1]
    lst3 = list(Runge_KuttaIV.calculate(1, t0, list(Wartosc_Poczatkowa_1), list(Wartosc_Poczatkowa_2), 500_000))

    t0, zeros1, zeros2 = lst3[-1]
    lst4 = list(Runge_KuttaIV.calculate(1, t0, list(Wartosc_Poczatkowa_1), list(Wartosc_Poczatkowa_2), 609_785))

    lst5 = lst + lst2 + lst3 + lst4

    Argument = [t for t, _, _ in lst5]
    RK4Method = [krotka for _, krotka, _ in lst5]
    Predkosci = [krotka for _, _, krotka in lst5]

    t0 = Argument
    xs1 = [x[0] for x in RK4Method]
    ys1 = [x[1] for x in RK4Method]
    xs2 = [x[2] for x in RK4Method]
    ys2 = [x[3] for x in RK4Method]
    vxs1 = [x[0] for x in Predkosci]
    vys1 = [x[1] for x in Predkosci]
    vxs2 = [x[2] for x in Predkosci]
    vys2 = [x[3] for x in Predkosci]
    #########################################################################################

    df = pd.DataFrame({'t0': t0,'x(moon)': xs1,'y(moon)': ys1,'vx(moon)':vxs1,'vy(moon)':vys1,'x(apollo)': xs2,'y(apollo)': ys2,'vx(apollo)':vxs2,'vy(apollo)':vys2})
    writer = pd.ExcelWriter(r"D:\dane.xlsx", engine="xlsxwriter")
    df.to_excel(writer, sheet_name='Sheet1')
    writer.save()
    #########################################################################################

    circle = plp.Circle((0, 0), promien_Ziemi, color='b')
    fig, ax = plp.subplots()
    ax.add_artist(circle)
    RK1,= plp.plot(xs1, ys1)
    RK2,= plp.plot(xs2, ys2)
    plp.axis("equal")
    plp.show()
