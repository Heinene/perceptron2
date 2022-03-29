import math
import matplotlib.pyplot as mp

def Func(t):
    x = math.exp(-0.1*t*t)
    return x

class Neiron:
    def __init__(self, a, b, N, n, p):       
        self._dec = []
        self._w = [0.0] * (p + 1)
        self._n = n
        self._po = []
        self._p = p
        self._N = N
        for i in range(2 * N):
            x_i = a + i * (b - a) / N
            self._po.append(x_i)
            self._dec.append(Func(self._po[i]))
        self._pred = self._dec[:self._p] + [0.0] * (len(self._dec) - self._p)
        self._eps = 0

    def Stud(self, M):        
        self._w = [0.0] * (self._p + 1)
        k = 0
        eps = 0.0
        while k < M:
            for index in range(self._p, self._N):
                net = self._w[0]
                for i in range(1, self._p + 1):
                    net += self._w[i] * self._dec[index - self._p + i - 1]
                self._pred[index] = net
                delta = self._dec[index] - self._pred[index]
                self._w[0] += self._n * delta * 1
                for j in range(1, self._p + 1):
                    self._w[j] += self._n * delta * self._dec[index - self._p + j - 1]
            eps = 0.0
            try:
                for i in range(self._p, self._N):
                    eps += math.pow(self._dec[i] - self._pred[i], 2)
            except OverflowError:
                print('Not learned' , self._p)
                return
            eps = math.sqrt(eps)
            if eps < 0.0001:
                break
            k += 1
        print('P = {}, n = {:.2}, eras = {}, epsilon = {:.6}, w = '.format(self._p, self._n, k, eps), end='')
        for i in range(len(self._w) - 1):
            print('{:.4}'.format(self._w[i]), end=', ')
        print('{:.4}'.format(self._w[-1]))
        self._eps = eps

    def Preds(self):       
        x = self._dec[:self._N] + self._pred[self._N:]
        for i in range(self._N, self._N * 2):
            net = self._w[0]
            for j in range(1, self._p + 1):
                net += self._w[j] * x[i - self._p + j - 1]
            x[i] = net
            self._pred[i] = net

    def g_dec(self):
        return self._dec

    def g_po(self):
        return self._po

    def g_preds(self):
        x = self._dec[:self._N] + self._pred[self._N:]
        return x

    def get_eps(self):
        return self._eps

def main():
    nei = Neiron(-5, 5, 20, 0.06, 16)
    nei.Stud(1000)
    nei.Preds()
    predicted = nei.g_preds()
    decision = nei.g_dec()
    x = nei.g_po()
    mp.xlabel('x')
    mp.ylabel('y')
    mp.grid()
    mp.plot(x, decision)
    mp.show()
    mp.xlabel('x')
    mp.ylabel('y')
    mp.grid()
    mp.plot(x, decision)
    mp.plot(x, predicted)
    mp.show()
    eps_mas = []
    m_mas = []
    for m in range(1, 20):
        nei1 = Neiron(-5, 5, 20, 0.06, 16)
        nei1.Stud(250 * m)
        nei1.Preds()
        eps_mas.append(nei1.get_eps())
        m_mas.append(250 * m)
    mp.xlabel('M')
    mp.ylabel('eps')
    mp.grid()
    mp.plot(m_mas, eps_mas)
    mp.show()
    eps_mas = []
    n_array = []
    for n in range(1, 21):
        nei2 = Neiron(-5, 5, 20, 0.04 * n, 4)
        nei2.Stud(5000)
        nei2.Preds()
        eps_mas.append(nei2.get_eps())
        n_array.append(0.05 * n)
    mp.xlabel('n')
    mp.ylabel('eps')
    mp.grid()
    mp.plot(n_array, eps_mas)
    mp.show()
    eps_mas = []
    p_array = []
    for p in range(2, 20):
        nei3 = Neiron(-5, 5, 20, 0.1, p)
        nei3.Stud(5000)
        nei3.Preds()
        eps_mas.append(nei3.get_eps())
        p_array.append(p)
    mp.xlabel('p')
    mp.ylabel('eps')
    mp.grid()
    mp.plot(p_array, eps_mas)
    mp.show()

if __name__ == '__main__':
    main()

