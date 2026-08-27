import numpy as np
from scipy.linalg import expm
import matplotlib.pyplot as plt

ket_0 = np.array([[1], [0]])
ket_1 = np.array([[0], [1]])

X = np.array([[0, 1], [1, 0]])
Y = np.array([[0, -1j], [1j, 0]])
Z = np.array([[1, 0], [0, -1]])
I = np.eye(2)

data = [
    (0.2,  -0.4803,  0.8919, -0.8919,  0.1224,  0.5328),
    (0.4,  -0.8952,  0.6984, -0.6984,  0.0748,  0.4330),
    (0.6,  -1.0513,  0.5214, -0.5214,  0.0339,  0.3109),
    (0.735,-1.0524,  0.3979, -0.3979, -0.0112,  0.1809),
    (0.8,  -1.0400,  0.3603, -0.3603, -0.0329,  0.1509),
    (1.0,  -0.9584,  0.2622, -0.2622, -0.0785,  0.0979),
    (1.2,  -0.8766,  0.1989, -0.1989, -0.1063,  0.0665),
    (1.4,  -0.8062,  0.1551, -0.1551, -0.1229,  0.0469),
    (1.6,  -0.7468,  0.1235, -0.1235, -0.1324,  0.0341),
    (1.8,  -0.6966,  0.0998, -0.0998, -0.1372,  0.0254),
    (2.0,  -0.6539,  0.0816, -0.0816, -0.1391,  0.0192),
    (2.5,  -0.5654,  0.0519, -0.0519, -0.1393,  0.0105),
    (3.0,  -0.5044,  0.0358, -0.0358, -0.1360,  0.0063),
    (4.0,  -0.4312,  0.0198, -0.0198, -0.1287,  0.0025)
]

def Hamiltonian(a0, a1, a2, a3, a4):
    return (a0 * np.kron(I, I) +
            a1 * np.kron(Z, I) +
            a2 * np.kron(I, Z) +
            a3 * np.kron(Z, Z) +
            a4 * (np.kron(Y, Y) + np.kron(X, X)))
            
H2State = np.kron(ket_1, ket_1)
dt = 0.5
t = 10

class QPE():
    def __init__(self, U, t, s_down):
        self.U = U
        self.t = t
        self.s_down = s_down
        self.ket_0 = np.array([[1], [0]])
        self.ket_1 = np.array([[0], [1]])
        self.H = np.array([[1, 1], [1, -1]]) / np.sqrt(2)
        self.I = np.eye(2)
        self.I2 = np.eye(4)
        self.U_swap = np.array([[1, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, 1]])
        self.U_swap3down = np.kron(self.I, self.U_swap) @ np.kron(self.U_swap, self.I)
        self.U_swap3up = np.kron(self.U_swap, self.I) @ np.kron(self.I, self.U_swap)

    #|0>^t
    def state_all(self):
        s_up = self.ket_0
        for i in range(self.t-1):
            s_up = np.kron(s_up, self.ket_0)
        return np.kron(s_up, self.s_down)

    #H^t
    def Hadamar_transform(self, state):
        H_all = self.H
        for i in range(self.t-1):
            H_all = np.kron(H_all, self.H)
        return np.kron(H_all, self.I2) @ state

    @staticmethod
    def Kron_Power(U, n):
        if n == 0:
            return 1
        result = U
        for _ in range(n-1):
            result = np.kron(result, U)
        return result

    @staticmethod
    def Matrix_Power(U, n):
        if n == 0:
            return 1
        result = U
        for _ in range(n-1):
            result = result @ U
        return result

    def Control_Operator(self, U, state):
        pr0 = self.ket_0 @ self.ket_0.conj().T
        pr1 = self.ket_1 @ self.ket_1.conj().T
        mas_C_U = []
        for j in range(self.t):
            C_U = np.kron(pr0, self.I2) + np.kron(pr1, self.Matrix_Power(U, 2**j))
            mas_C_U.append(C_U)
        
        for i in range(self.t):
            state = np.kron(self.Kron_Power(self.I, self.t-1), mas_C_U[i]) @ state
            for j in range(self.t-1):
                state = np.kron(np.kron(self.Kron_Power(self.I, self.t-2-j),
                                        np.kron(self.U_swap, self.Kron_Power(self.I, j))), self.I2) @ state
        return state

    def QFT_dag(self, state):
        QFT = []
        N = 2**self.t
        for i in range(N):
            row = []
            for j in range(N):
                row.append(np.exp((1j*2*np.pi*i*j)/N) / np.sqrt(N))
            QFT.append(row)
        return np.kron(np.array(QFT).conj().T, self.I2) @ state

    def measure(self, state_end):
        Pr_max = 0
        c = 0
        for i in range(2**self.t):
            mes = format(i, f'0{self.t}b')
            mas = [self.ket_0 if x == "0" else self.ket_1 for x in mes]
            state = mas[0]
            for j in range(self.t-1):
                state = np.kron(state, mas[j+1])
            Pr = 0
            for k in range(4):
                mes_k = format(k, f'02b')
                mas_k = [self.ket_0 if x == "0" else self.ket_1 for x in mes_k]
                state_down = np.kron(mas_k[0], mas_k[1])
                Pr += np.abs((np.kron(state, state_down).T @ state_end))**2
            if Pr > Pr_max:
                Pr_max = Pr
                c = i
        return c

    def qpe(self):
        state = self.state_all()
        state = self.Hadamar_transform(state)
        state = self.Control_Operator(self.U, state)
        state = self.QFT_dag(state)
        measured = self.measure(state)
        return measured / (2**self.t)

def estimate_energy(idx, dt, t, state):
    H = Hamiltonian(data[idx][1], data[idx][2], data[idx][3], data[idx][4], data[idx][5])
    U = expm(-1j * H * dt)
    qpe_obj = QPE(U, t, state)
    phase = qpe_obj.qpe()
    energy = -2.0 * np.pi * phase / dt
    return energy

bond_lengths = [data[i][0] for i in range(len(data))]
energy = [estimate_energy(j, dt, t, H2State) for j in range(len(data))]

plt.figure(figsize=(10,6))
plt.plot(bond_lengths, energy, "o-")
plt.xlabel("Длина связи(Å)")
plt.ylabel("Энергия основного состояния(Хартри)")
plt.title("Кривая энергии молекулярного водорода")
plt.grid(True)
plt.show()
