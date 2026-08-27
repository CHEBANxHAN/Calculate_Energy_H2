import qsharp
import matplotlib.pyplot as plt

qsharp.init()

with open("Calc_Energy.qs", "r", encoding="utf-8") as f:
    qsharp.eval(f.read())

bond_lengths = [qsharp.eval(f"Calc_Energy.H2Coeffs({i})[0]") for i in range(14)]
energy = [qsharp.eval(f"Calc_Energy.Calc({j})") for j in range(14)]

plt.figure(figsize=(10,6))
plt.plot(bond_lengths, energy, "o-")
plt.xlabel("Длина связи(Å)")
plt.ylabel("Энергия основного состояния(Хартри)")
plt.title("Кривая энергии молекулярного водорода")
plt.grid(True)
plt.show()
