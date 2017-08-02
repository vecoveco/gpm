
def FEstat():
    """

    Statistik zum Fererkundungs Tutorium

    """
    import numpy as np
    import matplotlib.pyplot as plt

    ff = 20

    Uebung = np.arange(1,7,1)
    Quickies = np.array([3,5,6,5,8,5])
    Aufgaben = np.array([2,3,4,2,4,8])
    Program = np.array([1,1,1,1,0,0])

    plt.plot(Uebung, Quickies+Aufgaben+Program,
             label='Anzahl der Aufgaben pro Uebung',
             lw=3)
    plt.plot(Uebung, Quickies, label='Quickies', lw=3)
    plt.plot(Uebung, Aufgaben, label='Aufgaben', lw=3)
    plt.plot(Uebung, Program, label='Programiereaufgabe',lw=3)

    plt.xlabel('Nummer der Uebung', fontsize=ff)
    plt.ylabel('Anzahl #', fontsize=ff)
    plt.legend(loc=2)
    plt.title('Verlauf der Uebungen in FE', fontsize=ff)
    plt.xticks(fontsize=ff)
    plt.yticks(fontsize=ff)

    plt.grid(lw=2, color='black')
    plt.show()

FEstat()