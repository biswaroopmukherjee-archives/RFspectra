import numpy as np
import matplotlib.pyplot as plt


class VirialUnitarity:
    def __init__(self, LogPoints=1e4, Order=3, ContactOrder=3, BetaMuRange=[-6,-0.3]):

        # Input Parameters
        self.LogPoints = LogPoints
        self.Order = Order
        self.ContactOrder = ContactOrder

        # Coefficients for the virial expansion
        b2ref = 3*np.sqrt(2)/8
        b3ref = -0.29095295; # Phys. Rev. Lett. 102, 160401 (2009)
        b4ref = 0.065; # Science 335, 563 (2012)

        # Coefficients for the contact virial expansion
        c2ref = 1/np.pi;
        c3ref = -0.141;

        # Set the orders
        if self.Order==1:
            self.b2 = 0
            self.b3 = 0
            self.b4 = 0
        elif self.Order==2:
            self.b2 = b2ref
            self.b3 = 0
            self.b4 = 0
        elif self.Order==3:
            self.b2 = b2ref
            self.b3 = b3ref
            self.b4 = 0
        elif self.Order==4:
            self.b2 = b2ref
            self.b3 = b3ref
            self.b4 = b4ref

        if self.ContactOrder==2:
            self.c2 = c2ref
            self.c3 = 0
        elif self.ContactOrder==3:
            self.c2 = c2ref
            self.c3 = c3ref

        # Generate EOS data
        BetaMu_start = BetaMuRange[0];
        BetaMu_stop = BetaMuRange[1];
        BetaMu_vec = np.linspace(BetaMu_start, BetaMu_stop, LogPoints);
        Z_vec = np.exp(BetaMu_vec);

        self.PTilde = 10*np.pi/(6*np.pi**2)**(2/3)* (Z_vec + self.b2 * Z_vec**2 + self.b3 * Z_vec**3 + self.b4 * Z_vec**4)/ (Z_vec + 2*self.b2 * Z_vec**2 + 3*self.b3 * Z_vec**3 + 4*self.b4 * Z_vec**4)**(5/3);

        self.KappaTilde = (6*np.pi**2)**(2/3)/(6*np.pi) * (Z_vec + 4*self.b2 * Z_vec**2 + 9*self.b3 * Z_vec**3 + 16*self.b4 * Z_vec**4)/(Z_vec + 2*self.b2 * Z_vec**2 + 3*self.b3 * Z_vec**3 + 4*self.b4 * Z_vec**4)**(1/3);

        self.TTilde = 4*np.pi / (6* np.pi**2 * (Z_vec + 2*self.b2 * Z_vec**2 + 3*self.b3 * Z_vec**3 + 4*self.b4 * Z_vec**4))**(2/3);

        self.CI_NkF = 48 * np.pi**4 * (self.c2 * Z_vec**2+ self.c3 * Z_vec**3)/ (6*np.pi**2*(Z_vec + 2*self.b2 * Z_vec**2 + 3*self.b3 * Z_vec**3 + 4*self.b4 * Z_vec**4))**(4/3);

    def plotCtilde(self, xlim=[0,2]):
        plt.plot(self.TTilde, self.CI_NkF, '-')
        plt.xlim(xlim)
        plt.show()
