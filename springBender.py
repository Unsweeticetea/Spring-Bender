'''Runs a DIY desktop continuous spring bender'''
import sys
#import pint
from sympy.physics.units import *
from sympy import pi, sqrt, N, Eq, symbols, solve, simplify
from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QApplication, QMainWindow, QWidget, QHeaderView
from PyQt5.QtGui import QPalette, QColor, QStandardItemModel, QStandardItem
from ui.ui_MainWindow import Ui_MainWindow
from ui.ui_OutputTable import  Ui_OutputTable

lbf=Quantity('lbf', abbrev='lbf')
lbf.set_global_relative_scale_factor((convert_to(pound*acceleration_due_to_gravity,newton))/newton, newton)
gf=Quantity('gf', abbrev='gf')
gf.set_global_relative_scale_factor((convert_to(gram*acceleration_due_to_gravity,newton))/newton, newton)

class MainWindow(QMainWindow, Ui_MainWindow):
    '''main window'''
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setupUi(self)
        self.setupUnits()
        self.connections()

    def connections(self):
        '''Programs all connections for main window'''
        #self.btn_Generate.clicked.connect(self.equations)
        self.btn_Generate.clicked.connect(self.test)

    def setupUnits(self):
        '''adds units to combo boxes'''
        self.cb_F.addItems(["lbf","gf"])
        self.cb_Y.addItems(["in","mm"])
        self.cb_L.addItems(["in","mm"])
        self.cb_d.addItems(["in","mm"])
        self.cb_OD.addItems(["in","mm"])
        self.cb_E.addItems(["Mpsi","GPa"])
        self.cb_G.addItems(["Mpsi","GPa"])
        self.cb_A.addItems(["kpsi*in^m","MPa*mm^m"])

    def equations(self):
        sibase = UnitSystem.get_unit_system("SI")._base_units
        designZeros = 0
        materialZeros = 0

        F, Y, L, E, G, A, m, ns, xi, d, OD = symbols('F Y L E G A m ns xi d OD')
        ssy, alpha, beta, C, D, kb, tau, ID, na, nt = symbols('ssy alpha beta C D kb tau ID na nt')
        
        #test and assign design properties
        if self.sb_F.value() != 0:
            F= self.sb_F.value()
            F*= lbf if self.cb_F.currentIndex()==0 else gf
        else:
            designZeros+=1
        if self.sb_Y.value() != 0:
            Y= self.sb_Y.value()
            Y*= inch if self.cb_Y.currentIndex()==0 else mm
        else:
            designZeros+=1
        if self.sb_L.value() != 0:
            L= self.sb_L.value()
            L*= inch if self.cb_L.currentIndex()==0 else mm
        else:
            designZeros+=1
        if self.sb_d.value() != 0:
            d= self.sb_d.value()
            d*= inch if self.cb_d.currentIndex()==0 else mm
        else:
            designZeros+=1
        if self.sb_OD.value() != 0:
            OD= self.sb_OD.value()
            OD*= inch if self.cb_OD.currentIndex()==0 else mm
        else:
            designZeros+=1
        
        #test and assign material properties
        if self.sb_E.value() != 0:
            E= self.sb_E.value()
            E*= mega*psi if self.cb_E.currentIndex()==0 else giga*Pa
        else:
            materialZeros+=1
        if self.sb_G.value() != 0:
            G= self.sb_E.value()
            G*= mega*psi if self.cb_G.currentIndex()==0 else giga*Pa
        else:
            materialZeros+=1
        if self.sb_m.value() != 0:
            m= self.sb_m.value()
        else:
            materialZeros+=1
        if self.sb_A.value() != 0:
            A= self.sb_A.value()
            A*= kilo*psi*inch**m if self.cb_A.currentIndex()==0 else mega*Pa*mm**m
        else:
            materialZeros+=1
        if self.sb_ns.value() != 0:
            ns= self.sb_ns.value()
        else:
            materialZeros+=1
        if self.sb_Xi.value() != 0:
            xi= self.sb_Xi.value()
        else:
            materialZeros+=1

        #tests if proper number of zero values are present
        if designZeros!=1 or materialZeros!=0:
            #print("invalid inputs.\nDesign:",designZeros,"\nMaterial:",materialZeros)
            return
        else:
            pass
            #print("valid inputs.\nDesign:",designZeros,"\nMaterial:",materialZeros)

        eqSsy=Eq(ssy,convert_to(.45*A/(d**m),sibase))
        ssy=solve(eqSsy, eqSsy.lhs)[0]
        print("ssy:",ssy)
        eqAlpha=Eq(alpha,convert_to(eqSsy.rhs/ns,sibase))
        alpha=solve(eqAlpha, eqAlpha.lhs)[0]
        print("alpha:",alpha)
        eqBeta=Eq(beta,convert_to((8*(1+xi)*F)/(N(pi)*(d**2)),sibase))
        beta=solve(eqBeta, eqBeta.lhs)[0]
        print("beta:",beta)
        eqC=Eq(C,convert_to(((2*eqAlpha.rhs-eqBeta.rhs)/(4*eqBeta.rhs))+sqrt((((2*eqAlpha.rhs-eqBeta.rhs)/(4*eqBeta.rhs))**2)-(3*eqAlpha.rhs)/(4*eqBeta.rhs)),sibase))
        C=solve(eqC,eqC.lhs)[0]
        print("C:",C)
        
        '''
        D=C*d
        kb=(4*C+2)/(4*C-3)
        tau=kb*8*F*(1+xi)*D/(N(pi)*d**3)
        OD=D+d
        ID=D-d
        na=(G*(d**4)*Y)/(8*(D**3)*F)
        nt=na+2'''
        #print(nsolve([ssy, alpha, beta, C, D, kb, tau, OD, ID, na, nt],[F, Y, L, E, G, A, m, ns, xi, d, OD]))
        #print(alpha,'\n', beta,'\n', C,'\n', D,'\n', kb,'\n', tau,'\n', OD,'\n', ID,'\n', na,'\n', nt)

    def test(self):
        self.sb_F.setValue(20)
        self.cb_F.setCurrentIndex(0)
        self.sb_Y.setValue(2)
        self.cb_Y.setCurrentIndex(0)
        self.sb_L.setValue(3.25)
        self.cb_L.setCurrentIndex(0)
        #determine OD
        if True:
            self.sb_d.setValue(.08)
            self.cb_d.setCurrentIndex(0)
            self.sb_OD.setValue(0)
            self.cb_OD.setCurrentIndex(0)
        else:
            self.sb_d.setValue(0)
            self.cb_d.setCurrentIndex(0)
            self.sb_OD.setValue(0.922679008)
            self.cb_OD.setCurrentIndex(0)
        #determine D
        self.equations()

if __name__ == "__main__":
    app = QApplication(sys.argv)
    win = MainWindow()
    app.setStyle("Fusion")
    palette = QPalette()
    palette.setColor(QPalette.Window, QColor(53, 53, 53))
    palette.setColor(QPalette.WindowText, Qt.white)
    palette.setColor(QPalette.Base, QColor(25, 25, 25))
    palette.setColor(QPalette.AlternateBase, QColor(53, 53, 53))
    palette.setColor(QPalette.ToolTipBase, Qt.black)
    palette.setColor(QPalette.ToolTipText, Qt.white)
    palette.setColor(QPalette.Text, Qt.white)
    palette.setColor(QPalette.Button, QColor(53, 53, 53))
    palette.setColor(QPalette.ButtonText, Qt.white)
    palette.setColor(QPalette.BrightText, Qt.red)
    palette.setColor(QPalette.Link, QColor(42, 130, 218))
    palette.setColor(QPalette.Highlight, QColor(42, 130, 218))
    palette.setColor(QPalette.HighlightedText, Qt.black)
    app.setPalette(palette)
    win.show()
    sys.exit(app.exec())
