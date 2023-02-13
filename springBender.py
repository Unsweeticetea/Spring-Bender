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

'''MPsi=Quantity('MegaPsi', abbrev='MPsi')
MPsi.set_global_relative_scale_factor(mega, psi)
KPsi=Quantity('MegaPsi', abbrev='KPsi')
KPsi.set_global_relative_scale_factor(kilo, psi)
GPa=Quantity('GigaPa', abbrev='GPa')
GPa.set_global_relative_scale_factor(giga, Pa)
MPa=Quantity('MegaPa', abbrev='MPa')
MPa.set_global_relative_scale_factor(mega, Pa)'''
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
            
        #if not (self.countZeros([F, Y, L, d, OD], 1) and self.countZeros([E, G, m, A, ns, xi], 0)):
        if designZeros!=1 or materialZeros!=0:
            print("invalid inputs")
            return
        #else:
            #print("valid inputs")

        eqSsy=Eq(ssy,.45*A/(d**m))
        ssy=solve(eqSsy)[0]
        eqAlpha=Eq(alpha,eqSsy.rhs/ns)
        alpha=solve(eqAlpha)[0]
        eqBeta=Eq(beta,convert_to((8*(1+xi)*F)/(N(pi)*(d**2)),psi)) #I hate hardcoding this into psi, but I don't know how to get it to automatically convert to apropriate pressure units
        beta=solve(eqBeta)[0]
        eqC=Eq(C,((2*eqAlpha.rhs-eqBeta.rhs)/(4*eqBeta.rhs))+sqrt((((2*eqAlpha.rhs-eqBeta.rhs)/(4*eqBeta.rhs))**2)-(3*eqAlpha.rhs)/(4*eqBeta.rhs)))
        C=solve(eqC)[0]
        print(ssy, alpha, beta, C)
        
        '''
        C=((2*alpha-beta)/(4*beta))+sqrt((((2*alpha-beta)/(4*beta))**2)-(3*alpha)/(4*beta))
        D=C*d
        kb=(4*C+2)/(4*C-3)
        tau=kb*8*F*(1+xi)*D/(N(pi)*d**3)
        OD=D+d
        ID=D-d
        na=(G*(d**4)*Y)/(8*(D**3)*F)
        nt=na+2'''
        #print(nsolve([ssy, alpha, beta, C, D, kb, tau, OD, ID, na, nt],[F, Y, L, E, G, A, m, ns, xi, d, OD]))
        #print(alpha,'\n', beta,'\n', C,'\n', D,'\n', kb,'\n', tau,'\n', OD,'\n', ID,'\n', na,'\n', nt)

    '''def equations(self):
        units = pint.UnitRegistry()

        F, Y, L, E, G, A, m, ns, xi, d, OD = symbols('F Y L E G A m ns xi d OD')
        ssy, alpha, beta, C, D, kb, tau, ID, na, nt = symbols('ssy alpha beta C D kb tau ID na nt')

        if self.sb_F.value() != 0:
            F= self.sb_F.value()
        F*=units.force_pound if self.cb_F.currentIndex()==0 else self.sb_F.value()*units.force_gram
        if self.sb_Y.value() != 0:
            Y= self.sb_Y.value()
        Y*=units.inch if self.cb_Y.currentIndex()==0 else self.sb_Y.value()*units.mm
        if self.sb_L.value() != 0:
            L= self.sb_L.value()
        L*=units.inch if self.cb_L.currentIndex()==0 else self.sb_L.value()*units.mm
        if self.sb_d.value() != 0:
            d= self.sb_d.value()
        d*=units.inch if self.cb_d.currentIndex()==0 else self.sb_d.value()*units.mm
        if self.sb_OD.value() != 0:
            OD= self.sb_OD.value()
        OD*=units.inch if self.cb_OD.currentIndex()==0 else self.sb_OD.value()*units.mm
        E= self.sb_E.value()*units.Mpsi if self.cb_E.currentIndex()==0 else self.sb_E.value()*units.gigapascal
        G= self.sb_G.value()*units.Mpsi if self.cb_G.currentIndex()==0 else self.sb_G.value()*units.gigapascal
        m= self.sb_m.value()
        A= self.sb_A.value()*1000*units.psi*units.inch**m if self.cb_A.currentIndex()==0 else self.sb_A.value()*units.MPa*units.mm**m
        ns= self.sb_ns.value()
        xi= self.sb_Xi.value()
        print(A, d, m)
        
        eq_ssy=Eq(.45*A/(d**m))
        print(eq_ssy)
        ssy=solve(eq_ssy,ssy)
        alpha=ssy/ns
        beta=(8*(1+xi)*F)/(N(pi)*(d**2))
        C=((2*alpha-beta)/(4*beta))+sqrt((((2*alpha-beta)/(4*beta))**2)-(3*alpha)/(4*beta))
        D=C*d
        kb=(4*C+2)/(4*C-3)
        tau=kb*8*F*(1+xi)*D/(N(pi)*d**3)
        OD=D+d
        ID=D-d
        na=(G*(d**4)*Y)/(8*(D**3)*F)
        nt=na+2
        #print(nsolve([ssy, alpha, beta, C, D, kb, tau, OD, ID, na, nt],[F, Y, L, E, G, A, m, ns, xi, d, OD]))
        #print(alpha,'\n', beta,'\n', C,'\n', D,'\n', kb,'\n', tau,'\n', OD,'\n', ID,'\n', na,'\n', nt)
        print(d,'\n', OD)'''

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

    def show_tabs(self):
        '''Adds or removes advanced tabs'''
        if self.chk_advanced.isChecked():
            self.tabs.addTab(self.tab_Material, "Material")
            self.tabs.addTab(self.tab_Properties, "Properties")
        else:
            self.tabs.removeTab(3)
            self.tabs.removeTab(2)

    def send_row(self, d, f_max, y_max, l0, G, E, A, m, ns, xi):
        '''Creates a row from given inputs'''
        d=d*getattr(u, self.unitD.checkedButton().objectName()[5:])
        #computes expressions
        ssy=.45*A/(d**m)
        alpha=ssy/ns
        beta=u.convert_to((8*(1+xi)*f_max)/(N(pi)*(d**2)),u.psi)
        C=u.convert_to(((2*alpha-beta)/(4*beta))+sqrt(((2*alpha-beta)/(4*beta))**2-((3*alpha)/(4*beta))),1)
        D=C*d
        kb=(4*C+2)/(4*C-3)
        tau=u.convert_to(kb*8*f_max*(1+xi)*D/(N(pi)*d**3),u.psi)
        OD=D+d
        ID=D-d
        na=u.convert_to((G*(d**4)*y_max)/(8*(D**3)*f_max),1)
        nt=na+2
        if self.chk_ground.isChecked():
            p=((l0-2*d))/na
            ls=d*nt
        else:
            p=((l0-3*d))/na
            ls=d*(nt+1)
        lcr=u.convert_to((N(pi)*D/.5)*(((2*(E-G))/(2*G+E))**.5),u.inch)

        #creates items for all output variables
        i_d=QStandardItem(str(d))
        i_ssy=QStandardItem(str(ssy))
        i_alpha=QStandardItem(str(alpha))
        i_beta=QStandardItem(str(beta))
        i_C=QStandardItem(str(C))
        i_D=QStandardItem(str(D))
        i_kb=QStandardItem(str(kb))
        i_tau=QStandardItem(str(tau))
        i_OD=QStandardItem(str(OD))
        i_ID=QStandardItem(str(ID))
        i_na=QStandardItem(str(na))
        i_nt=QStandardItem(str(nt))
        i_p=QStandardItem(str(p))
        i_ls=QStandardItem(str(ls))
        i_lcr=QStandardItem(str(lcr))

        #puts items into a new row
        self.output_table.output_model.appendRow([\
            i_d, \
            i_ssy, \
            i_alpha, \
            i_beta, \
            i_C, \
            i_D, \
            i_kb, \
            i_tau, \
            i_OD, \
            i_ID, \
            i_na, \
            i_nt, \
            i_p, \
            i_ls, \
            i_lcr])

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
