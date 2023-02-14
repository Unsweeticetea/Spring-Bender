'''Runs a DIY desktop continuous spring bender'''
import sys
from sympy.physics.units import *
from sympy import pi, sqrt, N, Eq, symbols, solve, nsolve, solve_poly_system
from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QApplication, QMainWindow
from PyQt5.QtGui import QPalette, QColor
from ui.ui_MainWindow import Ui_MainWindow

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
        missingItem = None

        F, Y, L, E, G, A, m, ns, xi, d, OD = symbols('F Y L E G A m ns xi d OD')
        ssy, alpha, beta, C, D, kb, tau, ID, na, nt, p, Lcs, Lcr = symbols('ssy alpha beta C D kb tau ID na nt p Lcs Lcr')
        
        #test and assign design properties
        if self.sb_F.value() != 0:
            F= self.sb_F.value()
            F*= lbf if self.cb_F.currentIndex()==0 else gf
        else:
            designZeros+=1
            missingItem=F
        if self.sb_Y.value() != 0:
            Y= self.sb_Y.value()
            Y*= inch if self.cb_Y.currentIndex()==0 else mm
        else:
            designZeros+=1
            missingItem=Y
        if self.sb_L.value() != 0:
            L= self.sb_L.value()
            L*= inch if self.cb_L.currentIndex()==0 else mm
        else:
            designZeros+=1
            missingItem=L
        if self.sb_d.value() != 0:
            d= self.sb_d.value()
            d*= inch if self.cb_d.currentIndex()==0 else mm
        else:
            designZeros+=1
            missingItem=d
        if self.sb_OD.value() != 0:
            OD= self.sb_OD.value()
            OD*= inch if self.cb_OD.currentIndex()==0 else mm
        else:
            designZeros+=1
            missingItem=OD
        
        #test and assign material properties
        if self.sb_E.value() != 0:
            E= self.sb_E.value()
            E*= mega*psi if self.cb_E.currentIndex()==0 else giga*Pa
        else:
            materialZeros+=1
        if self.sb_G.value() != 0:
            G= self.sb_G.value()
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

        #eq=Eq(,convert_to(,sibase))
        #=solve(eq,eq.lhs)[0]
        eqSsy=Eq(ssy,convert_to(.45*A/(d**m),sibase))
        ssy=solve(eqSsy, eqSsy.lhs)[0]
        print('ssy:',ssy)
        
        eqAlpha=Eq(alpha,convert_to(ssy/ns,sibase))
        alpha=solve(eqAlpha, eqAlpha.lhs)[0]
        print('alpha:',alpha)
        
        eqBeta=Eq(beta,convert_to((8*(1+xi)*F)/(N(pi)*(d**2)),sibase))
        beta=solve(eqBeta, eqBeta.lhs)[0]
        print('beta:', beta)
        
        eqC=Eq(C,convert_to(((2*alpha-beta)/(4*beta))+sqrt((((2*alpha-beta)/(4*beta))**2)-(3*alpha)/(4*beta)),sibase))
        C=solve(eqC,eqC.lhs)[0]
        print('C:', C)
        
        eqD=Eq(D,convert_to(C*d,sibase))
        D=solve(eqD,eqD.lhs)[0]
        print('D:', D)
        
        eqOD=Eq(OD,convert_to(D+d,sibase))
        OD=solve(eqOD,eqOD.lhs)[0]
        print('OD:', OD)

        eqkb=Eq(kb,convert_to((4*C+2)/(4*C-3),sibase))
        kb=solve(eqkb,eqkb.lhs)[0]
        print('kb:', kb)
        
        eqtau=Eq(tau,convert_to(kb*8*F*(1+xi)*D/(N(pi)*d**3),sibase))
        tau=solve(eqtau,eqtau.lhs)[0]
        print('tau:', tau)
        
        eqID=Eq(ID,convert_to(D-d,sibase))
        ID=solve(eqID,eqID.lhs)[0]
        print('ID:', ID)
        
        eqna=Eq(na,convert_to((G*(d**4)*Y)/(8*(D**3)*F),sibase))
        na=solve(eqna,eqna.lhs)[0]
        print('na:', na)
        
        eqnt=Eq(nt,convert_to(na+2,sibase))
        nt=solve(eqnt,eqnt.lhs)[0]
        print('nt:', nt)
        
        eqp=Eq(p,convert_to((L-3*d)/na,sibase))
        p=solve(eqp,eqp.lhs)[0]
        print('p:', p)
        
        eqLcs=Eq(Lcs,convert_to(nt*d,sibase))
        Lcs=solve(eqLcs,eqLcs.lhs)[0]
        print('Ls:', Lcs)
        
        eqLcr=Eq(Lcr,convert_to(2.63*D/.5,sibase))
        Lcr=solve(eqLcr,eqLcr.lhs)[0]
        print('Lcr:', Lcr)
        
        print("solving for:",missingItem)
        equations = [eqSsy, eqAlpha, eqBeta, eqC, eqD, eqkb, eqtau, eqOD, eqID, eqna, eqnt, eqp, eqLcs, eqLcr]
        #equations = [eqOD, eqID, eqSsy, eqBeta, eqAlpha, eqC, eqD, eqkb, eqtau, eqna, eqnt, eqp, eqLcs, eqLcr]
        #outputs = [F, Y, L, d, eqSsy.lhs, eqAlpha.lhs, eqBeta.lhs, eqC.lhs, eqD.lhs, eqkb.lhs, eqtau.lhs, eqOD.lhs, eqID.lhs, eqna.lhs, eqnt.lhs, eqp.lhs, eqLcs.lhs, eqLcr.lhs]
        #outputs = [eqSsy.lhs, eqAlpha.lhs, eqBeta.lhs, eqC.lhs, eqD.lhs, eqkb.lhs, eqtau.lhs, eqOD.lhs, eqID.lhs, eqna.lhs, eqnt.lhs, eqp.lhs, eqLcs.lhs, eqLcr.lhs]
        #answer = nsolve(equations, 'OD', 1)
        answer = solve(equations, missingItem)
        #answer = solve_poly_system(equations, outputs)
        #answer = solve(equations, missingItem)[0]
        #print(solve(eqBeta, eqBeta.lhs)[0])
        print(answer)
        #print(convert_to(answer,sibase))

        '''print('ssy:',ssy, "\t\t\t\t\texpected:", convert_to(130454.9381*psi,sibase),\
            '\nalpha:',alpha, "\t\t\t\t\texpected:", convert_to(108712.4484*psi,sibase),\
            '\nbeta:', beta, "\t\t\t\t\texpected:", convert_to(9151.409228*psi,sibase),\
            '\nC:', C, "\t\t\t\t\texpected:", convert_to(10.5334876,sibase),\
            '\nD:', D, "\t\t\t\t\texpected:", convert_to(0.842679008*inch,sibase),\
            '\nkb:', kb, "\t\t\t\t\texpected:", convert_to(1.127766299,sibase),\
            '\ntau:', tau, "\t\t\t\t\texpected:", convert_to(108712.4484*psi,sibase),\
            '\nOD:', OD, "\t\t\t\t\texpected:", convert_to(0.922679008*inch,sibase),\
            '\nID:', ID, "\t\t\t\t\texpected:", convert_to(0.762679008*inch,sibase),\
            '\nna:', na, "\t\t\t\t\texpected:", convert_to(10.05359304,sibase),\
            '\nnt:', nt, "\t\t\t\t\texpected:", convert_to(12.05359304,sibase),\
            '\np:', p, "\t\t\t\t\texpected:", convert_to(0.299395449*inch,sibase),\
            '\nLs:', Lcs, "\t\t\t\t\texpected:", convert_to(0.964287444*inch,sibase),\
            '\nLcr:', Lcr, "\t\t\t\t\texpected:", convert_to(4.432491582*inch,sibase))'''

    '''def equations2(self):
        sibase = UnitSystem.get_unit_system("SI")._base_units
        designZeros = 0
        materialZeros = 0
        missingItem = None

        F, Y, L, E, G, A, m, ns, xi, d, OD = symbols('F Y L E G A m ns xi d OD')
        ssy, alpha, beta, C, D, kb, tau, ID, na, nt, p, Lcs, Lcr = symbols('ssy alpha beta C D kb tau ID na nt p Lcs Lcr')
        
        #test and assign design properties
        if self.sb_F.value() != 0:
            F= self.sb_F.value()
            F*= lbf if self.cb_F.currentIndex()==0 else gf
        else:
            designZeros+=1
            missingItem=F
        if self.sb_Y.value() != 0:
            Y= self.sb_Y.value()
            Y*= inch if self.cb_Y.currentIndex()==0 else mm
        else:
            designZeros+=1
            missingItem=Y
        if self.sb_L.value() != 0:
            L= self.sb_L.value()
            L*= inch if self.cb_L.currentIndex()==0 else mm
        else:
            designZeros+=1
            missingItem=L
        if self.sb_d.value() != 0:
            d= self.sb_d.value()
            d*= inch if self.cb_d.currentIndex()==0 else mm
        else:
            designZeros+=1
            missingItem=d
        if self.sb_OD.value() != 0:
            OD= self.sb_OD.value()
            OD*= inch if self.cb_OD.currentIndex()==0 else mm
        else:
            designZeros+=1
            missingItem=OD
        
        #test and assign material properties
        if self.sb_E.value() != 0:
            E= self.sb_E.value()
            E*= mega*psi if self.cb_E.currentIndex()==0 else giga*Pa
        else:
            materialZeros+=1
        if self.sb_G.value() != 0:
            G= self.sb_G.value()
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

        #eq=Eq(,convert_to(,sibase))
        eqSsy=Eq(ssy,convert_to(.45*A/(d**m),sibase))
        eqAlpha=Eq(alpha,convert_to(eqSsy.rhs/ns,sibase))
        eqBeta=Eq(beta,convert_to((8*(1+xi)*F)/(N(pi)*(d**2)),sibase))
        eqC=Eq(C,convert_to(((2*eqAlpha.rhs-eqBeta.rhs)/(4*eqBeta.rhs))+sqrt((((2*eqAlpha.rhs-eqBeta.rhs)/(4*eqBeta.rhs))**2)-(3*eqAlpha.rhs)/(4*eqBeta.rhs)),sibase))
        eqD=Eq(D,convert_to(eqC.rhs*d,sibase))
        eqkb=Eq(kb,convert_to((4*eqC.rhs+2)/(4*eqC.rhs-3),sibase))
        eqtau=Eq(tau,convert_to(eqkb.rhs*8*F*(1+xi)*eqD.rhs/(N(pi)*d**3),sibase))
        eqOD=Eq(OD,convert_to(eqD.rhs+d,sibase))
        eqID=Eq(ID,convert_to(eqD.rhs-d,sibase))
        eqna=Eq(na,convert_to((G*(d**4)*Y)/(8*(eqD.rhs**3)*F),sibase))
        eqnt=Eq(nt,convert_to(eqna.rhs+2,sibase))
        eqp=Eq(p,convert_to((L-3*d)/eqna.rhs,sibase))
        eqLcs=Eq(Lcs,convert_to(eqnt.rhs*d,sibase))
        eqLcr=Eq(Lcr,convert_to(2.63*eqD.rhs/.5,sibase))

        print("solving for:",missingItem)
        equations = [eqSsy, eqAlpha, eqBeta, eqC, eqD, eqkb, eqtau, eqOD, eqID, eqna, eqnt, eqp, eqLcs, eqLcr]
        #equations = [eqOD, eqID, eqSsy, eqBeta, eqAlpha, eqC, eqD, eqkb, eqtau, eqna, eqnt, eqp, eqLcs, eqLcr]
        #outputs = [F, Y, L, d, eqSsy.lhs, eqAlpha.lhs, eqBeta.lhs, eqC.lhs, eqD.lhs, eqkb.lhs, eqtau.lhs, eqOD.lhs, eqID.lhs, eqna.lhs, eqnt.lhs, eqp.lhs, eqLcs.lhs, eqLcr.lhs]
        outputs = [eqSsy.lhs, eqAlpha.lhs, eqBeta.lhs, eqC.lhs, eqD.lhs, eqkb.lhs, eqtau.lhs, eqOD.lhs, eqID.lhs, eqna.lhs, eqnt.lhs, eqp.lhs, eqLcs.lhs, eqLcr.lhs]
        #answer = nsolve(equations, 'OD', 1)
        #answer = solve(equations, missingItem, minimal = True)
        answer = solve_poly_system(equations, outputs)
        #answer = solve(equations, missingItem)[0]
        #print(solve(eqBeta, eqBeta.lhs)[0])
        print(answer)
        #print(convert_to(answer,sibase))
         
        #=solve(eq,eq.lhs)[0]
        ssy=solve(eqSsy, eqSsy.lhs)[0]
        alpha=solve(eqAlpha, eqAlpha.lhs)[0]
        beta=solve(eqBeta, eqBeta.lhs)[0]
        C=solve(eqC,eqC.lhs)[0]
        D=solve(eqD,eqD.lhs)[0]
        kb=solve(eqkb,eqkb.lhs)[0]
        tau=solve(eqtau,eqtau.lhs)[0]
        OD=solve(eqOD,eqOD.lhs)[0]
        ID=solve(eqID,eqID.lhs)[0]
        na=solve(eqna,eqna.lhs)[0]
        nt=solve(eqnt,eqnt.lhs)[0]
        p=solve(eqp,eqp.lhs)[0]
        Lcs=solve(eqLcs,eqLcs.lhs)[0]
        Lcr=solve(eqLcr,eqLcr.lhs)[0]
        print('ssy:',ssy, "\t\t\t\t\texpected:", convert_to(130454.9381*psi,sibase),\
            '\nalpha:',alpha, "\t\t\t\t\texpected:", convert_to(108712.4484*psi,sibase),\
            '\nbeta:', beta, "\t\t\t\t\texpected:", convert_to(9151.409228*psi,sibase),\
            '\nC:', C, "\t\t\t\t\texpected:", convert_to(10.5334876,sibase),\
            '\nD:', D, "\t\t\t\t\texpected:", convert_to(0.842679008*inch,sibase),\
            '\nkb:', kb, "\t\t\t\t\texpected:", convert_to(1.127766299,sibase),\
            '\ntau:', tau, "\t\t\t\t\texpected:", convert_to(108712.4484*psi,sibase),\
            '\nOD:', OD, "\t\t\t\t\texpected:", convert_to(0.922679008*inch,sibase),\
            '\nID:', ID, "\t\t\t\t\texpected:", convert_to(0.762679008*inch,sibase),\
            '\nna:', na, "\t\t\t\t\texpected:", convert_to(10.05359304,sibase),\
            '\nnt:', nt, "\t\t\t\t\texpected:", convert_to(12.05359304,sibase),\
            '\np:', p, "\t\t\t\t\texpected:", convert_to(0.299395449*inch,sibase),\
            '\nLs:', Lcs, "\t\t\t\t\texpected:", convert_to(0.964287444*inch,sibase),\
            '\nLcr:', Lcr, "\t\t\t\t\texpected:", convert_to(4.432491582*inch,sibase))'''


    def test(self):
        
        self.sb_F.setValue(20)
        self.cb_F.setCurrentIndex(0)
        self.sb_Y.setValue(2)
        self.cb_Y.setCurrentIndex(0)
        self.sb_L.setValue(3.25)
        self.cb_L.setCurrentIndex(0)
        self.sb_d.setValue(.08)
        self.cb_d.setCurrentIndex(0)
        self.sb_OD.setValue(0.922679008)
        self.cb_OD.setCurrentIndex(0)
        #determine OD
        if False:
            self.sb_OD.setValue(0)
            self.cb_OD.setCurrentIndex(0)
        elif False:
            self.sb_d.setValue(0)
            self.cb_d.setCurrentIndex(0)
        elif True:
            self.sb_Y.setValue(0)
            self.cb_Y.setCurrentIndex(0)
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
