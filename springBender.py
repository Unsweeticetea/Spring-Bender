import sys
import sympy.physics.units as u
from sympy import pi, Symbol, symbols, sqrt, N
from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QApplication, QMainWindow, QWidget, QTableView, QHeaderView
from PyQt5.QtGui import QPalette, QColor, QStandardItemModel, QStandardItem
from ui.ui_MainWindow import Ui_MainWindow
from ui.ui_OutputTable import  Ui_OutputTable

class MainWindow(QMainWindow, Ui_MainWindow):
    '''main window'''
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setupUi(self)
        self.tabs.removeTab(3)
        self.tabs.removeTab(2)
        self.tabs.setCurrentIndex(0)
        self.connections()
        self.output_table = OutputTable()

    def connections(self):
        '''Programs all connections for main window'''
        self.btn_Generate.clicked.connect(self.generate)
        self.chk_advanced.clicked.connect(self.show_tabs)

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
        print(E, G, D)
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

    def generate(self):
        '''gathers values and units then sends for calculations'''
        #Collect inputs
        #Design
        f_max = self.sb_F.value()*getattr(u, self.unitF.checkedButton().objectName()[5:])*u.convert_to((9.80665*u.m/u.s**2),u.N)
        y_max = self.sb_F.value()*getattr(u, self.unitY.checkedButton().objectName()[5:])
        l0 = self.sb_F.value()*getattr(u, self.unitL.checkedButton().objectName()[5:])
        #Wire Diameter
        dmin=int(self.sb_dmin.value()*(10**6))
        dstep=int(self.sb_dstep.value()*(10**6))
        dmax=int(self.sb_dmax.value()*(10**6))
        #Material
        if self.unitE.checkedButton().objectName()[5:]=="Metric":
            print("E Metric")
            E=self.sb_E.value()*u.pascal*1000000000
        else:
            print("E Imperial")
            E=self.sb_E.value()*u.psi*1000000
        
        if self.unitG.checkedButton().objectName()[5:]=="Metric":
            print("G Metric")
            G=self.sb_G.value()*u.pascal*1000000000
        else:
            print("G Imperial")
            G=self.sb_G.value()*u.psi*1000000
        
        m=self.sb_m.value()
        A=self.sb_A.value()
        if self.unitA.checkedButton().objectName()[5:]=="Metric":
            A*=u.psi*1000*(u.inch**m)
        else:
            A*=u.pascal*1000000*(u.mm**m)
        #Properties
        ns=self.sb_ns.value()
        xi=self.sb_Xi.value()
        
        #generates expressions and creates output table
        self.output_table = OutputTable()
        self.output_table.show()

        #populates output table
        for d in (d / (10**6) for d in range(dmin,dmax+1,dstep)):
            self.send_row(d, f_max, y_max, l0, G, E, A, m, ns, xi)

    def test(self):
        '''tests output generation with known values'''
        #Collect inputs
        self.chk_ground.setChecked(True)
        #Design
        f_max = 20*u.pound*u.convert_to((9.80665*u.m/u.s**2),u.N)
        y_max = 2*u.inch
        l0 = 3.25*u.inch
        #Wire Diameter
        dmin=int(.055*(10**6))
        dstep=int(.005*(10**6))
        dmax=int(.105*(10**6))
        #Material
        E=28.5*u.psi*1000000
        G=11.75*u.psi*1000000
        m=.145
        A=201*u.psi*1000*(u.inch**m)
        #Properties
        ns=1.2
        xi=.15
        
        #generates expressions and creates output table
        self.output_table.show()

        #populates output table
        for d in (d / (10**6) for d in range(dmin,dmax+1,dstep)):
            self.send_row(d, f_max, y_max, l0, G, E, A, m, ns, xi)

class OutputTable(QWidget, Ui_OutputTable):
    '''Table for storing calculated values'''
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setupUi(self)
        self.output_model = QStandardItemModel(0,15)
        self.output_model.setHorizontalHeaderLabels(
            ["d","Ssy","α","β","C","D","K_B","τ_s","OD","ID","Na","Nt","p","Ls","L_cr"]
        )
        self.tableView.setModel(self.output_model)
        self.tableView.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)

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
