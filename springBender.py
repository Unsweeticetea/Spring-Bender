import sys
import sympy.physics.units as u
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

    def send_row(self, d, f_max, y_max, l0_max, R, E, A, m):
        i_d=QStandardItem()
        i_d.setData(d, role=Qt.DisplayRole)
        i_f=QStandardItem()
        i_f.setData("A", role=Qt.DisplayRole)
        self.table.model.appendRow([i_d,i_f])
        return

    def generate(self):
        '''gathers values and units then sends for calculations'''
        #Collect inputs
        #Design
        f_max = self.sb_F.value()*getattr(u, self.unitF.checkedButton().objectName()[5:])*u.acceleration_due_to_gravity
        y_max = self.sb_F.value()*getattr(u, self.unitY.checkedButton().objectName()[5:])
        l0_max = self.sb_F.value()*getattr(u, self.unitL.checkedButton().objectName()[5:])
        #Wire Diameter
        dmin=self.sb_dmin.value()
        dstep=self.sb_dstep.value()
        dmax=self.sb_dmax.value()
        #Material
        if self.unitR.checkedButton().objectName()[5:]=="Metric":
            R=self.sb_R.value()*u.pascal*1000000000
        else:
            R=self.sb_R.value()*u.psi*1000000
        if self.unitE.checkedButton().objectName()[5:]=="Metric":
            E=self.sb_E.value()*u.pascal*1000000000
        else:
            E=self.sb_E.value()*u.psi*1000000
        m=self.sb_m.value()
        A=self.sb_A.value()
        if self.unitA.checkedButton().objectName()[5:]=="Metric":
            A*=u.psi*1000*(u.inch**m)
        else:
            A*=u.pascal*1000000*(u.mm**m)
        #Properties
        ns=self.sb_ns.value()
        Xi=self.sb_Xi.value()
        self.table = OutputTable()
        self.table.show()
        for d in (d / (10**6) for d in range(int(dmin*(10**6)),int(dmax*(10**6))+1,int(dstep*(10**6)))):
            self.send_row(d, f_max, y_max, l0_max, R, E, A, m)

class OutputTable(QWidget, Ui_OutputTable):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setupUi(self)
        self.model = QStandardItemModel(0,15)
        self.model.setHorizontalHeaderLabels(
            ["d","Ssy","α","β","C","D","K_B","τ_s","OD","ID","Na","Nt","p","Ls","L_cr",]
        )
        self.tableView.setModel(self.model)
        self.tableView.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)

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