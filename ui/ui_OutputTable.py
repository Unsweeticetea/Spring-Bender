# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'C:\Users\iamze\Github\Spring-Bender\ui\ui_OutputTable.ui'
#
# Created by: PyQt5 UI code generator 5.15.7
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_OutputTable(object):
    def setupUi(self, OutputTable):
        OutputTable.setObjectName("OutputTable")
        OutputTable.resize(789, 406)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.MinimumExpanding, QtWidgets.QSizePolicy.MinimumExpanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(OutputTable.sizePolicy().hasHeightForWidth())
        OutputTable.setSizePolicy(sizePolicy)
        self.gridLayout = QtWidgets.QGridLayout(OutputTable)
        self.gridLayout.setObjectName("gridLayout")
        self.tableView = QtWidgets.QTableView(OutputTable)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.tableView.sizePolicy().hasHeightForWidth())
        self.tableView.setSizePolicy(sizePolicy)
        self.tableView.setSortingEnabled(True)
        self.tableView.setObjectName("tableView")
        self.gridLayout.addWidget(self.tableView, 0, 0, 1, 1)

        self.retranslateUi(OutputTable)
        QtCore.QMetaObject.connectSlotsByName(OutputTable)

    def retranslateUi(self, OutputTable):
        _translate = QtCore.QCoreApplication.translate
        OutputTable.setWindowTitle(_translate("OutputTable", "Output"))
