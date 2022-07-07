# from PIL import Image
# import PySpin
# import time
# import myspincam
# import pulse_generator as pulser
# import os
# import experiment_cmds
import sys
import PWM_Acquisition
from PyQt5 import QtGui
from PyQt5.QtWidgets import QApplication, QWidget, QLineEdit, QPushButton, QTextEdit, QVBoxLayout, QMainWindow
from PyQt5.QtGui import QIcon

if __name__ == "__main__":
    sys.setrecursionlimit(10000)
    app = QApplication(sys.argv)
    MainWindow = QMainWindow()
    ui = PWM_Acquisition.Ui_Dialog()
    ui.setupUi(MainWindow)
    MainWindow.show()
    sys.exit(app.exec())

import moveFiles

