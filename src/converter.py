import sys
from PyQt5.QtWidgets import QApplication
from PyQt5.QtWidgets import QMainWindow
from PyQt5.QtWidgets import QFileDialog

from src.flowtocora import LinHybridFlowStarToCORA, NonLinHybridFlowToCORA, LinContFlowToCORA, NonLinContFlowToCORA
from src.coratoflow import CORAtoFlowStar
from src.Ui_MainWindow import Ui_MainWindow

infile_set = False


class MainWindow:

    def __init__(self):
        self.main_win = QMainWindow()
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self.main_win)
        self.ui.infile_button.clicked.connect(self.__browseInput)
        self.ui.outfile_button.clicked.connect(self.__browseOutput)
        self.ui.okay_button.clicked.connect(self.__enableConverter)
        self.ui.convert_button.clicked.connect(self.__startConverter)




    def show(self):
        self.main_win.show()

    def __browseInput(self):
        global infile_set
        filePath = QFileDialog.getOpenFileName(parent=self.main_win, caption="Select Input File",directory='', filter="Models (*.m *.model)")

        if filePath[0] != "":
            infile_set = True
            self.ui.input_line.setText(filePath[0])
            self.ui.input_line.setStyleSheet("color: black;")
            path = ""

            if '.model' in str(filePath[0]):
                path += str(filePath[0]).replace('.model', '.m')
            else:
                path += str(filePath[0]).replace('.m', '.model')

            self.ui.output_line.setText(path)
            self.ui.output_line.setStyleSheet("color: black;")

    def __browseOutput(self):
        filePath = QFileDialog.getOpenFileName(parent=self.main_win, caption="Select Output File",directory='', filter="*.*")
        self.ui.output_line.setText(filePath[0])
        self.ui.output_line.setStyleSheet("color: black;")

    def __enableConverter(self):
        global infile_set
        self.ui.tabs.setCurrentWidget(self.ui.files_widget)
        if infile_set:
            self.ui.convert_button.setEnabled(True)

    def __startConverter(self):
        infile_path = self.ui.input_line.text()
        outfile_path = self.ui.output_line.text()

        options = self.__getOptions()

        if '.model' in infile_path:
            print("Converting Flow* file to CORA file...")
            converter = FlowStarToCORA()
        else:
            print("Converting CORA file to Flow* file")
            converter = CORAtoFlowStar()

        converter.convert(infile_path, outfile_path, options)

    def __getOptions(self):
        system = self.ui.systemGroup.checkedButton()
        taylor = self.ui.taylor.text()
        zonotope = self.ui.zonotope.text()
        polytope = self.ui.polytope.text()
        guard = self.ui.guard.text()
        reduction = self.ui.reduction.text()
        en = self.ui.enclosure_button.isChecked()
        hy = self.ui.enclosure_button.isChecked()
        ori = self.ui.origin_button.isChecked()

        enclosure = str(0)
        hyperplane = str(0)
        origin = str(0)

        if en:
            enclosure = str(1)
        if hy:
            hyperplane = str(1)
        if ori:
            origin = str(1)

        options = {'system': system,
                   'taylor': taylor,
                   'zonotope': zonotope,
                   'polytope': polytope,
                   'guard': guard,
                   'reduction': reduction,
                   'enclosure': enclosure,
                   'hyperplane': hyperplane,
                   'origin': origin}

        return options


if __name__ == '__main__':
    # app = QApplication(sys.argv)
    # main_win = MainWindow()
    # main_win.show()
    #sys.exit(app.exec_())
    options = {'system': 'linear hybrid',
               'taylor': '10',
               'zonotope': '20',
               'polytope': '10',
               'guard': 'polytope',
               'reduction': 'girard',
               'enclosure': '5',
               'hyperplane': '0',
               'origin': '0'}

    if options['system'] == 'linear hybrid':
        conv = LinHybridFlowStarToCORA()
    elif options['system'] == 'non-linear hybrid':
        conv = NonLinHybridFlowToCORA()
    elif options['system'] == 'linear continuous':
        conv = LinContFlowToCORA()
    else:
        conv = NonLinContFlowToCORA()

    conv.convert('/Users/marta/Desktop/Hybrid linear/rod_reactor.model', '/Users/marta/Desktop/Hybrid linear/switching_5.m',options)



