from Tkinter import *
import tkFileDialog
import matplotlib.pyplot as plot
import numpy as np
import pickle

def loadDataFromFile(inputFilename):
    '''
        Loads a saved file with captured curves. Returns
        data as two vectors of: ([[curves]], [voltages])
    '''

    with open(inputFilename, 'rb') as input:
        _loadedData1 = pickle.load(input)
        _loadedData2 = pickle.load(input)

    return _loadedData1, _loadedData2


def askForInputFile(fileFilter="*.*"):
    root = Tk()
    root.withdraw()
    _filename = tkFileDialog.askopenfilename(filetypes = [("Open file ...", fileFilter)])
    return _filename


if __name__ == "__main__":

    print "-- Tube capture utility --\nFile reading test program\nIgnacio D. Simon 2016\n"

    # Ask for input file
    _filename = askForInputFile(fileFilter="*.cur")
    if _filename == "":
        print "No input file selected. Quitting."
        quit()

    # Load data from file
    [_readCurves, _readVoltages] = loadDataFromFile(_filename)

    for i in range(0, len(_readVoltages)):
        _voltages = []
        _currents = []
        for j in range(0, len(_readCurves[i])):
           _voltages.extend([_readCurves[i][j][0]])
           _currents.extend([_readCurves[i][j][1]])

        print _voltages, _currents
        print '\n'
        plot.plot(_voltages, _currents, 'x-')

    plot.grid()
    plot.xlabel('Vak [V]')
    plot.ylabel('Iak [mA]')
    plot.title('Captured curves from loaded file:\n"%s"' % _filename)
    plot.show()



















