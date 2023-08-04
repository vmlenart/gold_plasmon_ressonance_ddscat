import os
from numpy import *
import matplotlib.pyplot as plt
from itertools import islice
import fileinput
import sys

def replaceAll(file,searchExp,replaceExp):
    for line in fileinput.input(file, inplace=1):
        if searchExp in line:
            line = line.replace(searchExp,replaceExp)
        sys.stdout.write(line)

def seek_to_line(ark, n):
    for ignored_line in islice(ark, n-1):
        pass

d = 2 # dipole spacing
n = 101 # number of points of interpolation
serie = 1

l = 0.191

for i in range(900):
    replaceAll("ddscat.par",str(l-0.001),str(l))

    os.system('./ddscat7.1/src/ddscat')
    print 'ddscat concluded for %f' %l

    os.system('./ddscat7.1/src/ddfield')
    print 'ddfield concluded for %f' %l

    x = []
    y = []
    Exr = []
    Eyr = []
    Ezr = []
    Exi = []
    Eyi = []
    Ezi = []

    myfile = open('ddfield.E', 'rb')
    seek_to_line(myfile,24)

    for line in myfile:
        row = line.split()
        x.append(float(row[0]))
        y.append(float(row[1]))
        Exr.append(float(row[3]))
        Eyr.append(float(row[5]))
        Ezr.append(float(row[7]))
        Exi.append(float(row[4]))
        Eyi.append(float(row[6]))
        Ezi.append(float(row[8]))

    x = array(x)
    y = array(y)
    Exr = array(Exr)
    Exi = array(Exi)
    Eyr = array(Eyr)
    Eyi = array(Eyi)
    Ezr = array(Ezr)
    Ezi = array(Ezi)

    Ex = sqrt(Exr**2 + Exi**2)
    Ey = sqrt(Eyr**2 + Eyi**2)
    Ez = sqrt(Ezr**2 + Ezi**2)

    E = sqrt(Ex**2 + Ey**2 + Ez**2)

    X = d*x.reshape(n,n)
    Y = d*y.reshape(n,n)
    Z = E.reshape(n,n)

    plt.pcolor(X, Y, Z)
    cb = plt.colorbar()
    plt.axis([-100,100,-100,100])
    plt.xlabel('$x$ nm')
    plt.ylabel('$y$ nm')
    cb.set_label('$E/E_0$')
    plt.title('$\lambda$ = %f' %l + '$\mu$m')
    filename = str('E\%03d'%serie) + '.png'
    plt.savefig(filename)
    plt.clf()


    x = []
    y = []
    Bxr = []
    Byr = []
    Bzr = []
    Bxi = []
    Byi = []
    Bzi = []

    myfile = open('ddfield.B', 'rb')
    seek_to_line(myfile,24)

    for line in myfile:
        row = line.split()
        x.append(float(row[0]))
        y.append(float(row[1]))
        Bxr.append(float(row[3]))
        Byr.append(float(row[5]))
        Bzr.append(float(row[7]))
        Bxi.append(float(row[4]))
        Byi.append(float(row[6]))
        Bzi.append(float(row[8]))

    x = array(x)
    y = array(y)
    Bxr = array(Bxr)
    Bxi = array(Bxi)
    Byr = array(Byr)
    Byi = array(Byi)
    Bzr = array(Bzr)
    Bzi = array(Bzi)

    Bx = sqrt(Bxr**2 + Bxi**2)
    By = sqrt(Byr**2 + Byi**2)
    Bz = sqrt(Bzr**2 + Bzi**2)

    B = sqrt(Bx**2 + By**2 + Bz**2)

    X = d*x.reshape(n,n)
    Y = d*y.reshape(n,n)
    Z = B.reshape(n,n)

    plt.pcolor(X, Y, Z)
    cb = plt.colorbar()
    plt.axis([-100,100,-100,100])
    plt.xlabel('$x$ nm')
    plt.ylabel('$y$ nm')
    cb.set_label('$B/B_0$')
    plt.title('$\lambda$ = %f' %l + '$\mu$m')
    filename = str('B\%03d'%serie) + '.png'
    plt.savefig(filename)
    plt.clf()

    serie+=1
    l+=0.001

