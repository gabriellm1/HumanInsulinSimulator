import math
import matplotlib.pyplot as plt
import numpy as np



# bolustime -> tempo que foi feito o bolus
# insulindu -> duração insulina ativa(parâmetro do tipo de insulina)
# Walsh curves
def IOB(bolustime,insulindur):
    if(bolustime<=0.0):
        tot=100.0
    elif (bolustime>=insulindur*60.0):
        tot=0.0
    else:
        if(insulindur==3):
            tot=-3.203e-7*math.pow(bolustime,4)+1.354e-4*math.pow(bolustime,3)-1.759e-2*math.pow(bolustime,2)+9.255e-2*bolustime+99.951
        elif (insulindur==4):
            tot=-3.31e-8*math.pow(bolustime,4)+2.53e-5*math.pow(bolustime,3)-5.51e-3*math.pow(bolustime,2)-9.086e-2*bolustime+99.95
        elif (insulindur==5):
            tot=-2.95e-8*math.pow(bolustime,4)+2.32e-5*math.pow(bolustime,3)-5.55e-3*math.pow(bolustime,2)+4.49e-2*bolustime+99.3
        elif (insulindur==6):
            tot=-1.493e-8*math.pow(bolustime,4)+1.413e-5*math.pow(bolustime,3)-4.095e-3*math.pow(bolustime,2)+6.365e-2*bolustime+99.7
    return tot

#
#
#
#
#
def intIOB(x1,x2,insulindur,bolustime):

    nn=50 #nn needs to be even
    ii=1

    #initialize with first and last terms of simpson series
    dx=(x2-x1)/nn
    intebolustimeral=IOB((bolustime-x1),insulindur)+IOB(bolustime-(x1+nn*dx),insulindur)
    integral=0
    while(ii<nn-2):
      integral = integral + 4*IOB(bolustime-(x1+ii*dx),insulindur)+2*IOB(bolustime-(x1+(ii+1)*dx),insulindur)
      ii=ii+2

    integral=integral*dx/3.0
    return(integral)


#scheiner gi curves fibolustime 7-8 from Think Like a Pancreas, fit with a trianbolustimele shaped absorbtion rate curve
#see basic math pdf on repo for details
#bolustime is time in minutes,ct is carb type
#
#  |   -----------------
#  |  /
#  | /
#  |/
#  |---------------------
#
def cob(bolustime,carbtype):

    if(bolustime<=0):
        tot=0.0
    elif (bolustime>=carbtype):
        tot=1.0
    elif ((bolustime>0) and (bolustime<=carbtype/2.0)):
        tot=2.0/math.pow(carbtype,2)*math.pow(bolustime,2)
    else:
        tot=-1.0+4.0/carbtype*(bolustime-math.pow(bolustime,2)/(2.0*carbtype))
    return(tot)

def deltatempBGI(bolustime,dbdt,sensf,insulindur,t1,t2):
    return -dbdt*sensf*((t2-t1)-1/100*intIOB(t1,t2,insulindur,bolustime))

def deltaBGC(bolustime,sensf,carbratio,carbamount,carbtype):
    return sensf/carbratio*carbamount*cob(bolustime,carbtype)

def deltaBGI(bolustime,bolus,sensf,insulindur):
    return -bolus*sensf*(1-IOB(bolustime,insulindur)/100.0)

def deltaBG(bolustime,sensf,carbratio,carbamount,carbtype,bolus,insulindur):
    return deltaBGI(bolustime,bolus,sensf,insulindur)+deltaBGC(g,sensf,carbratio,carbamount,carbtype)


if __name__ == '__main__':

    simlengthseg = 60*60
    dt = 0.1
    n = simlengthseg / dt

    simbg = [initialBG]*n
    simbgc = [0.0]*n
    simbgi = [0.0]*n

    eventCarb = [0.0]*n
    eventCarb[20/dt] = 20 #gramas
