import numpy as np
import matplotlib.pyplot as plt

def rcosine(beta, Tbaud, oversampling, Nbauds, Norm):
    """ Respuesta al impulso del pulso de caida cosenoidal """
    t_vect = np.arange(-0.5*Nbauds*Tbaud, 0.5*Nbauds*Tbaud, float(Tbaud)/oversampling)
    
    y_vect = []
    for t in t_vect:
        y_vect.append(np.sinc(t/Tbaud)*(np.cos(np.pi*beta*t/Tbaud)/(1-(4.0*beta*beta*t*t/(Tbaud*Tbaud)))))

    y_vect =np.array(y_vect)

    
    if(Norm):
        return (t_vect, y_vect/y_vect.sum())
    else:
        return (t_vect,y_vect)

def eyediagram(data, n, offset, period):
    span     = (2*n)
    segments = int(len(data)/span)
    xmax     = (n-1)*period
    xmin     = -(n-1)*period
    x        = list(np.arange(-n,n,)*period)
    xoff     = offset


    plt.figure()
    plt.title("Diagrama de Ojos de la senal a transmitir")
    for i in range(0,segments-1):
        plt.plot(x,data[(i*span+xoff):((i+1)*span+xoff)],'b')
        plt.hold(True)
        plt.grid(True)
      
    plt.xlim(xmin, xmax)

def rrcosine(beta, T, oversampling, Nbauds):

    t_vect = np.arange(-0.5*Nbauds*T, 0.5*Nbauds*T, float(T)/oversampling)  
    y_vect = []
    A = []
    B = []
    C = []

    for t in t_vect:
        A.append(np.sin(np.pi*(t/T)*(1-beta)))
        B.append((4*beta*(t/T))*np.cos(np.pi*(t/T)*(1+beta)))
        C.append((np.pi*(t/T))*(1-(4*beta*(t/T))**2))

    A = np.array(A)
    B = np.array(B)
    C = np.array(C)
    
    g=(A+B)/C
    g0=(1-beta)+((4*beta)/np.pi)
    g1=(beta/np.sqrt(2))*((1+2/np.pi)*np.sin(np.pi/(4*beta))+(1-2/np.pi)*np.cos(np.pi/(4*beta)))

    p=np.zeros(len(t_vect))
    p=g

    for ptr in range(len(t_vect)):
        if (t_vect[ptr]==0):
            p[ptr]=g0
        elif (t_vect[ptr]==(oversampling/(4*beta)) or t_vect[ptr]==(-oversampling/(4*beta))):
            p[ptr]=g1
                   

    return (t,p)


