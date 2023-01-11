# -*- coding: utf-8 -*-
##############################################################
#
#  Identify peaks in SAXS values
#  @Anderson Ferreira Sepulveda
#  #Doctorate in Biosystems Program
#  #Federal University of ABC
#
#      Mar 2018
##############################################################
import matplotlib.pyplot as plt
import math
import pylab as plb
import numpy as np
from scipy.optimize import curve_fit



    
def main():
    
    
    arq = r'C:\Users\ander\Documents\Aryane_amostras_M\Manipulated\P30_IPM_LEC_CUR1_40C.txt'
    
    q = []
    X = []   
    SIH =[]
    
    ##################################################################
    #Choose the best distance between peak and parameter calculated
    #################################################################
    var = 0.025
    
    ################################################################
    #   If "math domain error" or "RuntimeError", change X.append()
    #   put or remove '#'
    ##############################################################
    Q = open(arq, 'r')
    i = Q.readline()
    for i in Q:
        d = i.split()
        q.append(float(d[0]))
        #q.append(math.log(float(d[0]),10))
        X.append(float(d[1]))
        #X.append(math.log(float(d[1]) + 0.0001,10))
        #X.append(float(d[1]))
        
    
    
    peaks = []
    min_peaks = []
    
    ###################################################################
    #  div is a parameter which divides the date in peaces
    #  If you see few identified peaks, choose bigger div        
    
    div = 21
    
    ###################################################################
    lar = len(X)//div    
    j = 0
    NX = []
    q_max =[]
    X_max =[]
    q_min =[]
    X_min = []
    for i in range(0,div):
        NX = []        
        k = j        
        while j < k + lar:            
            NX.append(X[j])                        
            j += 1        
        if len(NX) != 0:
            peaks.append(max(NX))
            min_peaks.append(min(NX))
            
    #Max intensity and Max q
    for i in range(0,len(peaks)):
        for j in range(0,len(X)):
            if peaks[i] == X[j]:
                q_max.append(q[j])
                X_max.append(j)
    
    #Min intensity and Min q
    curv_ini = 2                                                         #choose first min peak
    curv_fin = 5                                                 #choose second min peak
    for i in range(0,len(X)-1):
        if X[i] == min_peaks[curv_ini]:
            X_min.append(X[i])
            q_min.append(q[i])
            while X[i+1] != min_peaks[curv_fin]:
                X_min.append(X[i+1])
                q_min.append(q[i+1])
                i += 1
            break    
    for i in range(0,len(X)):
        if X[i] == min_peaks[curv_fin]:
            X_min.append(X[i])
            q_min.append(q[i])
    
    
     
    first = 4
    second = 5
    second = len(peaks)             
    Mpos = peaks.index(max(peaks[first:second]))
    QMax = q_max[Mpos]                              #First peak of your data'''
    
    print(Mpos)
    print(QMax)
    
    
    ###########################################################
    # the Gaussian curves
    ##########################################################
    x = np.array(q_min)
    y = np.array(X_min)
    
    mean = sum(x*y)/sum(y)
    sigma = np.sqrt(sum(y * (x - mean)**2) / sum(y))
    
    def Gaussian(x,A,xc,sigma):
        '''
        Gaussian curve
        '''
        #print xc
        return A*np.exp(-(x - xc)**2/(2*sigma**2))
    
    def peak_Gaussian(x,A,xc,sigma):
        '''
        Identify peak in Gaussian
        '''
        #print xc
        return xc, A*np.exp(-(x - xc)**2/(2*sigma**2))
    
    #popt, pcov = curve_fit(Gaussian,x,y,p0=[max(y),mean,sigma])
    popt, pcov = curve_fit(Gaussian,x,y,p0=[max(y),QMax,sigma])
    
    E, curv = peak_Gaussian(x,*popt)
    
    qMax = QMax #E                                                                  #First peak of the gaussian fit
    
    #error_q = abs(qMax - QMax)
    error_q = sigma/2
    
    print(qMax)
    print(sigma)
                    
    ##################################################################
    #   supramolecular organization
    ##################################################################
    
    
    #Lamellar 1
    L1 = qMax
    
    #Lamellar 2
    Lam2 = qMax 
    
    #Hexagonal main peak
    HMP = qMax
    
    #Pm3n Cubic
    Pm3n = qMax
            
    #Pn3m Cubic
    Pn3m = qMax
    
    #Fm3m
    Fm3m = qMax
    
    #Fd3m
    Fd3m = qMax
    
    #la3d cubic
    la3d = qMax
    
    #lm3m
    lm3m = qMax
    
    #Cubic 203, 227
    Cubic = qMax
    
    #FCC
    fcc = qMax
    
    #BCC
    bcc = qMax
    
    #HCP
    hcp = qMax
    
    #####################################################
    #  Ribbon phase
    #  Put #'s at each line to inactivate function
    #####################################################
    rph_a = 50.                                                                 #angstrons
    rph_b = 95.                                                                 #angstrons
    
    #####################################################
    #  Ripple phase
    #  Put #'s at each line to inactivate function
    #####################################################
    rip_a = 55.3                                                                #angstrons
    rip_b = 85.3                                                                #angstrons
    rip_gamma = 110                                                             #degrees
    
    ###############################################################################
    #  Calculated parameters
    ##############################################################################
    
    #Inverse hexagonal
    SIH = []
    SH = []
    h_IH = [1.,1.,2.,2.,2.,2.,3.,4.,3.,4.]
    k_IH = [0.,1.,0.,1.,0.,2.,1.,0.,2.,1.]
    for i in range(len(h_IH)):
        SH.append(math.sqrt((h_IH[i])**2 + h_IH[i]*k_IH[i] + (k_IH[i])**2))        
    SIH = [HMP, HMP, HMP]
    print(SH)
    for i in range(len(SH)):
        SIH.append(HMP*SH[i])
        SIH.append(HMP*SH[i])
        SIH.append(HMP*SH[i])
    
    d10_hmp = 2*math.pi/HMP
    d11_hmp = d10_hmp/math.sqrt(3)
    
    error_d10_hmp = (2*math.pi/(HMP)**2)*error_q
    error_d11_hmp = (1/3.)*error_q  
    
    
    #Lamellar 1
    SL1 = []
    SL = [1., 2., 3., 4., 5., 6.]
    SL1 = [L1, L1, L1]
    for i in range(len(SL)):
        SL1.append(L1*SL[i])
        SL1.append(L1*SL[i])
        SL1.append(L1*SL[i])
        
    d10_SL1 = 2*math.pi/L1
    
    error_d10_SL1 = (2*math.pi/(L1)**2)*error_q
    
    #Lamellar 2
    SL2 = []
    SL2 = [Lam2, Lam2, Lam2]
    for i in range(len(SL)):
        SL2.append(Lam2*SL[i])
        SL2.append(Lam2*SL[i])
        SL2.append(Lam2*SL[i])
    
    d10_SL2 = 2*math.pi/Lam2
    
    error_d10_SL2 = (2*math.pi/Lam2**2)*error_q
    
    #Pm3n    
    q_Pm3n = []
    h_Pm3n = [1.,1.,1.,2.,2.,2.,2.,2.,3.,3.,2.,3.,3.,4.,4.,4.,3.,4.,3.,4.,4.,3.]
    k_Pm3n = [0.,1.,1.,0.,1.,1.,2.,2.,1.,1.,2.,2.,2.,0.,1.,1.,3.,2.,3.,2.,3.,3.]
    l_Pm3n = [0.,0.,1.,0.,0.,1.,0.,1.,0.,1.,2.,0.,1.,0.,0.,1.,1.,0.,2.,2.,1.,3.]
    for i in range(len(h_Pm3n)):
        q_Pm3n.append(math.sqrt((h_Pm3n[i])**2 + (k_Pm3n[i])**2 + (l_Pm3n[i])**2))
    Pm3n_cubic = [Pm3n,Pm3n,Pm3n]
    for i in range(1,len(h_Pm3n)):
        if i != 2 and i != 7 and i != 9 and i != 16 and i != 21:
            Pm3n_cubic.append(Pm3n*(q_Pm3n[i]/q_Pm3n[1]))
            Pm3n_cubic.append(Pm3n*(q_Pm3n[i]/q_Pm3n[1]))
            Pm3n_cubic.append(Pm3n*(q_Pm3n[i]/q_Pm3n[1]))
            
    a_Pm3n = 2*math.pi*math.sqrt(2)/Pm3n
    
    error_a_Pm3n = (2*math.pi*math.sqrt(2)/(Pm3n)**2)*error_q

        
    #Pn3m
    q_Pn3m = []
    h_Pn3m = [1.,1.,1.,2.,2.,2.,2.,2.,3.,3.,2.,3.,3.,4.,4.,4.,3.,4.,3.,4.,4.,3.]
    k_Pn3m = [0.,1.,1.,0.,1.,1.,2.,2.,1.,1.,2.,2.,2.,0.,1.,1.,3.,2.,3.,2.,3.,3.]
    l_Pn3m = [0.,0.,1.,0.,0.,1.,0.,1.,0.,1.,2.,0.,1.,0.,0.,1.,1.,0.,2.,2.,1.,3.]
    for i in range(len(h_Pn3m)):
        q_Pn3m.append(math.sqrt((h_Pn3m[i])**2 + (k_Pn3m[i])**2 + (l_Pn3m[i])**2))
    P_cubic = [Pn3m,Pn3m,Pn3m]
    for i in range(1,len(h_Pn3m)):
        if i != 4 and i != 7 and i != 9 and  i != 17 and i != 21:
            P_cubic.append(Pn3m*(q_Pn3m[i]/q_Pn3m[1]))
            P_cubic.append(Pn3m*(q_Pn3m[i]/q_Pn3m[1])) 
            P_cubic.append(Pn3m*(q_Pn3m[i]/q_Pn3m[1]))
                     
    a_Pn3m = 2*math.pi*math.sqrt(2)/Pn3m
    
    error_a_Pn3m = (2*math.pi*math.sqrt(2)/(Pn3m)**2)*error_q
    
    #Fm3m
    q_Fm3m = []
    h_Fm3m = [1.,1.,1.,2.,2.,2.,2.,2.,3.,3.,2.,3.,3.,4.,4.,4.,3.,4.,3.,4.,4.,3.]
    k_Fm3m = [0.,1.,1.,0.,1.,1.,2.,2.,1.,1.,2.,2.,2.,0.,1.,1.,3.,2.,3.,2.,3.,3.]
    l_Fm3m = [0.,0.,1.,0.,0.,1.,0.,1.,0.,1.,2.,0.,1.,0.,0.,1.,1.,0.,2.,2.,1.,3.]
    for i in range(len(h_Fm3m)):
        q_Fm3m.append(math.sqrt((h_Fm3m[i])**2 + (k_Fm3m[i])**2 + (l_Fm3m[i])**2))
    Fm3m_cubic = [Fm3m, Fm3m, Fm3m]
    for i in range(len(h_Fm3m)):
           if i != 0 and i != 1 and i != 4 and i != 5 and i != 7 and i != 8 and i != 11 and i != 12 and i != 14 and i != 15 and i != 18 and i != 20:
               Fm3m_cubic.append(Fm3m*(q_Fm3m[i]/q_Fm3m[2]))
               Fm3m_cubic.append(Fm3m*(q_Fm3m[i]/q_Fm3m[2]))
               Fm3m_cubic.append(Fm3m*(q_Fm3m[i]/q_Fm3m[2]))
               
    a_Fm3m = 2*math.pi*math.sqrt(3)/Fm3m
    
    error_a_Fm3m = (2*math.pi*math.sqrt(3)/(Fm3m)**2)*error_q
    
    #Fd3m
    q_Fd3m = []
    h_Fd3m = [1,1,1,2,2,2,2,2,3,3,2,3,3,4,4,4,3,4,3,4,4,3]
    k_Fd3m = [0,1,1,0,1,1,2,2,1,1,2,2,2,0,1,1,3,2,3,2,3,3]
    l_Fd3m = [0,0,1,0,0,1,0,1,0,1,2,0,1,0,0,1,1,0,2,2,1,3]
    for i in range(len(h_Fm3m)):
        q_Fd3m.append(math.sqrt((h_Fd3m[i])**2 + (k_Fd3m[i])**2 + (l_Fd3m[i])**2))
    Fd3m_cubic = [Fd3m, Fd3m, Fd3m]
    for i in range(1,len(h_Fd3m)):
           if i != 0 and i != 1 and i != 3 and i != 4 and i != 5 and i != 7 and i != 8 and i != 11 and i != 12 and i != 14 and i != 15 and i != 17 and i != 18 and i != 20:
               Fd3m_cubic.append(Fd3m*(q_Fd3m[i]/q_Fd3m[2]))
               Fd3m_cubic.append(Fd3m*(q_Fd3m[i]/q_Fd3m[2]))
               Fd3m_cubic.append(Fd3m*(q_Fd3m[i]/q_Fd3m[2]))
               
    a_Fd3m = 2*math.pi*math.sqrt(3)/Fd3m
    
    error_a_Fd3m = (2*math.pi*math.sqrt(3)/(Fd3m)**2)*error_q
                                                                              
    #la3d
    q_la3d = []
    h_la3d = [1,1,2,2,2,2,2,3,3,2,3,3,4,4,4,3,4,3,4,4,3,5,6,5,4,5]
    k_la3d = [1,1,0,1,1,2,2,1,1,2,2,2,0,1,1,3,2,3,2,3,3,2,1,4,4,4]
    l_la3d = [0,1,0,0,1,0,1,0,1,2,0,1,0,0,1,1,0,2,2,1,3,1,1,1,4,3]
    for i in range(len(h_la3d)):
        q_la3d.append(math.sqrt((h_la3d[i])**2 + (k_la3d[i])**2 + (l_la3d[i])**2))
    L_cubic = [la3d, la3d,la3d]
    for i in range(len(h_la3d)):
        if i != 0 and i != 1 and i != 2 and i != 3 and i != 6 and i != 7 and i != 8 and i != 9 and i != 10 and i != 13 and i != 14 and i != 15 and i != 20:
            L_cubic.append(la3d*(q_la3d[i]/q_la3d[4]))
            L_cubic.append(la3d*(q_la3d[i]/q_la3d[4]))
            L_cubic.append(la3d*(q_la3d[i]/q_la3d[4]))
            
    a_la3d = 2*math.pi*math.sqrt(6)/la3d
    
    error_a_la3d = (2*math.pi*math.sqrt(6)/(la3d)**2)*error_q
    
    #lm3m
    q_lm3m = []
    h_lm3m = [1,1,1,2,2,2,2,2,3,3,2,3,3,4,4,4,3,4,3,4,4,3]
    k_lm3m = [0,1,1,0,1,1,2,2,1,1,2,2,2,0,1,1,3,2,3,2,3,3]
    l_lm3m = [0,0,1,0,0,1,0,1,0,1,2,0,1,0,0,1,1,0,2,2,1,3]
    for i in range(len(h_lm3m)):
        q_lm3m.append(math.sqrt((h_lm3m[i])**2 + (k_lm3m[i])**2 + (l_lm3m[i])**2))
    l_cubic = [lm3m, lm3m,lm3m]
    for i in range(1,len(h_lm3m)):
        if i != 2 and i != 4 and i != 7 and i != 9 and i != 11 and i != 14 and i != 16 and i != 21:
            l_cubic.append(lm3m*(q_lm3m[i]/q_lm3m[1]))
            l_cubic.append(lm3m*(q_lm3m[i]/q_lm3m[1]))
            l_cubic.append(lm3m*(q_lm3m[i]/q_lm3m[1]))
            
    a_lm3m = 2*math.pi*math.sqrt(2)/lm3m
    
    error_a_lm3m = (2*math.pi*math.sqrt(2)/(lm3m)**2)*error_q
            
    #Cubic 203, 207
    Cubic_2 = []
    q_Cubic = []
    h_Cubic = [1.,1.,1.,2.,2.,2.,2.,2.,3.,3.,2.,3.,3.,4.,4.,4.,3.,4.,3.,4.,4.,3.]
    k_Cubic = [0,1,1,0,1,1,2,2,1,1,2,2,2,0,1,1,3,2,3,2,3,3]
    l_Cubic = [0,0,1,0,0,1,0,1,0,1,2,0,1,0,0,1,1,0,2,2,1,3]
    for i in range(len(h_Cubic)):
        q_Cubic.append(math.sqrt((h_Cubic[i])**2 + (k_Cubic[i])**2 + (l_Cubic[i])**2))
    Cubic_2 = [Cubic,Cubic,Cubic]
    for i in range(1,len(h_Cubic)):
        Cubic_2.append(Cubic*q_Cubic[i])
        Cubic_2.append(Cubic*q_Cubic[i])
        Cubic_2.append(Cubic*q_Cubic[i])
        
    a_cubic = 2*math.pi/Cubic
    
    error_a_cubic = (2*math.pi/(Cubic)**2)*error_q
    
    #FCC
    q_FCC = []
    h_FCC = [1.,2.,2.,2.,2.,2.,3.,3.,2.,3.,3.,4.,4.,4.,3.,4.,3.,4.,4.,3.]
    k_FCC = [1.,0.,1.,1.,2.,2.,1.,1.,2.,2.,2.,0.,1.,1.,3.,2.,3.,2.,3.,3.]
    l_FCC = [1.,0.,0.,1.,0.,1.,0.,1.,2.,0.,1.,0.,0.,1.,1.,0.,2.,2.,1.,3.]
    q_FCC = [2*math.pi*math.sqrt((h_FCC[0])**2 + (k_FCC[0])**2 + (l_FCC[0])**2)/fcc]
    for i in range(1,len(h_FCC)):
        q_FCC.append(math.sqrt((h_FCC[i])**2 + (k_FCC[i])**2 + (l_FCC[i])**2))
    FCC = [q_FCC[0],q_FCC[0],q_FCC[0]]
    for i in range(len(q_FCC)):
        if  i !=0 and i != 1 and i != 4 and i != 5 and i != 7 and i != 8 and i != 11 and i!= 12 and i != 14 and i != 15 and i != 18 and i != 20:
            FCC.append(FCC[0]*q_FCC[i])
            FCC.append(FCC[0]*q_FCC[i])
            FCC.append(FCC[0]*q_FCC[i])
           
    a_FCC = 2*math.pi*math.sqrt(3)/fcc
    
    error_a_FCC = (2*math.pi*math.sqrt(3)/fcc**2)*error_q
    
    #BCC
    h_BCC = [1.,1.,2.,2.,2.,2.,2.,3.,3.,2.,3.,3.,4.,4.,4.,3.,4.,3.,4.,4.,3.]
    k_BCC = [1.,1.,0.,1.,1.,2.,2.,1.,1.,2.,2.,2.,0.,1.,1.,3.,2.,3.,2.,3.,3.]
    l_BCC = [0.,1.,0.,0.,1.,0.,1.,0.,1.,2.,0.,1.,0.,0.,1.,1.,0.,2.,2.,1.,3.]    
    q_BCC = [2*math.pi*math.sqrt((h_BCC[0])**2 + (k_BCC[0])**2 + (l_BCC[0])**2)/bcc]
    for i in range(1,len(h_BCC)):
        q_BCC.append(math.sqrt((h_BCC[i])**2 + (k_BCC[i])**2 + (l_BCC[i])**2))
    BCC = [q_BCC[0],q_BCC[0],q_BCC[0]]
    for i in range(1,len(q_BCC)):
        if  i != 1 and i != 3 and i != 6 and i != 8 and i != 10 and i != 13 and i!= 15 and i != 20:
            BCC.append(BCC[0]*q_BCC[i])
            BCC.append(BCC[0]*q_BCC[i])
            BCC.append(BCC[0]*q_BCC[i])
            
    a_BCC = 2*math.pi*math.sqrt(2)/bcc
    
    error_a_BCC = (2*math.pi*math.sqrt(2)/bcc**2)*error_q
    
    #HCP
    h_HCP = [1,0,1,1,0,1]
    k_HCP = [0,0,1,1,0,1]
    l_HCP = [0,2,0,2,4,4]
    q_HCP = []
    HCP = []
    for i in range(len(h_HCP)):
        q_HCP.append(2*math.pi*math.sqrt(4*((h_HCP[i])**2 + h_HCP[i]*k_HCP[i] + (k_HCP[i])**2)/(3*(hcp)**2) + (l_HCP[i]/(1.633*hcp))**2))
    for i in range(len(q_HCP)):
        HCP.append(q_HCP[i])
        HCP.append(q_HCP[i])
        HCP.append(q_HCP[i])
    

        
    a_HCP = 2*math.pi/hcp
    c_HCP = 1.633*a_HCP
    
    error_a_HCP = (2*math.pi/hcp**2)*error_q
    error_c_HCP = 1.633*error_q
    
    #Ribbon phase (2D monoclinic - gamma = 90°)
    #h_Rph = [0.,1.,2.,4.,3.,2.,1.,0.,5.,4.,3.,2.,1.,0.,4.,2.,1.]
    #k_Rph = [2.,1.,0.,0.,1.,2.,3.,4.,1.,2.,3.,4.,5.,6.,4.,6.,7.]
    #l_Rph = []
    #RPH = []
    #for i in range(len(h_Rph)):
        #l_Rph.append(2*math.pi*math.sqrt(((h_Rph[i]/rph_a)**2 + (k_Rph[i]/rph_b)**2 - (2*h_Rph[i]*k_Rph[i]*math.cos(math.radians(90)))/(rph_a*rph_b))/(1-((math.cos(math.radians(90)))**2))))
    #for i in range(len(h_Rph)):
    #    RPH.append(l_Rph[i])
    #    RPH.append(l_Rph[i])
    #    RPH.append(l_Rph[i])
    
    
    #Ripple phase
    #h_Rip = [0,1,0,1,1,2,2,2,3,3,3,3,4,4]
    #k_Rip = [1,0,2,1,2,0,1,2,-1,0,1,2,-1,0]
    #l_Rip = []
    #RIP = []
    #for i in range(len(h_Rip)):
    #    l_Rip.append(2*math.pi*math.sqrt(((h_Rip[i]/rip_a)**2 + (k_Rip[i]/rip_b)**2 - (2*h_Rip[i]*k_Rip[i]*math.cos(math.radians(rip_gamma))/(rip_a*rip_b))/(1-((math.cos(math.radians(rip_gamma))))**2))))
    #for i in range(len(h_Rip)):
    #    RIP.append(l_Rip[i])
    #    RIP.append(l_Rip[i]) 
    #    RIP.append(l_Rip[i]) 
    
    ##############################################################################################
    # It Identifies peaks near calculated parameters
    ############################################################################################## 
    l1 = 0
    q_SL1 = []
    X_SL1 = []                                    
    for i in range(0,len(SL1)):            
        for j in range(2,len(q_max)):
            if SL1[i] >= (1-var)*q_max[j] and SL1[i] < (1+var)*q_max[j] and q_max[j] not in q_SL1: 
                l1 += 1
                q_SL1.append(q_max[j])
                X_SL1.append(X[X_max[j]])
    
    l2 = 0
    q_SL2 = []
    X_SL2 = []                                    
    for i in range(0,len(SL2)):            
        for j in range(2,len(q_max)):
            if SL2[i] >= (1-var)*q_max[j] and SL2[i] < (1+var)*q_max[j] and q_max[j] not in q_SL2: 
                l2 += 1
                q_SL2.append(q_max[j])
                X_SL2.append(X[X_max[j]])
    
    hexa = 0
    q_SIH =[]
    X_SIH = []            
    for i in range(0,len(SIH)):            
        for j in range(2,len(q_max)):
            if SIH[i] >= (1-var)*q_max[j] and SIH[i] < (1+var)*q_max[j] and q_max[j] not in q_SIH:               
                hexa += 1
                q_SIH.append(q_max[j])
                X_SIH.append(X[X_max[j]])
    mn_cubic = 0
    q_Pm3n_cubic = []
    X_Pm3n_cubic = []
    for i in range(0,len(Pm3n_cubic)):            
        for j in range(2,len(q_max)):
            if Pm3n_cubic[i] >= (1-var)*q_max[j] and Pm3n_cubic[i] < (1+var)*q_max[j] and q_max[j] not in q_Pm3n_cubic:                
                mn_cubic += 1
                q_Pm3n_cubic.append(q_max[j])
                X_Pm3n_cubic.append(X[X_max[j]])                       
    Lcubic = 0
    q_L_cubic = []
    X_L_cubic = []            
    for i in range(0,len(L_cubic)):            
        for j in range(2,len(q_max)):
            if L_cubic[i] >= (1-var)*q_max[j] and L_cubic[i] < (1+var)*q_max[j] and q_max[j] not in q_L_cubic:            
                Lcubic += 1
                q_L_cubic.append(q_max[j])
                X_L_cubic.append(X[X_max[j]])
    lcubic = 0 
    q_l_cubic = []
    X_l_cubic = []           
    for i in range(0,len(l_cubic)):            
        for j in range(2,len(q_max)):
            if l_cubic[i] >= (1-var)*q_max[j] and l_cubic[i] < (1+var)*q_max[j] and q_max[j] not in q_l_cubic:            
                lcubic += 1
                q_l_cubic.append(q_max[j])
                X_l_cubic.append(X[X_max[j]])
    Pcubic = 0
    q_P_cubic = []
    X_P_cubic = []            
    for i in range(0,len(P_cubic)):            
        for j in range(2,len(q_max)):
            if P_cubic[i] >= (1-var)*q_max[j] and P_cubic[i] < (1+var)*q_max[j] and q_max[j] not in q_P_cubic:            
                Pcubic += 1
                q_P_cubic.append(q_max[j])
                X_P_cubic.append(X[X_max[j]])
    mm_cubic = 0
    q_Fm3m_cubic = []
    X_Fm3m_cubic = []
    for i in range(0,len(Fm3m_cubic)):            
        for j in range(2,len(q_max)):
            if Fm3m_cubic[i] >= (1-var)*q_max[j] and Fm3m_cubic[i] < (1+var)*q_max[j] and q_max[j] not in q_Fm3m_cubic:            
                mm_cubic += 1
                q_Fm3m_cubic.append(q_max[j])
                X_Fm3m_cubic.append(X[X_max[j]])
    dm_cubic = 0
    q_Fd3m_cubic = []
    X_Fd3m_cubic = []            
    for i in range(0,len(Fd3m_cubic)):            
        for j in range(2,len(q_max)):
            if Fd3m_cubic[i] >= (1-var)*q_max[j] and Fd3m_cubic[i] < (1+var)*q_max[j] and q_max[j] not in q_Fd3m_cubic:                 
                dm_cubic += 1
                q_Fd3m_cubic.append(q_max[j])
                X_Fd3m_cubic.append(X[X_max[j]]) 
    cubic2 = 0
    q_Cubic_2 = []
    X_Cubic_2 = []                                                      
    for i in range(0,len(Cubic_2)):                       
        for j in range(2,len(q_max)):
            if Cubic_2[i] > (1-var)*q_max[j] and Cubic_2[i] < (1+var)*q_max[j] and q_max[j] not in q_Cubic_2:            
                cubic2 += 1
                q_Cubic_2.append(q_max[j])
                X_Cubic_2.append(X[X_max[j]])
    face = 0
    q_Face = []
    X_Face = []
    for i in range(0,len(FCC)):            
        for j in range(2,len(q_max)):
            if FCC[i] >= (1-var)*q_max[j] and FCC[i] < (1+var)*q_max[j] and q_max[j] not in q_Face:            
                face += 1
                q_Face.append(q_max[j])
                X_Face.append(X[X_max[j]])
                
    body = 0
    q_Body = []
    X_Body = []
    for i in range(0,len(BCC)):            
        for j in range(2,len(q_max)):
            if BCC[i] >= (1-var)*q_max[j] and BCC[i] < (1+var)*q_max[j] and q_max[j] not in q_Body:            
                body += 1
                q_Body.append(q_max[j])
                X_Body.append(X[X_max[j]])
    cp = 0
    q_hcp = []
    X_hcp = []
    for i in range(0,len(HCP)):            
        for j in range(2,len(q_max)):
            if HCP[i] >= (1-var)*q_max[j] and HCP[i] < (1+var)*q_max[j] and q_max[j] not in q_hcp:            
                cp += 1
                q_hcp.append(q_max[j])
                X_hcp.append(X[X_max[j]])
    #ripple = 0
    #q_ripple = []
    #X_ripple = []
    #for i in range(0,len(RIP)):            
    #    for j in range(2,len(q_max)):
    #        if RIP[i] >= (1-var)*q_max[j] and RIP[i] < (1+var)*q_max[j] and q_max[j] not in q_ripple:            
    #            ripple += 1 
    #            q_ripple.append(q_max[j])
    #            X_ripple.append(X[X_max[j]])
                                       
    ############################################################################
    # Plots
    # 
    ############################################################################
    plt.figure()
    plt.plot(q,X,'r-')
    plt.plot(x,Gaussian(x,*popt),'k-',label = 'Gaussian fit')    #Plot the gaussian curve 
    dimensions = []
    print(SL1)
    #Put q values                                                                      
    for i in range(0,len(SL1)):            
        for j in range(2,len(q_max)):
            if SL1[i] >= (1-var)*q_max[j] and SL1[i] < (1+var)*q_max[j]:                              
                plt.text(q_max[j], X[X_max[j]], r'q = '+str(q_max[j]), fontsize = 10)           
     
    for i in range(0,len(SL2)):            
        for j in range(2,len(q_max)):
            if SL2[i] >= (1-var)*q_max[j] and SL2[i] < (1+var)*q_max[j]:                              
                plt.text(q_max[j], X[X_max[j]], r'q = '+str(q_max[j]), fontsize = 10)
               
    for i in range(0,len(SIH)):            
        for j in range(2,len(q_max)):
            if SIH[i] >= (1-var)*q_max[j] and SIH[i] < (1+var)*q_max[j]:
                plt.text(q_max[j], X[X_max[j]], r'q = '+str(q_max[j]), fontsize = 10)
                
    
    for i in range(0,len(Pm3n_cubic)):            
        for j in range(2,len(q_max)):
            if Pm3n_cubic[i] >= (1-var)*q_max[j] and Pm3n_cubic[i] < (1+var)*q_max[j]:                         
                plt.text(q_max[j], X[X_max[j]], r'q = '+str(q_max[j]), fontsize = 10) 
                                   
                
    for i in range(0,len(L_cubic)):            
        for j in range(2,len(q_max)):
            if L_cubic[i] >= (1-var)*q_max[j] and L_cubic[i] < (1+var)*q_max[j]:                                
                plt.text(q_max[j], X[X_max[j]], r'q = '+str(q_max[j]), fontsize = 10)
             
    for i in range(0,len(l_cubic)):            
        for j in range(2,len(q_max)):
            if l_cubic[i] >= (1-var)*q_max[j] and l_cubic[i] < (1+var)*q_max[j]:                            
                plt.text(q_max[j], X[X_max[j]], r'q = '+str(q_max[j]), fontsize = 10)
                
    for i in range(0,len(P_cubic)):            
        for j in range(2,len(q_max)):
            if P_cubic[i] >= (1-var)*q_max[j] and P_cubic[i] < (1+var)*q_max[j]:                            
                plt.text(q_max[j], X[X_max[j]], r'q = '+str(q_max[j]), fontsize = 10)
    
    for i in range(0,len(Fm3m_cubic)):            
        for j in range(2,len(q_max)):
            if Fm3m_cubic[i] >= (1-var)*q_max[j] and Fm3m_cubic[i] < (1+var)*q_max[j]:            
                plt.text(q_max[j], X[X_max[j]], r'q = '+str(q_max[j]), fontsize = 10)
                
    for i in range(0,len(Fd3m_cubic)):            
        for j in range(2,len(q_max)):
            if Fd3m_cubic[i] >= (1-var)*q_max[j] and Fd3m_cubic[i] < (1+var)*q_max[j]:                           
                plt.text(q_max[j], X[X_max[j]], r'q = '+str(q_max[j]), fontsize = 10)
                                                       
    for i in range(0,len(Cubic_2)):                       
        for j in range(2,len(q_max)):
            if Cubic_2[i] > (1-var)*q_max[j] and Cubic_2[i] < (1+var)*q_max[j]:                                
                plt.text(q_max[j], X[X_max[j]], r'q = '+str(q_max[j]), fontsize = 10)
    
    for i in range(0,len(FCC)):            
        for j in range(2,len(q_max)):
            if FCC[i] >= (1-var)*q_max[j] and FCC[i] < (1+var)*q_max[j]:            
                plt.text(q_max[j], X[X_max[j]], r'q = '+str(q_max[j]), fontsize = 10)
    
    for i in range(0,len(BCC)):            
        for j in range(2,len(q_max)):
            if BCC[i] >= (1-var)*q_max[j] and BCC[i] < (1+var)*q_max[j]:            
                plt.text(q_max[j], X[X_max[j]], r'q = '+str(q_max[j]), fontsize = 10)
    
    for i in range(0,len(HCP)):            
        for j in range(2,len(q_max)):
            if HCP[i] >= (1-var)*q_max[j] and HCP[i] < (1+var)*q_max[j]:            
                plt.text(q_max[j], X[X_max[j]], r'q = '+str(q_max[j]), fontsize = 10)
    
    #for i in range(0,len(RIP)):            
    #    for j in range(2,len(q_max)):
    #        if RIP[i] >= (1-var)*q_max[j] and RIP[i] < (1+var)*q_max[j]:            
    #            plt.text(q_max[j], X[X_max[j]], r'q = '+str(q_max[j]), fontsize = 10)
    
    times = 3                      
    if l1 > times:
        plt.plot(q_max[first],X[X_max[first]], marker='|',markersize = 70,linewidth=4, color = 'm',label='lammelar 1')
        plt.text(q_max[first],X[X_max[first]], r'q = '+str(round(q_max[first],2)), fontsize = 10)
        print('Lammelar 1 phase dimension (nm):')
        print('d10 = ' + str(round(d10_SL1,2)) + ' (' + str(round(error_d10_SL1,2)) + ')')
        print(' ')
        dimensions.append('Lammelar 1 phase dimension (nm):\n')
        dimensions.append('d10 = ' + str(round(d10_SL1,2)) + ' (' + str(round(error_d10_SL1,2)) + ')\n')
        dimensions.append(' \n')
        for i in range(len(q_SL1)):
            plt.plot(q_SL1[i],X_SL1[i],marker='|',markersize = 70, linewidth=4,color='m')
    if l2 > times:
        plt.plot(q_max[first],X[X_max[first]], marker='|',markersize = 70,linewidth=4, color = 'm',label='lammelar 2')
        plt.text(q_max[first],X[X_max[first]], r'q = '+str(round(q_max[first],2)), fontsize = 10)
        print('Lammelar 2 phase dimension (nm):')
        print('d10 = ' + str(round(d10_SL2,2)) + ' (' + str(round(error_d10_SL2,2)) + ')')
        print(' ')
        dimensions.append('Lammelar 2 phase dimension (nm):\n')
        dimensions.append('d10 = ' + str(round(d10_SL2,2)) + ' (' + str(round(error_d10_SL2,2)) + ')\n')
        dimensions.append('\n ')
        for i in range(len(q_SL2)):
            plt.plot(q_SL2[i],X_SL2[i],marker='|',markersize = 70, linewidth=4,color='m')
    if hexa > times:
        plt.plot(q_max[first],X[X_max[first]], marker='|',markersize = 70,linewidth=4, color = 'r',label='hexagonal main peak')
        plt.text(q_max[first],X[X_max[first]], r'q = '+str(round(q_max[first],2)), fontsize = 10)
        print('Hexagonal phase dimensions (nm):')
        print('d10 = ' + str(round(d10_hmp,2)) + ' (' + str(round(error_d10_hmp,2)) + ')')
        print('d11 = ' + str(round(d11_hmp,2)) + ' (' + str(round(error_d11_hmp,2)) + ')')
        print(' ')
        dimensions.append('Hexagonal phase dimensions (nm):\n')
        dimensions.append('d10 = ' + str(round(d10_hmp,2)) + ' (' + str(round(error_d10_hmp,2)) + ')\n')
        dimensions.append('d11 = ' + str(round(d11_hmp,2)) + ' (' + str(round(error_d11_hmp,2)) + ')\n')
        dimensions.append('\n ')
        for i in range(len(q_SIH)):
            plt.plot(q_SIH[i],X_SIH[i],marker='|',markersize = 70,linewidth=4,color='r')
    if mn_cubic > times:
        plt.plot(q_max[first],X[X_max[first]], marker='|',markersize = 70,linewidth=4, color = 'B',label='Pm3n cubic')
        plt.text(q_max[first],X[X_max[first]], r'q = '+str(round(q_max[first],2)), fontsize = 10)        
        print('Pm3n cubic dimensions (nm):')
        print('a = ' + str(round(a_Pm3n,2)) + ' (' + str(round(error_a_Pm3n,2))+ ')')
        print(' ')
        dimensions.append('Pm3n cubic dimensions (nm):\n')
        dimensions.append('a = ' + str(round(a_Pm3n,2)) + ' (' + str(round(error_a_Pm3n,2))+ ')\n')
        dimensions.append(' \n')
        for i in range(len(q_Pm3n_cubic)):
            plt.plot(q_Pm3n_cubic[i],X_Pm3n_cubic[i],marker='|',markersize = 70,linewidth=4,color='B')            
    if Lcubic > times:
        plt.plot(q_max[first],X[X_max[first]], marker='|',linewidth=4,markersize = 70, color = 'b',label='la3d cubic')        
        plt.text(q_max[first],X[X_max[first]], r'q = '+str(round(q_max[first],2)), fontsize = 10)
        print('la3d cubic dimensions (nm):')
        print('a = ' + str(round(a_la3d,2)) + ' (' + str(round(error_a_la3d,2)) + ')')
        print(' ')
        dimensions.append('la3d cubic dimensions (nm):\n')
        dimensions.append('a = ' + str(round(a_la3d,2)) + ' (' + str(round(error_a_la3d,2)) + ')\n')
        dimensions.append('\n ')
        for i in range(len(q_L_cubic)):
            plt.plot(q_L_cubic[i],X_L_cubic[i],marker='|', markersize = 70,linewidth=4,color='b')            
    if lcubic > times:
        plt.plot(q_max[first],X[X_max[first]], marker='|',markersize = 70,linewidth=4, color = 'k',label='lm3m cubic')        
        plt.text(q_max[first],X[X_max[first]], r'q = '+str(round(q_max[first],2)), fontsize = 10)
        print('lm3m cubic dimension (nm):')
        print('a = ' + str(round(a_lm3m,2)) + ' (' + str(round(error_a_lm3m,2)) + ')')
        print(' ')
        dimensions.append('lm3m cubic dimension (nm):\n')
        dimensions.append('a = ' + str(round(a_lm3m,2)) + ' (' + str(round(error_a_lm3m,2)) + ')\n')
        dimensions.append(' \n')
        for i in range(len(q_l_cubic)):
            plt.plot(q_l_cubic[i],X_l_cubic[i],marker='|',markersize = 70, linewidth=4,color='k')            
    if Pcubic > times:
        plt.plot(q_max[first],X[X_max[first]], marker='|',markersize = 70,linewidth=4, color = 'c',label='Pn3m cubic')        
        plt.text(q_max[first],X[X_max[first]], r'q = '+str(round(q_max[first],2)), fontsize = 10)
        print('Pn3m cubic dimension (nm):')
        print('a = ' + str(round(a_Pn3m,2)) + ' (' + str(round(error_a_Pn3m,2)) + ')')
        print(' ')
        dimensions.append('Pn3m cubic dimension (nm):\n')
        dimensions.append('a = ' + str(round(a_Pn3m,2)) + ' (' + str(round(error_a_Pn3m,2)) + ')\n')
        dimensions.append('\n ')
        for i in range(len(q_P_cubic)):
            plt.plot(q_P_cubic[i],X_P_cubic[i],marker='|',markersize = 70, linewidth=4,color='c')
    if mm_cubic > times:
        plt.plot(q_max[first],X[X_max[first]], marker='|',markersize = 70,linewidth=4, color = 'g',label='Fm3m cubic') 
        plt.text(q_max[first],X[X_max[first]], r'q = '+str(round(q_max[first],2)), fontsize = 10)       
        print('Fm3m cubic dimension (nm):')
        print('a = ' + str(round(a_Fm3m,2)) + ' (' + str(round(error_a_Fm3m,2)) + ')')
        print(' ')
        dimensions.append('Fm3m cubic dimension (nm):\n')
        dimensions.append('a = ' + str(round(a_Fm3m,2)) + ' (' + str(round(error_a_Fm3m,2)) + ')\n')
        dimensions.append('\n ')
        for i in range(len(q_Fm3m_cubic)):
            plt.plot(q_Fm3m_cubic[i],X_Fm3m_cubic[i],marker='|',markersize = 70, linewidth=4,color='g')
    if dm_cubic > times:
        plt.plot(q_max[first],X[X_max[first]], marker='|',markersize = 70,linewidth=4, color = 'M',label='Fd3m cubic') 
        plt.text(q_max[first],X[X_max[first]], r'q = '+str(round(q_max[first],2)), fontsize = 10)       
        print('Fd3m cubic dimension (nm):')
        print('a = ' + str(round(a_Fd3m,2)) + ' (' + str(round(error_a_Fd3m,2)) + ')')
        print(' ')
        dimensions.append('Fd3m cubic dimension (nm):\n')
        dimensions.append('a = ' + str(round(a_Fd3m,2)) + ' (' + str(round(error_a_Fd3m,2)) + ')\n')
        dimensions.append(' \n')
        for i in range(len(q_Fd3m_cubic)):
            plt.plot(q_Fd3m_cubic[i],X_Fd3m_cubic[i],marker='|',markersize = 70, linewidth=4,color='M')
    if cubic2 > times:
        plt.plot(q_max[first],X[X_max[first]], marker='|',markersize = 70,linewidth=4, color = 'y',label='cubic')
        plt.text(q_max[first],X[X_max[first]], r'q = '+str(round(q_max[first],2)), fontsize = 10)        
        print('Cubic dimension (nm):')
        print('a = ' + str(round(a_cubic,2)) + ' (' + str(round(error_a_cubic,2)) + ')')
        print(' ')
        dimensions.append('Cubic dimension (nm):\n')
        dimensions.append('a = ' + str(round(a_cubic,2)) + ' (' + str(round(error_a_cubic,2)) + ')\n')
        dimensions.append('\n ')
        for i in range(len(q_Cubic_2)):
            plt.plot(q_Cubic_2[i],X_Cubic_2[i],marker='|',markersize = 70, linewidth=4,color='y')
    if face > times:
        plt.plot(q_max[first],X[X_max[first]], marker='|',markersize = 70,linewidth=4, color = 'R',label='FCC') 
        plt.text(q_max[first],X[X_max[first]], r'q = '+str(round(q_max[first],2)), fontsize = 10)       
        print('FCC dimension (nm):')
        print('a = ' + str(round(a_FCC,2)) + ' (' +str(round(error_a_FCC,2)) + ')')
        print(' ')
        dimensions.append('FCC dimension (nm):\n')
        dimensions.append('a = ' + str(round(a_FCC,2)) + ' (' +str(round(error_a_FCC,2)) + ')\n')
        dimensions.append('\n ')
        for i in range(len(q_Face)):
            plt.plot(q_Face[i],X_Face[i],marker='|',markersize = 70, linewidth=4,color='R')
    if body > times:
        plt.plot(q_max[first],X[X_max[first]], marker='|',markersize = 60,linewidth=4, color = 'C',label='BCC')
        plt.text(q_max[first],X[X_max[first]], r'q = '+str(round(q_max[first],2)), fontsize = 10)
        print('BCC dimension (nm):')
        print('a = ' + str(round(a_BCC,2)) + ' (' + str(round(error_a_BCC,2)) + ')')
        print(' ')
        dimensions.append('BCC dimension (nm):\n')
        dimensions.append('a = ' + str(round(a_BCC,2)) + ' (' + str(round(error_a_BCC,2)) + ')\n')
        dimensions.append('\n ')
        for i in range(len(q_Body)):
            plt.plot(q_Body[i],X_Body[i],marker='|',markersize = 70, linewidth=4,color='C')
    if cp > times:
        plt.plot(q_max[first],X[X_max[first]], marker='|',markersize = 70,linewidth=4, color = 'G',label='HCP')
        plt.text(q_max[first],X[X_max[first]], r'q = '+str(round(q_max[first],2)), fontsize = 10)
        print('HCP dimension (nm):')
        print('a = ' + str(round(a_HCP,2)) + ' (' + str(round(error_a_HCP,2)) + ')')
        print('c = ' + str(round(c_HCP,2)) + str(round(error_c_HCP,2)) + ')')
        print(' ')
        dimensions.append('HCP dimension (nm):\n')
        dimensions.append('a = ' + str(round(a_HCP,2)) + ' (' + str(round(error_a_HCP,2)) + ')\n')
        dimensions.append('c = ' + str(round(c_HCP,2)) + str(round(error_c_HCP,2)) + ')\n')
        dimensions.append(' \n')
        for i in range(len(q_hcp)):
            plt.plot(q_hcp[i],X_hcp[i],marker='|',markersize = 70, linewidth=4,color='G')
    #if ripple > times:
    #    plt.plot(q_max[1],X[X_max[1]], marker='|',markersize = 70,linewidth=4, color = 'Y',label='ripple phase')
    #    plt.text(q_max[1],X[X_max[1]], r'q = '+str(round(qMax,2)) + ' (' + str(round(error_q,2)) + ')', fontsize = 10)       
    #    for i in range(len(q_ripple)):
    #        plt.plot(q_ripple[i],X_ripple[i],marker='|',markersize = 70, linewidth=4,color='Y')
    plt.legend(fontsize=10)       
    #plt.xlim(0,1.)
    #for i in range(len(dimensions)):
        #plt.text(q[len(q)-50], min(X)*(-0.3*i + 2.8),dimensions[i], fontsize = 10) 
    #plt.title('PL407 30$\%$ + TR 46.12245$^\circ$C')
    plt.title('P30 + IPM + LEC 40 $^\circ$C')
    #plt.title('P30 + IPM + LDC + LEC 25$^\circ$C')
    plt.xlabel('q (A$^-$$^1$)')
    plt.ylabel('Intensity (arbitrary units)')
    plt.xscale('log')
    plt.yscale('log')
    plt.show()
    
    #plt.figure()
    #plt.plot(q,X,'r-')
    #if l1 != 0:
    #    for i in range(len(SL1)):           
    #        plt.axvline(x = SL1[i], color = 'g')            
    #if l2 != 0:
    #    for i in range(len(SL2)):
    #        plt.axvline(x = SL2[i], color = 'm')
    #if hexa != 0:
    #    for i in range(len(SIH)):
    #        plt.axvline(x = SIH[i], color = 'r')
    #if Ori != 0:
    #    for i in range(len(OCA)):
    #        plt.axvline(x = OCA[i], color = 'k')
    #if la3d != 0:
    #    for i in range(len(L_cubic)):
    #        plt.axvline(x = L_cubic[i], color = 'b')
    #if lm3m != 0:
    #    for i in range(len(l_cubic)):
    #        plt.axvline(x = l_cubic[i], color = 'k')
    #if Pn3m != 0:
    #    for i in range(len(P_cubic)):
    #        plt.axvline(x = P_cubic[i], color = 'c')
    #if Cubic != 0:
    #    for i in range(len(Cubic_2)):
    #        plt.axvline(x = Cubic_2[i], color = 'y')    
    #plt.title('SAXS')
    #plt.xlabel('q (nm$^-$$^1$)')
    #plt.ylabel('Intensity (arbitrary values)')
    #plt.show()
main()