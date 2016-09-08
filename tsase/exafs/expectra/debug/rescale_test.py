#!/usr/bin/env python


import numpy

def main():

    exp_chi_file='exp_chi.dat'
#    try:
    k_exp, chi_exp = read_chi(exp_chi_file) 
#    except:
#        k_exp, chi_exp = load_chi_dat(exp_chi_file)

    chi='chi.dat'
    try:
        k, chi = read_chi(chi) 
    except:
        k, chi = load_chi_dat(chi)

    dk   = k_exp[1] - k_exp[0]
    kmin =  0.8
    kmax =  10.0

    k_exp, chi_exp = rescale_chi_calc(k, chi_exp, k_exp, kmax)

#    window = hanning_window_origin(k_exp, kmin, kmax, dk)
#    chi_exp *= window
#    window = hanning_window_origin(k_exp, kmin, kmax, dk)
#    chi *= window

    f = open('test_1.dat', 'w')
    for i in range(len(k_exp)):
        f.write("%6.3f %16.8e\n" % (k_exp[i], chi_exp[i]))
    f.close()

def rescale_chi_calc(k_std, chi_src, k_src, kmax):
    """
    k_std..........k values used as a standard for the rescaling
    chi_src........chi values required to be rescaled
    k_src..........k values corresponding to chi_src
    """
    k_temp = []
    chi_temp = []
    #reset chi_calc based on k_exp
    #tell if k_exp starts from a smaller value
#     try:
#         result = compareValue(k_exp[0],k_cacl[0])
#     except MyValidationError as exception:
#         print exception.message
    i = 0   
    while ( 0 <= i < len(k_std) and k_std[i] <= kmax):
        for j in range(1,len(k_src)):
            if k_src[j-1] < k_std[i] and k_std[i] < k_src[j]:
                chi_temp.append(numpy.interp(k_std[i],
                                           [k_src[j-1],k_src[j]],
                                           [chi_src[j-1],chi_src[j]]))
                k_temp.append(k_std[i])

            elif k_std[i] == k_src[j-1]:
                chi_temp.append(chi_src[j-1])
                k_temp.append(k_std[i])
        i += 1
    return k_temp, chi_temp

def read_chi(filename):
    f = open(filename)
    ks = []
    chis = []
    for line in f:
        if line.startswith('#'):
            continue
        fields = [ float(field) for field in line.split() ]
        k = fields[0]                                                           
        chi = fields[1]
        ks.append(k)
        chis.append(chi)
    f.close()
    ks = numpy.array(ks)
    chis = numpy.array(chis)
    return ks, chis

def load_chi_dat(filename):
    f = open(filename)
    chi_section = False
    k = []
    chi = []
    for line in f:
        line = line.strip()
        if len(line) == 0:
            continue

        fields = line.split()
        if fields[0] == "k" and fields[1] == "chi" and fields[2] == "mag":
            chi_section = True
            continue

        if chi_section:
            k.append(float(fields[0]))
            chi.append(float(fields[1]))
    return numpy.array(k), numpy.array(chi)

def hanning_window(k, kmin, kmax, dk):
    condlist = []
    condlist.append((kmin-dk/2.0 <= k) & (k <  kmin+dk/2.0))
    condlist.append((kmin+dk/2.0 <= k) & (k <= kmax-dk/2.0))
    condlist.append((kmax-dk/2.0 <= k) & (k <= kmax+dk/2.0))

    funclist = []
    funclist.append(lambda x: numpy.sin(numpy.pi*(x-kmin+dk/2.0)/(2.0*dk))**2.0)
    funclist.append(1.0)
    funclist.append(lambda x: numpy.cos(numpy.pi*(x-kmax+dk/2.0)/(2.0*dk))**2.0)

    return numpy.piecewise(k, condlist, funclist) 

def hanning_window_1(k, kmin, kmax, dk):
    condlist = []
    condlist.append((kmin-dk/2.0 <= k) & (k <  kmin+dk/2.0))
    condlist.append((kmin+dk/2.0 <= k) & (k <= kmax-dk/2.0))
    condlist.append((kmax-dk/2.0 <= k) & (k <= kmax+dk/2.0))

    funclist = []
#    funclist.append(lambda x: x)
    funclist.append(1.0)
    funclist.append(1.0)
    funclist.append(1.0)
    funclist.append(0.0)

    return numpy.piecewise(k, condlist, funclist) 

if __name__ == '__main__':
    main()


