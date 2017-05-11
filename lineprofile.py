import numpy as np





# Computes lineprofile
def lineprofile(fluxc, redsc):
    print "Computing line profile..."

    xarr = np.array([])
    yarr = np.array([])

    Ny_dense, Nx_dense = np.shape(fluxc)

    energ = 1.0
    for jj in range(Ny_dense):
        for ii in range(Nx_dense):
            fy  = fluxc[jj, ii]
            xii = redsc[jj, ii]

            if xii > 0.0:
                xarr = np.append(xarr, xii*energ)
                yarr = np.append(yarr, fy)


    xind = np.argsort(xarr)
    xarrs = xarr[xind]
    yarrs = yarr[xind]

    NN = len(xarrs)

    emin = np.min(xarrs)*0.99
    emax = np.max(xarrs)*1.01

    Nr = 100
    es = np.linspace(emin, emax, Nr)
    yy2 = np.zeros((Nr))

    xst = 0
    for ii in range(1,Nr):
        for jj in range(xst, NN):
            if es[ii-1] <= xarrs[jj] < es[ii]:
                yy2[ii] += yarrs[jj]
            elif xarrs[jj] >= es[ii]:
                xst = jj
                break

    #normalize
    des = np.diff(es)[1]
    yy2 = yy2 / np.sum(yy2*des)

        
    return es, yy2         




