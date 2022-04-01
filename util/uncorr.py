#!/ usr / bin / env python

def ind(data):
    #Import packages
    import ruptures as rpt 
    import numpy as np
    from statistics import stdev

    #Convert data to float
    raw = np.zeros(len(data))
    for i in range(len(data)):
        raw[i] = float(data[i])

    #Apply ruptures to find uncorrelated samples
    model = 'l1'
    algo = rpt.Binseg(model=model, min_size=10, jump=10).fit(raw)
    n = len(raw)
    sigma = stdev(raw)
    t_uncorr = algo.predict(epsilon=3 * n * sigma ** 2)

    return t_uncorr

def sort(data, t_uncorr):
    #import packages
    import numpy as np

    #Convert data to float
    raw = np.zeros(len(data))
    for i in range(len(data)):
        raw[i] = float(data[i])

    #Reduce to uncorrelated data
    num=len(t_uncorr)
    data_uncorr = np.zeros(num)
    n = 0
    for i in range(len(raw)):
        if i in t_uncorr:
            data_uncorr[n] = raw[i]
            n += 1

    return data_uncorr

def char(data, t_uncorr):
    #Reduce to uncorrelated data
    num=len(t_uncorr)
    data_uncorr = []
    for i in range(len(data)):
        if i in t_uncorr:
            data_uncorr.append(data[i])

    return data_uncorr

