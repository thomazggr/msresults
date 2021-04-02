def dorev(it):
    if it / abs(it) == 1:
        return -abs(it)
    elif it / abs(it) == -1:
        return abs(it)


def rvs(val):

    if type(val) == list:
        ingresp = []
        for it in val:
            ingresp.append(dorev(it))
        return ingresp
    elif type(val) == str:
        return
    else:
        return dorev(val)