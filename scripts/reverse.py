import pandas as pd


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
    elif type(val) == type(pd.DataFrame()):
        df = val

        for i in df.index:
            for j in df.columns:
                if df.at[i, j] / abs(df.at[i, j]) == 1:
                    df.at[i, j] = -abs(df.at[i, j])

                elif df.at[i, j] / abs(df.at[i, j]) == -1:
                    df.at[i, j] = abs(df.at[i, j])

    elif type(val) == str:
        return
    else:
        return dorev(val)