import json
from collections import Counter

import pandas as pd

gses = ["GSE33814", "GSE89632", "GSE164760"]

final = []

for i in gses:
    print(i)
    q = "`adj.p.val` <= 0.05 and (`logfc` >= 1 or `logfc` <= -1)"
    t_t = pd.read_csv(f"mRNA_automated/{i}/{i}.tsv", sep="\t")
    lowered = [value.lower() for value in t_t.columns]
    t_t.columns = lowered
    symbols = list(set(list(t_t.query(q)["gene.symbol"])))
    final.extend(symbols)
    print(len(symbols))

Counter(final)

counts = Counter(final)

with open("sample.json", "w") as outfile:
    json.dump(counts, outfile)

for x, z in zip(counts.keys(), counts.values()):
    if z >= 3:
        print(x)
    else:
        pass
