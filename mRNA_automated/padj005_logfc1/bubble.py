import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

kegg_data = pd.read_csv("results__kegg_full_table.txt", sep="\t")

print(kegg_data.head())
# quit()

sns.set_context("talk", font_scale=1.1)
plt.figure(figsize=(20,12))
sns.scatterplot(x="p.adjust", 
                y="ID",
                size="Count",
                sizes=(100,900),
                alpha=0.5,
                data=kegg_data)
# Put the legend out of the figure
plt.legend(bbox_to_anchor=(1.01, 1),borderaxespad=0)
# Put the legend out of the figure
#plt.legend(bbox_to_anchor=(1.01, 0.54),  borderaxespad=0.)
plt.xlabel("Adjusted P Value")
plt.ylabel("Pathway")
plt.title("KEGG Pathways")
plt.tight_layout()
plt.savefig("Bubble_plot_size_range_Seaborn_scatterplot.png",
                    format='png',dpi=150)