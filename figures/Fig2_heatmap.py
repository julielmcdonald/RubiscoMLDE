# importing the modules 
import pandas as pd
import numpy as np 
import seaborn as sns
import matplotlib.pyplot as plt 
from matplotlib.colors import ListedColormap
from matplotlib.colors import LinearSegmentedColormap

#INPUTS - enter prediction file, AA sequences here
#rbcL
df_US_L = pd.read_csv('/enrichment.csv')

rbcL_seq='MSPQTETKASVGFKAGVKEYKLTYYTPEYQTKDTDILAAFRVTPQPGVPPEEAGAAVAAESSTGTWTTVWTDGLTSLDRYKGRCYRIERVVGEKDQYIAYVAYPLDLFEEGSVTNMFTSIVGNVFGFKALRALRLEDLRIPPAYVKTFQGPPHGIQVERDKLNKYGRPLLGCTIKPKLGLSAKNYGRAVYECLRGGLDFTKDDENVNSQPFMRWRDRFLFCAEALYKAQAETGEIKGHYLNATAGTCEEMIKRAVFARELGVPIVMHDYLTGGFTANTSLAHYCRDNGLLLHIHRAMHAVIDRQKNHGIHFRVLAKALRMSGGDHIHSGTVVGKLEGERDITLGFVDLLRDDFVEQDRSRGIYFTQDWASLPGVLPVASGGIHVWHMPALTEIFGDDSVLQFGGGTLGHPWGNAPGAVANRVALEACVKARNEGRDLAQEGNEIIREACKWSPELAAACEVWKEIVFNFAAVDVLDK*'

AAs = [
    "P",
    "W",
    "F",
    "Y",
    "I",
    "L",
    "M",
    "V",
    "A",
    "G",
    "S",
    "T",
    "C",
    "N",
    "Q",
    "D",
    "E",
    "K",
    "R",
    "H",
]

def AA_numbering(string):
    amino_acid_sequence = string

    numbered_list = [f'{aa}{i}' for i, aa in enumerate(amino_acid_sequence, start=1)]

    return(numbered_list)

# Delete the rows in enrichment df where score is False
mask = df_US_L['Unnamed: 3'] == "FALSE"
df_US_L = df_US_L[~mask]


# Separate into large and small subunits
RbcL_US_df = df_US_L[df_US_L['selection'] == 'NtRbcL']

#LSU - no need to change code

def get_log_likelihood(mutation,df):
    result = df.loc[df["Unnamed: 1"] == mutation, "Unnamed: 2"]
    return result.values[0] if not result.empty else None

def LSU_heatmap(df,output_name):

    #unsupervised
    start_df = df['Unnamed: 1']
    start_df = start_df.str.slice(start=0,stop=1)
    end_df = df['Unnamed: 1']
    end_df = end_df.str.slice(start=-1)
    score = df['Unnamed: 2']

    columns=AA_numbering(rbcL_seq)
    rows=AAs

    # Create a DataFrame filled with NaN values
    plot_df = pd.DataFrame(np.nan, index=rows, columns=columns)

    #Value = 1 for each predicted mutation
    for string in df['Unnamed: 1']:
        starting_AA=string[0:-1]
        columns+=[starting_AA]
        ending_AA=string[-1]
        for col in plot_df.columns:
            if starting_AA==col:
                for index, row in plot_df.iterrows():
                    if index==ending_AA:
                        plot_df.loc[index,col]=float(get_log_likelihood(string,df_US_L))

    #plot_df.fillna(0, inplace=True)

    plot_df.to_csv(output_name)
    
    return (plot_df)

L1=LSU_heatmap(RbcL_US_df,'enrich_out.csv')

# Normalize colors
colors = [(0, "orchid"), (.61, "white"), (1, "green")]
custom_cmap = LinearSegmentedColormap.from_list("custom_cmap", colors, N=256)

fig, ax = plt.subplots(figsize=(90,15))
ax = sns.heatmap(L1, cmap = custom_cmap) #linewidths=0.5, linecolor='black' )
#cmap = custom_cmap

#color values with no selection data as background
ax.set_facecolor('gainsboro')
ax.grid(False)
ax.collections[0].cmap.set_bad('gainsboro')

ax.set(xlabel="Sequence", ylabel="Mutation")

#set AA x-axis range

ax.set_xlim(1,477)

#change to desired name of file

#plt.savefig('enrichment.png', dpi=300)


