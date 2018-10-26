from utils import *

def MakeGeneColumn():
    reader = csv.reader(open("/Users/jiayao/Work/BrainDisorders/data/expression/brainspan/gene_matrix/columns_metadata.csv", 'rb'))
    header = reader.next()
    idx_donor_id, idx_structure_id = header.index("donor_id"), header.index("structure_id")
    column_ids = []
    for row in reader:
        column_id = row[idx_donor_id] + "." + row[idx_structure_id]
        column_ids.append(column_id)
    return column_ids

def MakeGeneRow():
    reader = csv.reader(open("/Users/jiayao/Work/BrainDisorders/data/expression/brainspan/gene_matrix/rows_metadata.csv", 'rb'))
    header = reader.next()
    idx_ensembl_gene_id = header.index("ensembl_gene_id")
    res = []
    for row in reader:
        res.append(row[idx_ensembl_gene_id])
    return res

def LoadGeneExp():
    reader = csv.reader(open("/Users/jiayao/Work/BrainDisorders/data/expression/brainspan/gene_matrix/expression_matrix.csv", 'rb'))
    header = MakeGeneColumn()
    rower = MakeGeneRow()
    res = {}
    for i, row in enumerate(reader):
        for j, item in enumerate(row):
            if j == 0:
                continue
            #print i,j
            ID = rower[i] + "." + header[j-1]
            res[ID] = item
    return res

def MakeExonColumn():
    reader = csv.reader(open("/Users/jiayao/Work/BrainDisorders/data/expression/brainspan/exons_matrix/columns_metadata.csv", "rb"))
    header = reader.next()
    idx_donor_id, idx_structure_id = header.index("donor_id"), header.index("structure_id")
    column_ids = []
    for row in reader:
        column_id = row[idx_donor_id] + "." + row[idx_structure_id]
        column_ids.append(column_id)
    return column_ids
   
def MakeExonRow():
    reader = csv.reader(open("/Users/jiayao/Work/BrainDisorders/data/expression/brainspan/exons_matrix/rows_metadata.csv", 'rb'))
    header = reader.next()
    idx_ensembl_gene_id = header.index("ensembl_gene_id")
    res = []
    for row in reader:
        res.append(row[idx_ensembl_gene_id])
    return res

def NormExonExpByGene():
    ExomRPKM = csv.reader(open("/Users/jiayao/Work/BrainDisorders/data/expression/brainspan/exons_matrix/expression_matrix.csv"))
    NormRPKM = csv.writer(open("NormedExonExpMatrix.csv", 'wb'), delimiter=",")
    GeneValues = LoadGeneExp()
    header = MakeExonColumn()
    rower = MakeExonRow()
    for i,row in enumerate(ExomRPKM):
        new_row = [i+1]
        tmp = dict(zip(header, row))
        for j, item in enumerate(row):
            if j == 0:
                continue
            ID = rower[i] + "." + header[j-1]
            GeneExp = GeneValues[ID]
            try:
                if float(item) == 0:
                    new_row.append(0)
                else:
                    new_row.append(round(float(item)/float(GeneExp),6))
            except ZeroDivisionError:
                print i, j, item, ID
                exit()
        NormRPKM.writerow(new_row)

def makeAllGene():
    GeneExp = pd.read_csv("/Users/jiayao/Work/BrainDisorders/data/expression/brainspan/gene_matrix/expression_matrix.csv", header=None)
    GeneRow = pd.read_csv("/Users/jiayao/Work/BrainDisorders/data/expression/brainspan/gene_matrix/rows_metadata.csv")
    GeneCol = pd.read_csv("/Users/jiayao/Work/BrainDisorders/data/expression/brainspan/gene_matrix/columns_metadata.csv")
    GeneCol["Period"] = GeneCol.apply(lambda row: TemporalMap(row["age"])[0], axis=1)
    GeneCol["Description"] = GeneCol.apply(lambda row: TemporalMap(row["age"])[1], axis=1)
    GeneDat = [GeneExp, GeneRow, GeneCol]
    Regionsgt20 = ['OFC', 'VFC', 'HIP', 'ITC', 'AMY', 'DFC', 'STC', 'MFC', 'STR', 'IPC', 'V1C', 'S1C', 'A1C', 'M1C', 'CBC', 'MD']
    test_genes = GeneRow["gene_symbol"]
    #test_genes = GeneRow.head(100)["gene_symbol"]
    DisplayGeneSumExpViolin(test_genes, GeneDat, Regionsgt20)

def DisplayGeneSumExpViolin(Genes, GeneDat, structure_acronyms):
    GeneExp, GeneRow, GeneCol = GeneDat
    SelectedGenes = GeneRow[GeneRow["gene_symbol"].isin(Genes)]
    SelectGeneExp = GeneExp[GeneExp[0].isin(list(SelectedGenes["row_num"]))]
    gene_ids = list(SelectedGenes["row_num"])
    stages = {}
    for i, gene_id in enumerate(gene_ids):
        sys.stdout.write("\rLoading Genes {}".format(i))
        for j, structure_acronym in enumerate(structure_acronyms): # sum over structures/regions
            region_meta = GeneCol[GeneCol["structure_acronym"]==structure_acronym]
            for index, row in region_meta.iterrows():
                Period = row["Period"]
                if Period not in stages:
                    stages[Period] = []
                stages[Period].append(GeneExp.get_value(gene_id-1, row["column_num"]))
    dat = []
    for stage in Stages:
        #print stage
        dat.append(stages[stage])
    ax = sns.violinplot(data=dat)
    #ax = sns.boxplot(data=dat)
    ax.xaxis.set_major_formatter(plt.FuncFormatter(format_func2))
    plt.title(";".join(structure_acronyms))
    plt.ylim(0,20)
    plt.savefig("ALL.Genes.RPKM.Violin.pdf", format='pdf')


def main():
    makeAllGene()


main()
