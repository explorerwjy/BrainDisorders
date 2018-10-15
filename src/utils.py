import pandas as pd
import csv 
import re
import numpy as np
import math
import matplotlib as mpl
import matplotlib.pyplot as plt
import os
import sys
import seaborn as sns

# Files
protein_coding_gene_file = "/Users/jiayao/Work/Resources/protein-coding_gene.txt"
wigler_predicted_lgd = "/Users/jiayao/Work/BrainDisorders/data/functions/wigler-predicted-lgd.txt"
wigler_predicted_male = "/Users/jiayao/Work/BrainDisorders/data/functions/wigler-predicted-male.txt"
wigler_predicted_fem = "/Users/jiayao/Work/BrainDisorders/data/functions/wigler-predicted-fem.txt"
wigler_fam_info = "/Users/jiayao/Work/BrainDisorders/data/nature13908-s2/Supplementary_Table_1.xlsx" #IQ in it

# Entrez ID of genes in each category; No overlap
vitkup_psd = "/Users/jiayao/Work/BrainDisorders/data/functions/asd-vitkup-snv-cnv-psd.txt"
vitkup_sig_skel = "/Users/jiayao/Work/BrainDisorders/data/functions/asd-vitkup-snv-cnv-sig-skel.txt"
vitkup_channel = "/Users/jiayao/Work/BrainDisorders/data/functions/asd-vitkup-snv-cnv-channel.txt"
vitkup_chromatin = "/Users/jiayao/Work/BrainDisorders/data/functions/asd-vitkup-snv-cnv-chromatin.txt"

class H2Analysis:
    def __init__(self):
        pass
    def PlotQQ(P, N=None):
        font = {'family' : 'normal',
            'weight' : 'normal',
            'size'   : 15}
        mpl.rc('font', **font)
        P = np.array(P) 
        P.sort(reverse=True)
        Q = [-1*math.log(x,10) for x in P]
        top = int(Q[-1]) + 1
        
        N = len(P) if N == None else N
        P_null = [float(i+1)/N for i in xrange(N)]
        P_null.sort(reverse=True)
        Q_null = [-1*math.log(x,10) for x in P_null]
        colors = (0,0,0)
        plt.scatter(Q_null, Q, c=colors, alpha=0.5)
        plt.plot([0, top], [0, top], ls="-")
        plt.xlabel('Exp Q')
        plt.ylabel('Obs Q')
        plt.show()
    def _gen_data(fhs, columns, sep):
        """
        iterate over the files and yield chr, start, pvalue
        """
        for fh in fhs:
            for line in fh:
                if line[0] == "#": continue
                toks = line.split(sep)
                yield toks[columns[0]], int(toks[columns[1]]), float(toks[columns[2]])
    def chr_cmp(a, b):
        a = a.lower().replace("_", ""); b = b.lower().replace("_", "")
        achr = a[3:] if a.startswith("chr") else a
        bchr = b[3:] if b.startswith("chr") else b

        try:
            return cmp(int(achr), int(bchr))
        except ValueError:
            if achr.isdigit() and not bchr.isdigit(): return -1
            if bchr.isdigit() and not achr.isdigit(): return 1
            # X Y
            return cmp(achr, bchr)
    def chr_loc_cmp(alocs, blocs):
        return chr_cmp(alocs[0], blocs[0]) or cmp(alocs[1], blocs[1])
    def manhattan(fhs, columns, image_path, no_log, colors, sep, title, lines, ymax):
        xs = []
        ys = []
        cs = []
        colors = cycle(colors)
        xs_by_chr = {}

        last_x = 0
        data = sorted(_gen_data(fhs, columns, sep), cmp=chr_loc_cmp)

        for seqid, rlist in groupby(data, key=itemgetter(0)):
            color = colors.next()
            rlist = list(rlist)
            region_xs = [last_x + r[1] for r in rlist]
            xs.extend(region_xs)
            ys.extend([r[2] for r in rlist])
            cs.extend([color] * len(rlist))

            xs_by_chr[seqid] = (region_xs[0] + region_xs[-1]) / 2

            # keep track so that chrs don't overlap.
            last_x = xs[-1]

        xs_by_chr = [(k, xs_by_chr[k]) for k in sorted(xs_by_chr.keys(), cmp=chr_cmp)]

        xs = np.array(xs)
        ys = np.array(ys) if no_log else -np.log10(ys)

        plt.close()
        f = plt.figure()
        ax = f.add_axes((0.1, 0.09, 0.88, 0.85))

        if title is not None:
            plt.title(title)

        ax.set_ylabel('-log10(p-value)')
        if lines:
            ax.vlines(xs, 0, ys, colors=cs, alpha=0.5)
        else:
            ax.scatter(xs, ys, s=2, c=cs, alpha=0.8, edgecolors='none')

        # plot 0.05 line after multiple testing.
        ax.axhline(y=-np.log10(0.05 / len(data)), color='0.5', linewidth=2)
        plt.axis('tight')
        plt.xlim(0, xs[-1])
        plt.ylim(ymin=0)
        if ymax is not None: plt.ylim(ymax=ymax)
        plt.xticks([c[1] for c in xs_by_chr], [c[0] for c in xs_by_chr], rotation=-90, size=8.5)
        print >>sys.stderr, "saving to: %s" % image_path
        plt.savefig(image_path)
        #plt.show() 

def get_gene_entrez_symbol_map():
    res = {}
    reader = csv.reader(open(protein_coding_gene_file), delimiter="\t")
    header = reader.next()
    idx_symbol, idx_entrez_id = header.index("symbol"), header.index("entrez_id")
    for row in reader:
        res[row[idx_entrez_id]] = row[idx_symbol]
    return res
class DevStageExp:
    def __init__(self, DevStage):
        self.DevStage = DevStage
        self.exps = np.array([])
        self.mean_exp = 0
    def add_exp(self, value):
        self.exps = np.append(self.exps, float(value))
        #self.mean_exp = round(np.mean(self.exps), 4)
    def get_nonzero_mean(self):
        try:
            self.non_zero_exps = [x for x in self.exps if x]
            return sum(self.non_zero_exps)/len(self.non_zero_exps)
        except ZeroDivisionError:
            return 0
class BrainSpan:
    def __init__(self):
        # store all exp value at a developemt time point
        #Stages = map(str,range(1,16))
        #self.Stages = ["1", "2A", "2B", "3A", "3B"] + map(str,range(4,12))
        self.Stages = ["2A", "2B", "3A", "3B"] + map(str,range(4,12))
    def TemporalMap(self, age):
        num, unit = age.split()
        num = float(num)
        if unit == 'pcw':
            if num >= 4 and num < 8:
                return ("1","1", "Embryonic")
            elif num >= 8 and num <10:
                return ("2A","2", "Early prenatal")
            elif num >= 10 and num <13:
                return ("2B", "3", "Early prenatal")
            elif num >= 13 and num <16:
                return ("3A", "4", "Early mid-prenatal")
            elif num >= 16 and num <19:
                return ("3B", "5", "Early mid-prenatal")
            elif num >= 19 and num <24:
                return ("4", "6", "Late mid-prenatal")
            elif num >= 24 and num <38:
                return ("5", "7", "Late prenatal")
            else:
                print "Unexpected Value"
        elif unit == 'mos':
            if num>= 0 and num < 6:
                return ("6", "8", "Early infancy")
            elif num>= 6 and num < 12:
                return ("7", "9", "Late infancy")
        elif unit == 'yrs':
            if num >= 1 and num < 6:
                return ("8", "10", "Early childhood")
            elif num >= 6 and num < 12:
                return ("9", "11", "Late childhood")
            elif num >= 12 and num < 20:
                return ("10", "12", "Adolescence")
            elif num >= 20 and num < 40:
                return ("11", "13", "Young adulthood")
            elif num >= 40 and num < 60:
                return ("11", "14", "Young adulthood")
            elif num >= 60:
                return ("11", "15", "Young adulthood")
    def Insert(self, row, Dict):
        return Dict[row["Gene"]]
    # Return res1: sequence of exp in 12 dev stages res2: num of samples at each timepoint
    def DisplayStageExp(self, Dict, out="nonzero"):
        res1 = []
        res2 = []
        for stage in self.Stages:
            try:
                #res1.append(Dict[stage].mean_exp) # Exp avg
                stage, Dict[stage].get_nonzero_mean()
                if out == "all":
                    res1.append(Dict[stage].mean_exp)
                else:
                    res1.append(Dict[stage].get_nonzero_mean()) # Exp avg
                res2.append(len(Dict[stage].exps)) # Num of samples at taht time point
            except KeyError:
                res1.append(0)
                res2.append(0)
        return res1, res2
    def Look(self, Gene, structure_acronym, bp_exon_row_meta, bp_exon_col_meta, ExonExp, smooth=True, drop_low_exp=True, fontsize=6):
        Exons = bp_exon_row_meta[bp_exon_row_meta["gene_symbol"]==Gene]
        GeneExonExp = ExonExp[ExonExp[0].isin(list(Exons["row_num"]))]
        ID_start, ID_end = list(Exons["row_num"])[0], list(Exons["row_num"])[-1]
        region_exon_meta = bp_exon_col_meta[bp_exon_col_meta["structure_acronym"]==structure_acronym]
        plt.close('all')
        font = {'family' : 'normal',
            'weight' : 'normal',
            'size'   : fontsize}
        mpl.rc('font', **font)
        res = {}
        for exon_id in xrange(ID_start,ID_end+1):
            stages = {}
            for index,row in region_exon_meta.iterrows():
                Period = row["Period"]
                if Period not in stages:
                    stages[Period] = DevStageExp(Period)
                stages[Period].add_exp(GeneExonExp.get_value(exon_id-1, row["column_num"]))
            seq, Nsamples = DisplayStageExp(stages)
            if drop_low_exp:
                if seq.count(0) >= 4:
                    res[exon_id] = None
                    continue
            if smooth:
                res[exon_id] = self.smooth_func(seq)
            else:
                res[exon_id] = seq 
        fig, ax = plt.subplots(dpi=200)
        plt.title("Gene:{} Region:{}".format(Gene, structure_acronym))
        #largest = 0
        for exon_id in xrange(ID_start,ID_end+1):
            if res[exon_id] == None:
                continue
            #if min(np.array([math.log(x, 2) if x!=0 else 0 for x in res[exon_id]])) < -2:
            #if (np.array([math.log(x, 2) if x!=0 else 0 for x in res[exon_id]])[7]) > largest:
            #    largest = (np.array([math.log(x, 2) if x!=0 else 0 for x in res[exon_id]])[7])
            #    largest_id = exon_id
            ax.plot(range(2,14),[math.log(x, 2) if x!=0 else 0 for x in res[exon_id]])
        #print exon_id
        ax.grid(True)
        ax.axvline(x=7.5)
        #plt.ylim(-5, 2)
        #Stages = Stages + [0,0,0,0,0]
        plt.xticks(np.arange(2,14), Stages, rotation=20)
        ax.xaxis.set_major_formatter(plt.FuncFormatter(self.format_func))
        plt.show()
    # I forget what it does
    def Look2(self, Gene, structure_acronym, bp_exon_row_meta, bp_exon_col_meta, ExonExp, fontsize=6):
        Exons = bp_exon_row_meta[bp_exon_row_meta["gene_symbol"]==Gene]
        GeneExonExp = ExonExp[ExonExp[0].isin(list(Exons["row_num"]))]
        ID_start, ID_end = list(Exons["row_num"])[0], list(Exons["row_num"])[-1]
        region_exon_meta = bp_exon_col_meta[bp_exon_col_meta["structure_acronym"]==structure_acronym]
        plt.close('all')
        font = {'family' : 'normal',
            'weight' : 'normal',
            'size'   : fontsize}
        mpl.rc('font', **font)
        res = {}
        for exon_id in xrange(ID_start,ID_end+1):
            stages = {}
            for index,row in region_exon_meta.iterrows():
                Period = row["Period"]
                if Period not in stages:
                    stages[Period] = DevStageExp(Period)
                stages[Period].add_exp(GeneExonExp.get_value(exon_id-1, row["column_num"]))
            res[exon_id] = DisplayStageExp(stages)
        fig, ax = plt.subplots(dpi=200)
        plt.title("Gene:{} Region:{}".format(Gene, structure_acronym))
        for exon_id in xrange(ID_start,ID_end+1):
            #if np.array([math.log(x, 2) if x!=0 else 0 for x in res[exon_id]]).
            #ax.plot(range(2,14),[math.log(x, 2) if x!=0 else 0 for x in res[exon_id]])
            ax.plot(range(2,14),res[exon_id])
        ax.grid(True)
        ax.axvline(x=7.6)
        #Stages = Stages + [0,0,0,0,0]
        plt.xticks(np.arange(2,14), Stages, rotation=20)
        ax.xaxis.set_major_formatter(plt.FuncFormatter(self.format_func))
        plt.show()
    def example_plot(self, ax, structure_acronym, fontsize=6):
        ax.locator_params(nbins=3)
        #ax.set_xlabel('Dev Stages', fontsize=fontsize)
        ax.set_xlabel('Dev Stages', fontsize=fontsize)
        ax.set_ylabel('log2 fold change', fontsize=fontsize)
        ax.set_title(structure_acronym, fontsize=fontsize)
    # Display exon express in a series of regions, for a certain gene
    def LookGrid(self, Gene, structure_acronyms, bp_exon_row_meta, bp_exon_col_meta, ExonExp, smooth=True, drop_low_exp=True, fontsize=6):
        Exons = bp_exon_row_meta[bp_exon_row_meta["gene_symbol"]==Gene]
        GeneExonExp = ExonExp[ExonExp[0].isin(list(Exons["row_num"]))]
        ID_start, ID_end = list(Exons["row_num"])[0], list(Exons["row_num"])[-1]
        NN = int(math.sqrt(len(structure_acronyms))) 
        #print len(structure_acronyms), NN
        if NN * NN == len(structure_acronyms):
            XX, YY = NN, NN
        else:
            XX, YY = NN+1 , NN + 1
        plt.close('all')
        fig = plt.figure(dpi=800)
        font = {'family' : 'normal',
            'weight' : 'normal',
            'size'   : fontsize}
        mpl.rc('font', **font)
        from mpl_toolkits.axes_grid1 import Grid
        grid = Grid(fig, rect=111, nrows_ncols=(XX,YY), axes_pad=0.25, label_mode='L',)
        #print grid
        for i, structure_acronym in enumerate(structure_acronyms):
            sys.stdout.write("\r{}".format(i+1))
            region_exon_meta = bp_exon_col_meta[bp_exon_col_meta["structure_acronym"]==structure_acronym]
            res = {}
            for exon_id in xrange(ID_start,ID_end+1):
                stages = {}
                for index,row in region_exon_meta.iterrows():
                    Period = row["Period"]
                    if Period not in stages:
                        stages[Period] = DevStageExp(Period)
                    stages[Period].add_exp(GeneExonExp.get_value(exon_id-1, row["column_num"]))
                seq, Nsamples = DisplayStageExp(stages)
                if drop_low_exp:
                    if seq.count(0) >= 4:
                        res[exon_id] = None
                        continue
                res[exon_id] = smooth_func(seq) if smooth else seq
                #if smooth:
                #    #print DisplayStageExp(stages)[0]
                #    res[exon_id] = smooth_func(seq)
                #    #print res[exon_id]
                #else:
                #    res[exon_id] = seq
            ax = grid[i]
            example_plot(ax, structure_acronym)
            for exon_id in xrange(ID_start,ID_end+1):
                if res[exon_id] == None:
                    continue
                #ax.plot(range(2,14),res[exon_id])
                #ax.plot(range(2,14),[math.log(x, 2) if x!=0 else 0 for x in res[exon_id]])
                ax.plot(range(2,14),[math.log(x, 2) if x!=0 else 0 for x in res[exon_id]])
            ax.grid(True)
            ax.axvline(x=7.5)
            #plt.ylim(-2, 2)
            #plt.xticks(np.arange(2,14,2), Stages, rotation=20)
            #ax.xaxis.set_major_formatter(plt.FuncFormatter(self.format_func))
        plt.tight_layout(pad=4, w_pad=5, h_pad=10)
        #plt.xticks(np.arange(2,14,2), Stages, rotation=20)
        #plt.ylim(-2, 2)
        sys.stdout.write("\rplotting....")
        plt.show()
    def smooth_func(self, seq):
        seq.insert(0,0)
        seq.append(0)
        new = []
        i = 1
        while i < len(seq) - 1:
            nonzero = [ x for x in seq[i-1: i+2] if x ]
            try:
                new.append(sum(nonzero)/len(nonzero))
            except ZeroDivisionError:
                print "Three 0 data at stage",i
                new.append(0)
            i += 1
        return new
    def format_func(self, value, tick_number):
        try:
            return Stages[int(value)-2]
        except:
            #print value
            return 0
    class expvalue:
        def __init__(self):
            self.value = []
        def addvalue(self, value):
            self.value.append(value)
            self.mean = self.value.mean()
    def LookGridSumRegion(self, Gene, structure_acronyms, bp_exon_row_meta, bp_exon_col_meta, ExonExp, smooth=True, drop_low_exp=True, fontsize=6):
        Exons = bp_exon_row_meta[bp_exon_row_meta["gene_symbol"]==Gene]
        GeneExonExp = ExonExp[ExonExp[0].isin(list(Exons["row_num"]))]
        ID_start, ID_end = list(Exons["row_num"])[0], list(Exons["row_num"])[-1]
        plt.close('all')
        fig = plt.figure(dpi=800)
        font = {'family' : 'normal',
            'weight' : 'normal',
            'size'   : fontsize}
        mpl.rc('font', **font)
        from mpl_toolkits.axes_grid1 import Grid
        #grid = Grid(fig, rect=111, nrows_ncols=(XX,YY), axes_pad=0.25, label_mode='L',)
        res = {}
        for exon_id in xrange(ID_start, ID_end+1):
            stages = {}
            for j, structure_acronym in enumerate(structure_acronyms):
                region_exon_meta = bp_exon_col_meta[bp_exon_col_meta["structure_acronym"]==structure_acronym]
                for index, row in region_exon_meta.iterrows():
                    Period = row["Period"]
                    if Period not in stages:
                        stages[Period] = DevStageExp(Period)
                    stages[Period].add_exp(GeneExonExp.get_value(exon_id-1, row["column_num"]))
            #res[exon_id]
            seq = DisplayStageExp(stages)[0]
            if drop_low_exp:
                if seq.count(0) >= 4:
                    res[exon_id] = None
                    continue
            if smooth:
                res[exon_id] = smooth_func(seq)
            else:
                res[exon_id] = seq

        fig, ax = plt.subplots(dpi=200)
        plt.title("Gene:{} over Region:{}".format(Gene, ",".join(structure_acronyms)))
        for exon_id in xrange(ID_start,ID_end+1):
            if res[exon_id] == None:
                continue
            #if np.array([math.log(x, 2) if x!=0 else 0 for x in res[exon_id]]).
            ax.plot(range(2,14),[math.log(x, 2) if x!=0 else 0 for x in res[exon_id]])
            #ax.plot(range(2,14),res[exon_id])
            #print len(res[exon_id]), res[exon_id]
            ax.plot(range(2,14),res[exon_id])
        ax.grid(True)
        ax.axvline(x=7.5)
        #Stages = Stages + [0,0,0,0,0]
        plt.xticks(np.arange(2,14), Stages, rotation=20)
        ax.xaxis.set_major_formatter(plt.FuncFormatter(format_func))
        plt.show()
    def NoiseEstimate():
        pass
    def AssignVar2Exon(self, bp_exon_row_meta, VarFile):
        reader = csv.reader(open(VarFile, 'rb'))
        header = reader.next()
        res = {} # k:gene v:list of variants in that gene
        for row in reader:
            tmp = dict(zip(header, row))
            if int(tmp["GeneCount"]) < 2:
                continue
            if tmp["Gene"] not in res:
                res[tmp["Gene"]] = [tmp]
            else:
                res[tmp["Gene"]].append(tmp)
        bp_exon_row_meta["Vars"] = ""
        for i, row in bp_exon_row_meta.iterrows():
            sys.stdout.write("\r{}".format(i))
            start, end = int(row["start"]), int(row["end"])
            if row["gene_symbol"] not in res:
                continue
            for var in res[row["gene_symbol"]]:
                pos = int(var["Position"])
                if pos >= start and pos <= end:
                    ##bp_exon_row_meta["Vars"].append(var["cDnaVariant"])
                    row["Vars"] = row["Vars"] + ";" + var["cDnaVariant"] 
                    bp_exon_row_meta.at[i, "Vars"] = bp_exon_row_meta.get_value(i, "Vars") + ";" + var["cDnaVariant"] 
                    continue
    # Assign Category to a gene
    # Cates are dict of category:geneset
    def AssignFunc(self, gene, Cates): # Assign Category to each gene
        res = []
        for k,v in Cates.items():
            if gene in v:
                res.append(k)
        return ";".join(res)
    #def GetIQ(gene, FamID):
    # Assign Variants to row meta data file of Exons. 
    # Markers on 1)Male vs. Female; 2)IQ70+ vs. IQ70- 3)Functions
    def AssignVar2Exon2(self, bp_exon_row_meta, VarFile):
        VarDF = pd.read_excel(VarFile)
        # Load annotations to gene
        entrez_symbol_map = get_gene_entrez_symbol_map()
        wigler_predicted_lgd_genes = set([entrez_symbol_map[x.strip()] for x in file(wigler_predicted_lgd)])
        vitkup_channel_genes = set([entrez_symbol_map[x.strip()] for x in file(vitkup_channel)])
        vitkup_chromatin_genes = set([entrez_symbol_map[x.strip()] for x in file(vitkup_chromatin)])
        vitkup_psd_genes = set([entrez_symbol_map[x.strip()] for x in file(vitkup_psd)])
        vitkup_sig_skel_genes = set([entrez_symbol_map[x.strip()] for x in file(vitkup_sig_skel)])
        fam_info = pd.read_excel(wigler_fam_info)
        famID2NVIQ = dict(zip(list(fam_info["familyId"]), list(fam_info["probandNVIQ"])))
        famID2VIQ = dict(zip(list(fam_info["familyId"]), list(fam_info["probandVIQ"])))
        
        Cates = dict(zip(["channel", "chromatin", "psd", "sig_skel"],[vitkup_channel_genes,vitkup_chromatin_genes,vitkup_psd_genes,vitkup_sig_skel_genes]))
        res = {} # k:gene v:list of variants in that gene
        var_genders = {}
        var_IQ = {}
        for i, row in VarDF.iterrows():
            if row["effectGene"] not in wigler_predicted_lgd_genes:
                continue
            # Get Variant Info
            gene,var = row["effectGene"], row["location"]
            if gene not in res:
                res[gene] = [row["location"]]
            else:
                res[gene].append(row["location"])
            # Get Gender Info
            if var not in var_genders:
                var_genders[var] = row["inChild"][:2]
            else:
                var_genders[var].append(row["inChild"][:2])
            # Get IQ Info
            if var not in var_IQ:
                FamID = row["familyId"]
                var_IQ[var] = (famID2VIQ[FamID],famID2NVIQ[FamID])
        bp_exon_row_meta["NVIQ70"] = ""
        bp_exon_row_meta["VIQ70"] = ""
        bp_exon_row_meta["Gender"] = ""
        bp_exon_row_meta["Func"] = ""
        bp_exon_row_meta["Vars"] = ""
        bp_exon_row_meta["GeneHited"] = "F"
        bp_exon_row_meta["Last"] = "F"
        LAST_GENE = None
        for i, row in bp_exon_row_meta.iterrows():
            sys.stdout.write("\r{}".format(i))
            gene = row["gene_symbol"]
            if gene != LAST_GENE and gene != None:
                bp_exon_row_meta.at[i-1, "Last"] = "T" 
            LAST_GENE = gene
            if gene not in res:
                continue
            bp_exon_row_meta.at[i, "GeneHited"] = "T"
            start, end = int(row["start"]), int(row["end"])
            gene = row["gene_symbol"]
            # Assign Variant to Exon
            for var in res[gene]:
                pos = int(var.split(":")[1])
                #print pos
                if pos >= start and pos <= end: #Var In This Exon
                    bp_exon_row_meta.at[i, "Vars"] = bp_exon_row_meta.get_value(i, "Vars") + ";" + var if bp_exon_row_meta.get_value(i, "Vars") != "" else var
                    bp_exon_row_meta.at[i, "VIQ70"] = var_IQ[var][0]
                    bp_exon_row_meta.at[i, "NVIQ70"] = var_IQ[var][1]
                    bp_exon_row_meta.at[i, "Gender"] = bp_exon_row_meta.get_value(i, "Gender") + ";" + var_genders[var] if bp_exon_row_meta.get_value(i, "Gender") != "" else var_genders[var]
                    continue
            # Assign Other Info
            bp_exon_row_meta.at[i, "Func"] = self.AssignFunc(gene, Cates)
        bp_exon_row_meta.at[i, "Last"] = "T" #Last Exon in DF is last exon of that gene
    def LookMutationTargetedExon(self, Gene, structure_acronyms, bp_exon_row_meta, bp_exon_col_meta, ExonExp, GeneDat=None, smooth=True, drop_low_exp=True, fontsize=6):
        Exons = bp_exon_row_meta[bp_exon_row_meta["gene_symbol"]==Gene]
        GeneExonExp = ExonExp[ExonExp[0].isin(list(Exons["row_num"]))]
        exon_ids = list(Exons["row_num"])
        plt.close('all')
        fig = plt.figure(dpi=800)
        font = {'family' : 'normal',
            'weight' : 'normal',
            'size'   : fontsize}
        mpl.rc('font', **font)
        from mpl_toolkits.axes_grid1 import Grid
        res = {}
        for exon_id in exon_ids: # Plot each exon
            stages = {}
            for j, structure_acronym in enumerate(structure_acronyms): # sum over structures/regions
                region_exon_meta = bp_exon_col_meta[bp_exon_col_meta["structure_acronym"]==structure_acronym]
                for index, row in region_exon_meta.iterrows(): 
                    Period = row["Period"]
                    if Period not in stages:
                        stages[Period] = DevStageExp(Period)
                    stages[Period].add_exp(GeneExonExp.get_value(exon_id-1, row["column_num"]))
            seq = DisplayStageExp(stages)[0]
            if drop_low_exp:
                if seq.count(0) >= 4:
                    res[exon_id] = None
                    continue
            if smooth:
                res[exon_id] = smooth_func(seq)
            else:
                res[exon_id] = seq
        fig, ax = plt.subplots(dpi=200)
        plt.title("Gene:{} over Region:{}".format(Gene, ",".join(structure_acronyms)))
        for exon_id in exon_ids: 
            if res[exon_id] == None:
                continue
            ax.plot(range(2,14),[math.log(x, 2) if x!=0 else 0 for x in res[exon_id]], label=exon_id)
        if GeneDat != None:# Plot Gene exp over stages 
            GeneExp, GeneRow, GeneCol = GeneDat
            Genes = GeneRow[GeneRow["gene_symbol"]==Gene]
            Gene = GeneExp[GeneExp[0].isin(list(Genes["row_num"]))]
            gene_ids = list(Genes["row_num"])
            for gene_id in gene_ids:
                for j, structure_acronym in enumerate(structure_acronyms): # sum over structures/regions
                    region_meta = GeneCol[GeneCol["structure_acronym"]==structure_acronym]
                    for index, row in region_meta.iterrows():
                        Period = row["Period"]
                        if Period not in stages:
                            stages[Period] = DevStageExp(Period)
                        stages[Period].add_exp(GeneExp.get_value(gene_id-1, row["column_num"]))
                seq, Nsample = DisplayStageExp(stages)
                print seq
                if smooth:
                    res[gene_id] = smooth_func(seq)
                else:
                    res[gene_id] = seq
            #add_layout(LinearAxis(y_range_name="GeneRPKM"), 'right')
            ax2 = ax.twinx()
            for gene_id in gene_ids:
                print res[gene_id]
                ax2.plot(range(2,14), [math.log(x,2) if x !=0 else 0 for x in res[gene_id]], color="black") 
                #ax2.plot(range(2,14), res[gene_id], color="black") 
            ax2.set_ylabel("GeneRPKM")
        ax.legend()
        ax.grid(True)
        ax.axvline(x=7.5)
        #Stages = Stages + [0,0,0,0,0]
        plt.xticks(np.arange(2,14), Stages, rotation=20)
        ax.xaxis.set_major_formatter(plt.FuncFormatter(format_func))
        plt.show()
    def LookMutationTargetedGene(self, Gene, structure_acronyms, bp_exon_row_meta, bp_exon_col_meta, ExonExp, smooth=True, drop_low_exp=True, fontsize=6):
        Exons = bp_exon_row_meta[bp_exon_row_meta["gene_symbol"]==Gene]
        #print Exons
        GeneExonExp = ExonExp[ExonExp[0].isin(list(Exons["row_num"]))]
        #ID_start, ID_end = list(Exons["row_num"])[0], list(Exons["row_num"])[-1]
        exon_ids = list(Exons["row_num"])
        #print exon_ids
        plt.close('all')
        fig = plt.figure(dpi=800)
        font = {'family' : 'normal',
            'weight' : 'normal',
            'size'   : fontsize}
        mpl.rc('font', **font)
        from mpl_toolkits.axes_grid1 import Grid
        #grid = Grid(fig, rect=111, nrows_ncols=(XX,YY), axes_pad=0.25, label_mode='L',)
        res = {}
        for exon_id in exon_ids:
            stages = {}
            for j, structure_acronym in enumerate(structure_acronyms):
                region_exon_meta = bp_exon_col_meta[bp_exon_col_meta["structure_acronym"]==structure_acronym]
                for index, row in region_exon_meta.iterrows():
                    Period = row["Period"]
                    if Period not in stages:
                        stages[Period] = DevStageExp(Period)
                    stages[Period].add_exp(GeneExonExp.get_value(exon_id-1, row["column_num"]))
            #res[exon_id]
            seq = DisplayStageExp(stages)[0]
            if drop_low_exp:
                if seq.count(0) >= 4:
                    res[exon_id] = None
                    continue
            if smooth:
                res[exon_id] = smooth_func(seq)
            else:
                res[exon_id] = seq

        fig, ax = plt.subplots(dpi=200)
        plt.title("Gene:{} over Region:{}".format(Gene, ",".join(structure_acronyms)))
        for exon_id in exon_ids: 
            if res[exon_id] == None:
                continue
            #if np.array([math.log(x, 2) if x!=0 else 0 for x in res[exon_id]]).
            ax.plot(range(2,14),[math.log(x, 2) if x!=0 else 0 for x in res[exon_id]])
            #ax.plot(range(2,14),res[exon_id])
            #print len(res[exon_id]), res[exon_id]
            #ax.plot(range(2,14),res[exon_id])
        ax.grid(True)
        ax.axvline(x=7.5)
        #Stages = Stages + [0,0,0,0,0]
        plt.xticks(np.arange(2,14), Stages, rotation=20)
        ax.xaxis.set_major_formatter(plt.FuncFormatter(format_func))
        plt.show()
    def LoadingDat2SeqCrossRecordCrossRegion(self, record_ids, regions, row_meta, col_meta, Matrix, smooth=True, drop_low_exp=True):
        lengths = []
        stages = {}
        for i, _id in enumerate(record_ids):
            sys.stdout.write("\rLoading records:{}".format(i))
            lengths.append(row_meta.get_value(int(_id)-1, "exon length"))
            for j, region in enumerate(regions):
                region_col_meta = col_meta[col_meta["structure_acronym"]==region]
                for index, row in region_col_meta.iterrows():
                    Period = row["Period"]
                    if Period not in stages:
                        stages[Period] = DevStageExp(Period)
                    stages[Period].add_exp(Matrix.get_value(_id-1, row["column_num"]))
        sys.stdout.write("\n")
        seq = self.DisplayStageExp(stages)[0]
        if smooth:
            seq = self.smooth_func(seq)
        return seq, lengths
    def LookALLMutationTargetedExon2(self, exon_ids, structure_acronyms, bp_exon_row_meta, bp_exon_col_meta, ExonExp, smooth=True, drop_low_exp=True, fontsize=6, rd=False):
        plt.close('all')
        fig = plt.figure(dpi=800)
        stages = {}
        ExonLengths = []
        seq, lengths = LoadingDat2SeqCrossRecordCrossRegion(exon_ids, structure_acronyms, bp_exon_row_meta, bp_exon_col_meta, ExonExp, smooth=True) # Loading data needed
        fig, ax = plt.subplots(dpi=200)
        plt.title("TargetedExons over Region:{}".format(",".join(structure_acronyms)))
        ax.plot(range(2,14),[math.log(x, 2) if x!=0 else 0 for x in seq], label="TargetedExons")
        ax.grid(True)
        ax.axvline(x=7.5)
        if rd != False:
            if rd == "gt1000":
                rd_exon_ids = list(bp_exon_row_meta[bp_exon_row_meta["exon length"]>=1000].sample(len(exon_ids)*10)["row_num"])
            if rd == "lt1000":
                rd_exon_ids = list(bp_exon_row_meta[bp_exon_row_meta["exon length"]<1000].sample(len(exon_ids)*10)["row_num"])
            rd_seq, rd_length = LoadingDat2SeqCrossRecordCrossRegion(rd_exon_ids,structure_acronyms, bp_exon_row_meta, bp_exon_col_meta, ExonExp, smooth=True)
            ax.plot(range(2,14),[math.log(x, 2) if x!=0 else 0 for x in rd_seq], label="RandomExons")
        plt.xticks(np.arange(2,14), Stages, rotation=20)
        ax.xaxis.set_major_formatter(plt.FuncFormatter(self.format_func))
        ax.legend()
        plt.show()
        return ExonLengths
    def LookALLMutationTargetedExon(self, exon_ids_sets, structure_acronyms, bp_exon_row_meta, bp_exon_col_meta, ExonExp, smooth=True, drop_low_exp=True, fontsize=6, rd=False):
        plt.close('all')
        fig = plt.figure(dpi=800)
        stages = {}
        fig, ax = plt.subplots(dpi=200)
        plt.title("Exons Exp over Region:{}".format(",".join(structure_acronyms)))
        for SetName, (color, ExonIDs) in exon_ids_sets.items():
            print SetName
            ExonIDs = [int(x) for x in ExonIDs]
            seq, lengths = self.LoadingDat2SeqCrossRecordCrossRegion(ExonIDs, structure_acronyms, bp_exon_row_meta, bp_exon_col_meta, ExonExp, smooth=True) # Loading data needed
            if "Other" in SetName:
                ax.plot(range(2,14),seq, '--', label=SetName, color=color)
                #ax.plot(range(2,14),[math.log(x, 2) if x!=0 else 0 for x in seq], '--', label=SetName, color=color)
            else:
                ax.plot(range(2,14),seq, label=SetName, color=color)
                #ax.plot(range(2,14),[math.log(x, 2) if x!=0 else 0 for x in seq], label=SetName, color=color)
        ax.grid(True)
        ax.axvline(x=7.5)
        plt.xticks(np.arange(2,14), self.Stages, rotation=20)
        ax.xaxis.set_major_formatter(plt.FuncFormatter(self.format_func))
        ax.legend()
        plt.show()
    def LookALLMutationTargetedGenes(self, Genes, structure_acronyms, GeneDat, smooth=True, drop_low_exp=True, fontsize=6, ylim=None):
        GeneExp, GeneRow, GeneCol = GeneDat
        SelectedGenes = GeneRow[GeneRow["gene_symbol"].isin(Genes)]
        #GeneExonExp = ExonExp[ExonExp[0].isin(list(Exons["row_num"]))]
        gene_ids = list(SelectedGenes["row_num"])
        plt.close('all')
        fig = plt.figure(dpi=800)
        font = {'family' : 'normal',
            'weight' : 'normal',
            'size'   : fontsize}
        mpl.rc('font', **font)
        stages = {}
        for gene_id in gene_ids:
            for structure_acronym in structure_acronyms:
                #region_exon_meta = bp_exon_col_meta[bp_exon_col_meta["structure_acronym"]==structure_acronym]
                region_meta = GeneCol[GeneCol["structure_acronym"]==structure_acronym]
                for index, row in region_meta.iterrows():
                    #Period = row["Period"]
                    Period = row["Stage"]
                    if Period not in stages:
                        #stages[Period] = DevStageExp(Period)
                        stages[Period] = []
                    #stages[Period].add_exp(ExonExp.get_value(exon_id-1, row["column_num"]))
                    stages[Period].append(GeneExp.get_value(gene_id-1, row["column_num"]))
            #res[exon_id]
        #seq = DisplayStageExp(stages, out="all")[0]
        seq = []
        for stage in Stages:
            try:
                seq.append(np.mean(stages[stage]))
            except:
                print stage
                seq.append(0)
                #seq.append(0)
        if smooth:
            seq = smooth_func(seq)
        fig, ax = plt.subplots(dpi=200)
        plt.title("TargetedExons over Region:{}".format(",".join(structure_acronyms)))
        #ax.plot(range(2,14),[math.log(x, 2) if x!=0 else 0 for x in seq])
        #ax.plot(range(2,14),seq)
        #ax.plot(range(1,16),seq)
        ax.plot(range(1,16),[math.log(x, 2) if x!=0 else 0 for x in seq])
        ax.grid(True)
        ax.axvline(x=7.5)
        #Stages = Stages + [0,0,0,0,0]
        #plt.xticks(np.arange(2,14), Stages, rotation=20)
        #ax.xaxis.set_major_formatter(plt.FuncFormatter(format_func))
        if ylim != None:
            plt.ylim(ylim)
        plt.show()
    def AssignVar2Gene(self, bp_gene_row_meta, VarFile):
        reader = csv.reader(open(VarFile, 'rb'))
        header = reader.next()
        res = {} # k:gene v:list of variants in that gene
        for row in reader:
            tmp = dict(zip(header, row))
            if int(tmp["GeneCount"]) < 2:
                continue
            if tmp["Gene"] not in res:
                res[tmp["Gene"]] = [tmp]
            else:
                res[tmp["Gene"]].append(tmp)
        bp_gene_row_meta["Vars"] = ""
        for i, row in bp_gene_row_meta.iterrows():
            sys.stdout.write("\r{}".format(i))
            if row["gene_symbol"] not in res:
                continue
            else:
                for var in res[row["gene_symbol"]]:
                    bp_gene_row_meta.at[i, "Vars"] = bp_gene_row_meta.get_value(i, "Vars") + ";" + var["cDnaVariant"]
    def AssignVar2Gene2(self,bp_gene_row_meta, VarFile):
        df = pd.read_excel(VarFile, sheetname="Proband SNVs")
        

        res = {} # k:gene v:list of variants in that gene
        for row in reader:
            tmp = dict(zip(header, row))
            if tmp["Gene"] not in res:
                res[tmp["Gene"]] = [tmp]
            else:
                res[tmp["Gene"]].append(tmp)
        bp_gene_row_meta["Vars"] = ""
        for i, row in bp_gene_row_meta.iterrows():
            sys.stdout.write("\r{}".format(i))
            if row["gene_symbol"] not in res:
                continue
            else:
                for var in res[row["gene_symbol"]]:
                    bp_gene_row_meta.at[i, "Vars"] = bp_gene_row_meta.get_value(i, "Vars") + ";" + var["cDnaVariant"]
    def format_func2(self, value, tick_number):
        try:
            return Stages[int(value)]
        except:
            #print value
            return 0
    def DisplayGeneExpViolin(self, Gene, GeneDat, structure_acronyms):
        GeneExp, GeneRow, GeneCol = GeneDat
        SelectedGenes = GeneRow[GeneRow["gene_symbol"]==Gene]
        #SelectGeneExp = GeneExp[GeneExp[0].isin(list(SelectedGenes["row_num"]))]
        gene_ids = list(SelectedGenes["row_num"])
        res = {}
        for gene_id in gene_ids:
            res[gene_id] = {}
            stages = {}
            for j, structure_acronym in enumerate(structure_acronyms): # sum over structures/regions
                region_meta = GeneCol[GeneCol["structure_acronym"]==structure_acronym]
                for index, row in region_meta.iterrows():
                    Period = row["Period"]
                    if Period not in stages:
                        stages[Period] = []
                    stages[Period].append(GeneExp.get_value(gene_id-1, row["column_num"]))
            res[gene_id] = stages
        for gene_id in gene_ids:
            Xs, dat = [], []
            for stage in Stages:
                Xs.append(stage)
                dat.append(res[gene_id][stage])
            #ax = sns.boxplot(x=Xs, data=dat)
            ax = sns.violinplot(data=dat)
            ax.xaxis.set_major_formatter(plt.FuncFormatter(self.format_func2))
            plt.show()
    def DisplayGeneSumExpViolin(self, Genes, GeneDat, structure_acronyms, ylim=None):
        GeneExp, GeneRow, GeneCol = GeneDat
        SelectedGenes = GeneRow[GeneRow["gene_symbol"].isin(Genes)]
        #SelectGeneExp = GeneExp[GeneExp[0].isin(list(SelectedGenes["row_num"]))]
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
            try: 
                #dat.append([ math.log(2,x) for x in stages[stage]])
                dat.append(stages[stage])
            except:
                dat.append([0])
            #dat.append(np.mean(stages[stage]))
        plt.close('all')
        ax = sns.violinplot(data=dat)
        #ax = sns.boxplot(data=dat)
        #ax.xaxis.set_major_formatter(plt.FuncFormatter(format_func2))
        #plt.plot(xrange(12), dat)
        if ylim != None:
            plt.ylim(ylim)
        plt.show()

class HBT_DATA:
    def __init__(self):
        self.Stages = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15] 
    def extract_matrix_table(self, File, OutName):
        hand = open(File, 'rb')
        out = open(OutName, 'wb')
        flag_out = False
        for l in hand:
            if flag_out:
                out.write(l)
            if l.strip() == "!series_matrix_table_begin":
                flag_out = True
            if l.strip() == "!series_matrix_table_end":
                flag_out = False
                break
        hand.close()
        out.close()
    def get_char_ch1(self, string):
        return string.strip('"').split()[0]
    def extract_col_info(self, File, OutName):
        hand = open(File, 'rb')
        csv_writer = csv.writer(open(OutName, 'wb'), delimiter=",")
        for l in hand:
            llist = l.strip().split("\t")
            label, data = llist[0], llist[1:]
            if label == "!Sample_geo_accession":
                sample_id = [x.strip('"') for x in data]
            if label == "!Sample_source_name_ch1":
                source_name = [x.strip('"').split("_")[0] for x in data]
            if label == "!Sample_characteristics_ch1":
                if self.get_char_ch1(data[0]) == "region:":
                    region = [x.strip('"').split()[1] for x in data]
                elif self.get_char_ch1(data[0]) == "Sex:":
                    sex = [x.strip('"').split()[1] for x in data]
                elif self.get_char_ch1(data[0]) == "age:":
                    age = [x.strip('"').split()[1] for x in data]
                elif self.get_char_ch1(data[0]) == "Stage:":
                    stage = [x.strip('"').split()[1] for x in data]
            if label == "!Sample_description":
                sample_desc = [x.strip('"') for x in data]
            if label == "!series_matrix_table_begin":
                break
        csv_writer.writerow(["column_num", "sample_id", "donor_id", "region", "sex", "age", "stage", "sample_desc"])
        for i, (a, b, c, d, e, f, g) in enumerate(zip(sample_id, source_name, region, sex, age, stage, sample_desc)):
            csv_writer.writerow([i+1,a,b,c,d,e,f,g])
    def extract_row_info(self, File, OutName, platform):
        hand = open(File, 'rb')
        out = csv.writer(open(OutName, 'wb'), delimiter=",")
        flag_out = False
        plat = None
        for l in hand:
            if l.strip() == "^PLATFORM = "+platform:
                plat = platform
            if flag_out and plat == platform:
                out.writerow(l.split("\t"))
            if l.strip() == "!platform_table_begin":
                flag_out = True
            if l.strip() == "!platform_table_end":
                flag_out = False
                break
        hand.close()
        return
    def LookALLMutationTargetedGenesHBT(self, Genes, structure_acronyms, GeneDat, smooth=True, drop_low_exp=True, fontsize=6, ylim=None):
        GeneExp, GeneRow, GeneCol = GeneDat
        SelectedGenes = GeneRow[GeneRow["gene_symbol"].isin(Genes)]
        gene_ids = list(SelectedGenes.index)
        #print gene_ids
        plt.close('all')
        stages = {}
        for i, gene_id in enumerate(gene_ids):
            sys.stdout.write("\r{}".format(i))
            #print gene_id
            #print GeneExp.loc[gene_id, :]
            for structure_acronym in structure_acronyms:
                #print structure_acronym
                region_meta = GeneCol[GeneCol["region"]==structure_acronym]
                for index, row in region_meta.iterrows():
                    Period = row["stage"]
                    if Period not in stages:
                        stages[Period] = []
                    stages[Period].append(GeneExp.get_value(gene_id, row["sample_id"]))
        seq = []
        for stage in self.Stages:
            try:
                seq.append(np.mean(stages[stage]))
            except KeyError:
                seq.append(0)
        print seq
        if smooth:
            seq = smooth_func(seq)
        fig, ax = plt.subplots(dpi=200)
        plt.title("TargetedExons over Region:{}".format(",".join(structure_acronyms)))
        ax.plot(xrange(1,16),seq)
        ax.grid(True)
        ax.axvline(x=7.5)
        #plt.xticks(np.arange(2,14), Stages, rotation=20)
        #ax.xaxis.set_major_formatter(plt.FuncFormatter(format_func))
        if ylim != None:
            plt.ylim(ylim)
        plt.show()

    class DevStageExp:
        def __init__(self, DevStage):
            self.DevStage = DevStage
            self.exps = np.array([])
            self.mean_exp = 0
        def add_exp(self, value):
            self.exps = np.append(self.exps, float(value))
            #self.mean_exp = round(np.mean(self.exps), 4)
        def get_nonzero_mean(self):
            try:
                self.non_zero_exps = [x for x in self.exps if x]
                return sum(self.non_zero_exps)/len(self.non_zero_exps)
            except ZeroDivisionError:
                return 0

    # Return res1: sequence of exp in 12 dev stages res2: num of samples at each timepoint
