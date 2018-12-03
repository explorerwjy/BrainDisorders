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
import random
import bisect
import collections
import scipy.stats 

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

go_psd = "/Users/jiayao/Work/BrainDisorders/data/functions/go_synapse.list"
go_sig_skel = "/Users/jiayao/Work/BrainDisorders/data/functions/go_cyto_sig.list"
go_channel = "/Users/jiayao/Work/BrainDisorders/data/functions/go_channel.list"
go_chromatin = "/Users/jiayao/Work/BrainDisorders/data/functions/go_chromatin.list"

andy_psd_channel = "/Users/jiayao/Work/BrainDisorders/JW/Functional_Cluster_Gene_Lists/CLUST_NEURO_PSD_CHANN.txt"
andy_cytp_sign = "/Users/jiayao/Work/BrainDisorders/JW/Functional_Cluster_Gene_Lists/CLUST_NEURO_CYTO_SIGN.txt"
andy_snf_bromo = "/Users/jiayao/Work/BrainDisorders/JW/Functional_Cluster_Gene_Lists/CLUST_DEVEL_SNF_BROMO.txt"
andy_znf = "/Users/jiayao/Work/BrainDisorders/JW/Functional_Cluster_Gene_Lists/CLUST_DEVEL_ZN_FINGER.txt"


Stages = ["2A", "2B", "3A", "3B"] + map(str,range(4,12))

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
        self.median = 0
    def add_exp(self, value):
        self.exps = np.append(self.exps, float(value))
        #self.mean_exp = round(np.mean(self.exps), 4)
    def get_nonzero_mean(self):
        try:
            #self.non_zero_exps = [math.log(x, 10) for x in self.exps if x]
            self.non_zero_exps = [x for x in self.exps if x]
            N = len(self.non_zero_exps)
            return np.mean(self.non_zero_exps), math.sqrt(np.var(self.non_zero_exps))/math.sqrt(N), np.median(self.non_zero_exps)
        except ZeroDivisionError:
            return 0, 0, 0
class BrainSpan:
    def __init__(self):
        # store all exp value at a developemt time point
        #Stages = map(str,range(1,16))
        #self.Stages = ["1", "2A", "2B", "3A", "3B"] + map(str,range(4,12))
        self.Stages = ["2A", "2B", "3A", "3B"] + map(str,range(4,12))
        self.Descriptions = ["Early fetal", "Early fetal", "Early mid-fetal", "Early mid-fetal", "Late mid-fetal", "Late fetal", "Early infancy", "Late infancy", "Early childhood", "Late childhood", "Adolescence", "Adulthood"]
        self.Regions = []
        self.Regionsgt20 = ['OFC', 'VFC', 'HIP', 'ITC', 'AMY', 'DFC', 'STC', 'MFC', 'STR', 'IPC', 'V1C', 'S1C', 'A1C', 'M1C', 'CBC', 'MD']
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
    # Return res1: sequence of exp in 12 dev stages res2: num of samples at each timepoint
    def DisplayStageExp(self, Dict, out="nonzero"):
        res1 = []
        res2 = []
        res3 = []
        res4 = []
        for stage in self.Stages:
            try:
                #res1.append(Dict[stage].mean_exp) # Exp avg
                stage, Dict[stage].get_nonzero_mean()
                if out == "all":
                    res1.append(Dict[stage].mean_exp)
                else:
                    mean, std, median = Dict[stage].get_nonzero_mean()
                    res1.append(mean) # Exp avg
                    res3.append(std)
                    res4.append(median)
                res2.append(len(Dict[stage].exps)) # Num of samples at that time point
            except KeyError:
                res1.append(0)
                res2.append(0)
                res3.append(0)
                res4.append(0)
        return res1, res2, res3, res4
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
                #print "Three 0 data at stage",i
                new.append(0)
            i += 1
        return new
    def format_func(self, value, tick_number):
        try:
            return self.Stages[int(value)-2]
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
            seq = self.DisplayStageExp(stages)[0]
            if drop_low_exp:
                if seq.count(0) >= 4:
                    res[exon_id] = None
                    continue
            if smooth:
                res[exon_id] = self.smooth_func(seq)
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
        plt.xticks(np.arange(2,14), self.Stages, rotation=20)
        ax.xaxis.set_major_formatter(plt.FuncFormatter(self.format_func))
        plt.show()
    def LookGeneExonSumOverRegionwithTargetUntarget(self, Gene, structure_acronyms, bp_exon_row_meta, bp_exon_col_meta, ExonExp, GeneDat=None, smooth=True, drop_low_exp=True, space = 'log', fontsize=6):
        Exons = bp_exon_row_meta[bp_exon_row_meta["gene_symbol"]==Gene]
        TargetedExons = Exons[Exons["Vars"]!=""]
        UnTargetedExons = Exons[Exons["Vars"]==""]
        #print TargetedExons
        plt.close('all')
        fig, ax = plt.subplots(dpi=100)
        plt.title("Gene:{} over Region:{}".format(Gene, ",".join(structure_acronyms)))
        res = {}
        
        for exon_id in list(TargetedExons["row_num"]):
            exon_id = int(exon_id)
            seq, error_var, median = self.LoadingDat2SeqCrossRecordCrossRegion([exon_id], structure_acronyms, bp_exon_row_meta, bp_exon_col_meta, ExonExp, space = space,  smooth=True)
            ax.plot(range(2,14), seq, label=str(exon_id))
        for exon_id in list(UnTargetedExons["row_num"]):
            exon_id = int(exon_id)
            seq, error_var, median = self.LoadingDat2SeqCrossRecordCrossRegion([exon_id], structure_acronyms, bp_exon_row_meta, bp_exon_col_meta, ExonExp, smooth=True, space = space, shutup=True)
            ax.plot(range(2,14), seq, '--', alpha = 0.3, color='grey') 
        if GeneDat != None:# Plot Gene exp over stages 
            stages = {}
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
                        exp = GeneExp.get_value(gene_id-1, row["column_num"])
                        if space == 'log':
                            exp = math.log(exp+1, 2)
                        stages[Period].add_exp(exp)
                seq, Nsample, median, var = self.DisplayStageExp(stages)
                #print seq
                if smooth:
                    res[gene_id] = self.smooth_func(seq)
                else:
                    res[gene_id] = seq
            #add_layout(LinearAxis(y_range_name="GeneRPKM"), 'right')
            #ax2 = ax.twinx()
            for gene_id in gene_ids:
                #print res[gene_id]
                #ax2.plot(range(2,14), [math.log(x,2) if x !=0 else 0 for x in res[gene_id]], color="black") 
                #ax2.plot(range(2,14), res[gene_id], color="black") 
                ax.plot(range(2,14), res[gene_id], color="black") 
            #ax2.set_ylabel("GeneRPKM")
        ax.grid(True)
        ax.axvline(x=7.5)
        #Stages = Stages + [0,0,0,0,0]
        plt.xticks(np.arange(2,14), self.Stages, rotation=20)
        ax.xaxis.set_major_formatter(plt.FuncFormatter(self.format_func))
        plt.legend()
        ax.legend()
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
    # Assign Variants to row meta data file of Exons. 
    # Markers on 1)Male vs. Female; 2)IQ70+ vs. IQ70- 3)Functions
    def AssignVar2Exon2(self, bp_exon_row_meta, VarFile, IntersectionWithPredicted = True, ProSib="Pro"):
        VarDF = pd.read_excel(VarFile)
        # Load annotations to gene
        entrez_symbol_map = get_gene_entrez_symbol_map()
        wigler_predicted_lgd_genes = set([entrez_symbol_map[x.strip()] for x in file(wigler_predicted_lgd)])
        #vitkup_channel_genes = set([entrez_symbol_map[x.strip()] for x in file(vitkup_channel)])
        #vitkup_chromatin_genes = set([entrez_symbol_map[x.strip()] for x in file(vitkup_chromatin)])
        #vitkup_psd_genes = set([entrez_symbol_map[x.strip()] for x in file(vitkup_psd)])
        #vitkup_sig_skel_genes = set([entrez_symbol_map[x.strip()] for x in file(vitkup_sig_skel)])
        channel_genes = set([entrez_symbol_map.get(x.strip(), None) for x in file(go_channel)])
        chromatin_genes = set([entrez_symbol_map.get(x.strip(), None) for x in file(go_chromatin)])
        psd_genes = set([entrez_symbol_map.get(x.strip(), None) for x in file(go_psd)])
        sig_skel_genes = set([entrez_symbol_map.get(x.strip(), None) for x in file(go_sig_skel)])
        fam_info = pd.read_excel(wigler_fam_info)
        famID2NVIQ = dict(zip(list(fam_info["familyId"]), list(fam_info["probandNVIQ"])))
        famID2VIQ = dict(zip(list(fam_info["familyId"]), list(fam_info["probandVIQ"])))
        
        Cates = dict(zip(["channel", "chromatin", "psd", "sig_skel"],[channel_genes,chromatin_genes,psd_genes,sig_skel_genes]))
        res = {} # k:gene v:list of variants in that gene
        var_genders = {}
        var_IQ = {}
        for i, row in VarDF.iterrows():
            if IntersectionWithPredicted:
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
                if ProSib == 'Pro':
                    var_genders[var] = row["inChild"][:2]
                else:
                    var_genders[var] = row["inChild"][-2:]
                #var_genders[var] = [row["inChild"][:2]]
            else:
                if ProSib == 'Pro':
                    var_genders[var] = var_genders[var] + ";" + row["inChild"][:2]
                else:
                    var_genders[var] = var_genders[var] + ";" + row["inChild"][-2:]
            # Get IQ Info
            if var not in var_IQ:
                FamID = row["familyId"]
                var_IQ[var] = (famID2VIQ[FamID],famID2NVIQ[FamID])
        bp_exon_row_meta_rep = bp_exon_row_meta.copy(deep=True)
        bp_exon_row_meta_rep["NVIQ70"] = ""
        bp_exon_row_meta_rep["VIQ70"] = ""
        bp_exon_row_meta_rep["Gender"] = ""
        bp_exon_row_meta_rep["Func"] = ""
        bp_exon_row_meta_rep["Vars"] = ""
        bp_exon_row_meta_rep["GeneHited"] = "F"
        bp_exon_row_meta_rep["Last"] = "F"
        LAST_GENE = None
        for i, row in bp_exon_row_meta_rep.iterrows():
            sys.stdout.write("\r{}".format(i))
            gene = row["gene_symbol"]
            if gene != LAST_GENE and gene != None:
                bp_exon_row_meta_rep.at[i-1, "Last"] = "T" 
            LAST_GENE = gene
            if gene not in res:
                continue
            bp_exon_row_meta_rep.at[i, "GeneHited"] = "T"
            start, end = int(row["start"]) - 2 , int(row["end"]) + 2
            gene = row["gene_symbol"]
            # Assign Variant to Exon
            for var in res[gene]:
                pos = int(var.split(":")[1])
                #print pos
                if pos >= start and pos <= end: #Var In This Exon
                    bp_exon_row_meta_rep.at[i, "Vars"] = bp_exon_row_meta_rep.get_value(i, "Vars") + ";" + var if bp_exon_row_meta_rep.get_value(i, "Vars") != "" else var
                    bp_exon_row_meta_rep.at[i, "VIQ70"] = var_IQ[var][0]
                    bp_exon_row_meta_rep.at[i, "NVIQ70"] = var_IQ[var][1]
                    bp_exon_row_meta_rep.at[i, "Gender"] = bp_exon_row_meta_rep.get_value(i, "Gender") + ";" + var_genders[var] if bp_exon_row_meta_rep.get_value(i, "Gender") != "" else var_genders[var]
                    continue
            # Assign Other Info
            bp_exon_row_meta_rep.at[i, "Func"] = self.AssignFunc(gene, Cates)
        bp_exon_row_meta_rep.at[i, "Last"] = "T" #Last Exon in DF is last exon of that gene
        return bp_exon_row_meta_rep
    def AssignVar2Exon3(self, bp_exon_row_meta, VarFile, IntersectionWithPredicted = True, ProSib="Pro"):
        VarDF = pd.read_excel(VarFile)
        # Load annotations to gene
        entrez_symbol_map = get_gene_entrez_symbol_map()
        wigler_predicted_lgd_genes = set([entrez_symbol_map[x.strip()] for x in file(wigler_predicted_lgd)])
        channel_psd_genes = set([x.strip() for x in file(andy_psd_channel)])
        chromatin_genes = set([x.strip() for x in file(andy_snf_bromo)])
        cyto_sign = set([x.strip() for x in file(andy_cytp_sign)])
        tf_rbp = set([x.strip() for x in file(andy_znf)])
        fam_info = pd.read_excel(wigler_fam_info)
        famID2NVIQ = dict(zip(list(fam_info["familyId"]), list(fam_info["probandNVIQ"])))
        famID2VIQ = dict(zip(list(fam_info["familyId"]), list(fam_info["probandVIQ"])))
        
        Cates = dict(zip(["Channel&PSD", "Chromatin", "Cytp&Sig", "TF&RBP"],[channel_psd_genes,chromatin_genes,cyto_sign,tf_rbp]))
        res = {} # k:gene v:list of variants in that gene
        var_genders = {}
        var_IQ = {}
        for i, row in VarDF.iterrows():
            if IntersectionWithPredicted:
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
                if ProSib == 'Pro':
                    var_genders[var] = row["inChild"][:2]
                else:
                    var_genders[var] = row["inChild"][-2:]
                #var_genders[var] = [row["inChild"][:2]]
            else:
                if ProSib == 'Pro':
                    var_genders[var] = var_genders[var] + ";" + row["inChild"][:2]
                else:
                    var_genders[var] = var_genders[var] + ";" + row["inChild"][-2:]
            # Get IQ Info
            if var not in var_IQ:
                FamID = row["familyId"]
                var_IQ[var] = (famID2VIQ[FamID],famID2NVIQ[FamID])
        bp_exon_row_meta_rep = bp_exon_row_meta.copy(deep=True)
        bp_exon_row_meta_rep["NVIQ70"] = ""
        bp_exon_row_meta_rep["VIQ70"] = ""
        bp_exon_row_meta_rep["Gender"] = ""
        bp_exon_row_meta_rep["Func"] = ""
        bp_exon_row_meta_rep["Vars"] = ""
        bp_exon_row_meta_rep["GeneHited"] = "F"
        bp_exon_row_meta_rep["Last"] = "F"
        LAST_GENE = None
        for i, row in bp_exon_row_meta_rep.iterrows():
            sys.stdout.write("\r{}".format(i))
            gene = row["gene_symbol"]
            if gene != LAST_GENE and gene != None:
                bp_exon_row_meta_rep.at[i-1, "Last"] = "T" 
            LAST_GENE = gene
            if gene not in res:
                continue
            bp_exon_row_meta_rep.at[i, "GeneHited"] = "T"
            start, end = int(row["start"]) - 2 , int(row["end"]) + 2
            gene = row["gene_symbol"]
            # Assign Variant to Exon
            for var in res[gene]:
                pos = int(var.split(":")[1])
                #print pos
                if pos >= start and pos <= end: #Var In This Exon
                    bp_exon_row_meta_rep.at[i, "Vars"] = bp_exon_row_meta_rep.get_value(i, "Vars") + ";" + var if bp_exon_row_meta_rep.get_value(i, "Vars") != "" else var
                    bp_exon_row_meta_rep.at[i, "VIQ70"] = var_IQ[var][0]
                    bp_exon_row_meta_rep.at[i, "NVIQ70"] = var_IQ[var][1]
                    bp_exon_row_meta_rep.at[i, "Gender"] = bp_exon_row_meta_rep.get_value(i, "Gender") + ";" + var_genders[var] if bp_exon_row_meta_rep.get_value(i, "Gender") != "" else var_genders[var]
                    continue
            # Assign Other Info
            bp_exon_row_meta_rep.at[i, "Func"] = self.AssignFunc(gene, Cates)
        bp_exon_row_meta_rep.at[i, "Last"] = "T" #Last Exon in DF is last exon of that gene
        return bp_exon_row_meta_rep
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
    def LoadingDat2SeqWithLenCrossRecordCrossRegion(self, record_ids, regions, row_meta, col_meta, Matrix, smooth=True, drop_low_exp=True):
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
    def LoadingNormDat2SeqCrossRecordCrossRegion(self, record_ids, regions, row_meta, col_meta, Matrix, smooth=True, drop_low_exp=True):
        res = {}
        ALL_RPKM_X_L = [0]*13 # by Timepoint
        ALL_L = 0
        Numrecords = len(record_ids)
        for i, _id in enumerate(record_ids):
            sys.stdout.write("\rLoading records:{}".format(i))
            #lengths.append(row_meta.get_value(int(_id)-1, "exon length"))
            exonL = int(row_meta.get_value(int(_id)-1, "exon length"))
            stages = {} #Exom RPKM at each stage
            for j, region in enumerate(regions):
                region_col_meta = col_meta[col_meta["structure_acronym"]==region]
                for index, row in region_col_meta.iterrows():
                    Period = row["Period"]
                    if Period not in stages:
                        stages[Period] = DevStageExp(Period)
                    stages[Period].add_exp(Matrix.get_value(_id-1, row["column_num"]))
            seq = self.DisplayStageExp(stages)[0]
            if drop_low_exp:
                if seq.count(0) >= 4:
                    print _id
                    Numrecords -= 1
                    continue
            #print seq
            res[_id] = self.smooth_func(seq) if smooth else seq
            ALL_RPKM_X_L = [x+(y*exonL) for x,y in zip(ALL_RPKM_X_L, res[_id])]
            ALL_L += exonL
            #print ALL_RPKM_X_L, ALL_L
        records_mean = []
        records_mean_norm = []
        denominators = []
        for i, stage in enumerate(self.Stages):
            tmp1 = 0
            tmp2 = 0
            denominator = ALL_RPKM_X_L[i] / ALL_L
            denominators.append(denominator)
            expsi = []
            for j, _id in enumerate(record_ids):
                if _id in res:
                    #pprint stage, res[_id][i]
                    tmp1 += res[_id][i]
                    res[_id][i] = res[_id][i] / denominator
                    tmp2 += res[_id][i]
                    expsi.append(res[_id][i])
            print np.var(expsi),
            records_mean.append(tmp1/Numrecords)
            records_mean_norm.append(tmp2/Numrecords)
            print
        print records_mean, records_mean_norm, denominators
        plt.close('all')
        fig, ax = plt.subplots()
        ax.plot(xrange(2,14), records_mean, label='avg')
        ax.plot(xrange(2,14), denominators, label='norm term')
        ax2 = ax.twinx()
        ax2.plot(xrange(2,14), records_mean_norm, label='avg_nromed')
        ax.legend()
        ax2.legend()
        plt.show()
        sys.stdout.write("\n")
        print records_mean 
        return records_mean 
    def LoadingDat2SeqCrossRecordCrossRegion(self, record_ids, regions, row_meta, col_meta, Matrix, smooth=True, drop_low_exp=True, space = 'norm', shutup=True):
        stages = {}
        for i, _id in enumerate(record_ids):
            if not shutup:
                sys.stdout.write("\rLoading records:{}".format(i))
            for j, region in enumerate(regions):
                region_col_meta = col_meta[col_meta["structure_acronym"]==region]
                for index, row in region_col_meta.iterrows():
                    Period = row["Period"]
                    if Period not in stages:
                        stages[Period] = DevStageExp(Period)
                    exp = Matrix.get_value(_id-1, row["column_num"])
                    if space == 'log':
                        exp = math.log(exp+1, 2)
                    stages[Period].add_exp(exp)
        if not shutup:
            sys.stdout.write("\n")
        seq, Nsample, error_var, median = self.DisplayStageExp(stages)
        if smooth:
            seq = self.smooth_func(seq)
        return seq, error_var, median 
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
        #ax.legend()
        plt.show()
        return ExonLengths
    def LookALLMutationTargetedExon(self, exon_ids_sets, structure_acronyms, bp_exon_row_meta, bp_exon_col_meta, ExonExp, smooth=True, drop_low_exp=True, fontsize=6, rd=False, method = "mean", space = 'log'):
        plt.close('all')
        fig = plt.figure(dpi=800)
        stages = {}
        fig, ax = plt.subplots(dpi=200)
        plt.title("Exons Exp over Region:{}".format(",".join(structure_acronyms)))
        seq1, seq2 = None, None

        for SetName, (color, ExonIDs) in exon_ids_sets.items():
            print SetName
            ExonIDs = [int(x) for x in ExonIDs]
            seq, error_var, median = self.LoadingDat2SeqCrossRecordCrossRegion(ExonIDs, structure_acronyms, bp_exon_row_meta, bp_exon_col_meta, ExonExp, smooth=False, space = space, shutup=False) # Loading data needed
            #seq = [math.log(x, 10) for x in seq]
            if "Other" in SetName:
                if method == "mean":
                    ax.errorbar(range(2,14), seq, yerr = error_var, linestyle = '--', label=SetName, color=color)
                    seq2 = seq
                elif method == "median":
                    ax.plot(range(2,14), median, '--', label=SetName, color=color)
                    seq2 = median
            else:
                if method == "mean":
                    ax.errorbar(range(2,14), seq, yerr = error_var, label=SetName, color=color)
                    seq1 = seq
                elif method == "median":
                    ax.plot(range(2,14),median, label=SetName, color=color)
                    seq1 = median
        ax.grid(True)
        ax.axvline(x=7.5)
        plt.xticks(np.arange(2,14), self.Stages, rotation=20)
        ax.xaxis.set_major_formatter(plt.FuncFormatter(self.format_func))
        ax.legend()
        plt.show()
        pre, post = self.Bias(seq1, seq2)
        print pre, post, pre-post
    # with Nromlization within exon sets.
    def LookALLMutationTargetedExon(self, exon_ids_pairs, structure_acronyms, bp_exon_row_meta, bp_exon_col_meta, ExonExp, smooth=True, drop_low_exp=True, fontsize=6, rd=False, method = "mean", space = 'log'):
        plt.close('all')
        fig = plt.figure(dpi=800)
        stages = {}
        fig, ax = plt.subplots(dpi=200)
        plt.title("Exons Exp over Region:{}".format(",".join(structure_acronyms)))
        for (SetName, color, TExonIDs, UExonIDs) in exon_ids_pairs:
            seq1, seq2 = None, None
            TExonIDs = [int(x) for x in TExonIDs]
            UExonIDs = [int(x) for x in UExonIDs]
            Tseq, Terror_var, Tmedian = self.LoadingDat2SeqCrossRecordCrossRegion(TExonIDs, structure_acronyms, bp_exon_row_meta, bp_exon_col_meta, ExonExp, smooth=False, space = space, shutup=False) # Loading data needed
            Useq, Uerror_var, Umedian = self.LoadingDat2SeqCrossRecordCrossRegion(UExonIDs, structure_acronyms, bp_exon_row_meta, bp_exon_col_meta, ExonExp, smooth=False, space = space, shutup=False) # Loading data needed
            Tseq = [2 ** x for x in Tseq]
            Useq = [2 ** x for x in Useq]
            #Terror_var = self.converterror(Tseq, Terror_var)
            #Uerror_var = self.converterror(Useq, Uerror_var)
            #seq = [math.log(x, 10) for x in seq]
            if method == "mean":
                ax.errorbar(range(2,14), Tseq, yerr = Terror_var, label=SetName+"_Targeted", color=color)
                ax.errorbar(range(2,14), Useq, yerr = Uerror_var, linestyle = '--', label=SetName+"_Untargeted", color=color)
                seq1 = Tseq
                seq2 = Useq
            elif method == "median":
                ax.plot(range(2,14), Tmedian, '--', label=SetName, color=color)
                ax.plot(range(2,14), Umedian, '--', label=SetName, color=color)
                seq1 = Tmedian
                seq2 = Umedian
            pre, post = self.Bias(seq1, seq2)
            print SetName, pre, post, pre-post
        ax.grid(True)
        ax.axvline(x=7.5)
        plt.xticks(np.arange(2,14), self.Stages, rotation=20)
        ax.xaxis.set_major_formatter(plt.FuncFormatter(self.format_func))
        ax.legend()
        plt.show()
    def converterror(self, means, stderrs):
        upper, lower = [], []
        for mean, error in zip(means, stderrs):
            upper.append( (2**(mean+error)) - (2**(mean)) )
            lower.append( (2**(mean)) - (2**(mean-error)) )
        #print upper, lower
        return upper, lower
    # with Nromlization within exon sets.
    def ReadExonExpValuesEachTimePointAvgAcrossRegion(self, exon_ids, structure_acronyms, bp_exon_row_meta, bp_exon_col_meta, ExonExp, smooth=True, drop_low_exp=True):
        res = {}
        for ExonID in exon_ids:
            seq = self.LoadingDat2SeqCrossRecordCrossRegion([ExonID], structure_acronyms, bp_exon_row_meta, bp_exon_col_meta, ExonExp, smooth=True) # Loading data needed
            res[ExonID] = seq
        return res
    def LookALLMutationTargetedExon3(self, exon_ids_sets, structure_acronyms, bp_exon_row_meta, bp_exon_col_meta, ExonExp, smooth=True, drop_low_exp=True, fontsize=6, rd=False):
        plt.close('all')
        fig = plt.figure(dpi=800)
        stages = {}
        fig, ax = plt.subplots(dpi=200)
        plt.title("Exons Exp over Region:{}".format(",".join(structure_acronyms)))
        for SetName, (color, ExonIDs) in exon_ids_sets.items():
            print SetName
            ExonIDs = [int(x) for x in ExonIDs]
            seq = self.LoadingNormDat2SeqCrossRecordCrossRegion(ExonIDs, structure_acronyms, bp_exon_row_meta, bp_exon_col_meta, ExonExp, smooth=True) # Loading data needed
            seq = [math.log(x,10) for x in seq]
            if "Other" in SetName:
                ax.plot(range(2,14), seq, '--', label=SetName, color=color)
                #ax.plot(range(2,14),[math.log(x, 2) if x!=0 else 0 for x in seq], '--', label=SetName, color=color)
            else:
                ax.plot(range(2,14), seq, label=SetName, color=color)
                #ax.plot(range(2,14),[math.log(x, 2) if x!=0 else 0 for x in seq], label=SetName, color=color)
        ax.grid(True)
        ax.axvline(x=7.5)
        plt.xticks(np.arange(2,14), self.Stages, rotation=20)
        ax.xaxis.set_major_formatter(plt.FuncFormatter(self.format_func))
        ax.legend()
        plt.show()
    def LoadSummarizedExonExp(self, record_ids, bp_exon_row_meta, ExonExp, drop_low_exp = True):
        stages = {}
        drops = []
        for stage in self.Stages:
            stages[stage] = []
        N = 0
        Ndrop = 0
        for i, _id in enumerate(record_ids):
            dat = ExonExp.iloc[_id-1, 1:13].values
            if drop_low_exp and dat.mean() < 1 :
                Ndrop += 1
                drops.append(_id)
                continue
            for j, stage in enumerate(self.Stages):
                #stages[stage].append(dat[j])
                stages[stage].append(math.log(dat[j]+1, 2))
            N += 1
        log2seq = []
        stderrs = []
        medians = []
        for i, stage in enumerate(self.Stages):
            log2seq.append(np.mean(stages[stage]))
            stderrs.append(math.sqrt( np.var(stages[stage])/N ))
            medians.append(np.median(stages[stage]))
        print Ndrop
        return log2seq, stderrs, medians, drops
    def Bias(self, seq1, seq2):
        diff = []
        for i, stage in enumerate(self.Stages):
            diff.append(seq1[i]-seq2[i])
        pre, post = sum(diff[0:6]), sum(diff[6:12])
        return pre, post
    def LookALLMutationTargetedExon4(self, title, exon_ids_sets, bp_exon_row_meta, ExonExp, method = "mean", smooth=True, drop_low_exp=True, fontsize=6):
        plt.close('all')
        fig, ax = plt.subplots(dpi=120)
        plt.title("Exons Exp {}".format(title))
        seq1, seq2 = None, None
        ALL_drop = []
        for SetName, (color, ExonIDs) in exon_ids_sets.items():
            print SetName
            ExonIDs = [int(x) for x in ExonIDs]
            log2seq, stderrs, medians, drops = self.LoadSummarizedExonExp(ExonIDs, bp_exon_row_meta, ExonExp, drop_low_exp=drop_low_exp) # Loading data needed
            ALL_drop.extend(drops)
            if "Other" in SetName:
                if method == "median":
                    ax.plot(range(2,14), medians, '--', label=SetName, color=color)
                    seq2 = medians
                else:
                    ax.errorbar(range(2,14), log2seq, yerr = stderrs, linestyle = '--', label=SetName, color=color)
                    seq2 = log2seq
            else:
                if method == "median":
                    ax.plot(range(2,14), medians, label=SetName, color=color)
                    seq1 = medians
                else:
                    ax.errorbar(range(2,14), log2seq , yerr = stderrs, label=SetName, color=color)
                    seq1 = log2seq
        ax.grid(True)
        ax.axvline(x=7.5)
        plt.xticks(np.arange(2,14), self.Stages, rotation=20)
        ax.xaxis.set_major_formatter(plt.FuncFormatter(self.format_func))
        ax.legend()
        plt.show()
        pre, post = self.Bias(seq1, seq2)
        print pre, post, pre-post
        return ALL_drop

    def LookALLMutationTargetedGenes(self, Gene_id_sets, structure_acronyms, GeneDat, smooth=True, drop_low_exp=True, fontsize=6, ylim=None):
        plt.close('all')
        GeneExp, GeneRow, GeneCol = GeneDat
        plt.close('all')
        stages = {}
        fig, ax = plt.subplots(dpi=200)
        plt.title("TargetedExons over Region:{}".format(",".join(structure_acronyms)))
        for SetName, (color, GeneIDs) in Gene_id_sets.items():
            print SetName
            GeneIDs = [int(x) for x in GeneIDs]
            seq = self.LoadingDat2SeqCrossRecordCrossRegion(GeneIDs, structure_acronyms, GeneCol, GeneCol, GeneExp, smooth=True) # Loading data needed
            #ax.plot(range(1,16),[math.log(x, 2) if x!=0 else 0 for x in seq])
            ax.plot(range(2,14), seq, label=SetName, color=color)
        ax.grid(True)
        ax.axvline(x=7.5)
        plt.legend()
        #plt.xticks(np.arange(2,14), Stages, rotation=20)
        ax.xaxis.set_major_formatter(plt.FuncFormatter(self.format_func))
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

    def LoadGeneSetData2Fil(self, FilName, Genes, Regions, Exon_row_meta, Exon_col_meta, Matrix):
        #Regions_indice
        writer = csv.writer(open(FilName, 'wb'))
        print "Total Num of Genes:", len(Genes)
        for i, gene in enumerate(Genes):
            ExonIDs = Exon_row_meta[Exon_row_meta["gene_symbol"]==gene]["row_num"].values
            for ExonID in ExonIDs:
                ExonID = int(ExonID)
                for i, stage in enumerate(self.Stages):
                    tmp = [gene, ExonID, stage]
                    cols = Exon_col_meta[(Exon_col_meta["Period"]==stage) & (Exon_col_meta["structure_acronym"].isin(Regions))]["column_num"].values
                    for col in cols:
                        exp = Matrix.get_value(ExonID-1, col)
                        exp = round(math.log(exp+1, 2), 6)
                        tmp.append(exp)
                    writer.writerow(tmp)
            #if i % 10 == 0:
            #sys.stdout.write("\rLoad {} genes".format(i))

    def LoadGeneSetDataFromFil(self, FilName):
        res = {}
        fin = open(FilName, 'rb')
        for l in fin:
            llist = l.strip().split(",")
            Gene, Exon, time = llist[:3]
            Exon = int(Exon)
            exps = map(float, llist[3:])
            if Gene not in res:
                res[Gene] = {}
            if Exon not in res[Gene]:
                res[Gene][Exon] = {}
            if time not in res[Gene][Exon]:
                res[Gene][Exon][time] = None
            res[Gene][Exon][time] = exps
        return res

class GeneExonSet(BrainSpan):
    def __init__(self, expdict, Name="geneset", Color='black'):
        self.GeneSetName = Name
        self.GeneSetColor = Color
        self.Genes = {} # a list of genes involved
        self.expdict = expdict
        self.genes = []
        self.TargetedExon = []
        self.UntargetedExon = []
    def addGene(self, Genesymbol, TargetedExon, UntargetedExon, TargetedExonLength, UntargetedExonLength):
        self.Genes[Genesymbol] = TargetedGene(Genesymbol, TargetedExon, UntargetedExon, TargetedExonLength, UntargetedExonLength)
        self.genes.append(Genesymbol)
        self.TargetedExon.extend(TargetedExon)
        self.UntargetedExon.extend(UntargetedExon)
        #self.Genes[Genesymbol].addExp(ExpDict)

    def Reduce(self, logscale = False):
        self.Tseq, self.Terr, self.Useq, self.Uerr = [], [], [], []
        for stage in Stages:
            tmp1 = []
            tmp2 = []
            for gene in self.Genes:
                for ExonID in self.Genes[gene].TargetedExons:
                    tmp1.extend(self.expdict[gene][ExonID][stage])
                for ExonID in self.Genes[gene].UntargetedExons:
                    tmp2.extend(self.expdict[gene][ExonID][stage])
            self.Tseq.append(np.mean(tmp1))
            self.Terr.append( math.sqrt( np.var(tmp1) / len(tmp1) ) )
            self.Useq.append(np.mean(tmp2))
            self.Uerr.append( math.sqrt( np.var(tmp2) / len(tmp2) ) )
        if logscale == True:
            Tseq, Useq, Terr, Uerr = self.Tseq, self.Useq, self.Terr, self.Uerr
        else:
            Tseq, Useq, Terr, Uerr = [2 ** x for x in self.Tseq], [2 ** x for x in self.Useq], self.converterror(self.Tseq, self.Terr), self.converterror(self.Useq, self.Uerr)
        return (Tseq, Terr, Useq, Uerr)
    def Reduce2(self):
        self.Tseq, self.Terr, self.Useq, self.Uerr = [], [], [], []
        for stage in Stages:
            tmp1 = []
            tmp2 = []
            for gene in self.Genes:
                for ExonID in self.Genes[gene].TargetedExons:
                    tmp1.extend(self.expdict[gene][ExonID][stage])
                for ExonID in self.Genes[gene].UntargetedExons:
                    tmp2.extend(self.expdict[gene][ExonID][stage])
            self.Tseq.append(np.mean(tmp1))
            self.Terr.append( math.sqrt( np.var(tmp1) / len(tmp1) ) )
            self.Useq.append(np.mean(tmp2))
            self.Uerr.append( math.sqrt( np.var(tmp2) / len(tmp2) ) )

    def Plot(self, title="", dpi=200, color='black'):
        plt.close('all')
        fig, ax = plt.subplots(dpi=dpi)
        plt.title(title)
        Tseq, Useq, Terr, Uerr = [2 ** x for x in self.Tseq], [2 ** x for x in self.Useq], self.converterror(self.Tseq, self.Terr), self.converterror(self.Useq, self.Uerr)
        #Tseq, Useq, Terr, Uerr = self.Tseq, self.Useq, self.Terr, self.Uerr
        ax.errorbar(range(2,15), Tseq, yerr = Terr, label=title+"_Targeted", color=color)
        ax.errorbar(range(2,15), Useq, yerr = Uerr, linestyle = '--', label=title+"_Untargeted", color=color)
        pre, post = self.Bias(Tseq, Useq)
        print pre, post, pre-post
        ax.grid(True)
        ax.axvline(x=7.5, color="black")
        plt.xticks(np.arange(2,15), Stages, rotation=20)
        #ax.xaxis.set_major_formatter(plt.FuncFormatter(self.format_func))
        ax.legend()
        plt.show()
    def converterror(self, means, stderr):
        upper, lower = [], []
        for mean, error in zip(means, stderr):
            upper.append( (2**(mean+error)) - (2**(mean)) )
            lower.append( (2**(mean)) - (2**(mean-error)) )
        #print upper, lower
        return upper, lower
    def wilcoxonTest(self):
        self.genePreBias = []    
        self.genePostBias = []    
        for g, gene in self.Genes.items():
            pre, post = gene.reduce(self.expdict)
            self.genePreBias.append(pre)
            self.genePostBias.append(post)
        return self.genePreBias, self.genePostBias

    def Permute(self, plot=False):
        Tseq, Useq = [], []
        tmp1, tmp2 = {}, {}
        TExons, UExons = [],[]
        for Gene in self.genes:
            GeneTargetedExons, GeneUntargetedExons = self.Genes[Gene].Permute()
            #print Gene, len(GeneTargetedExons), len(GeneUntargetedExons)
            #print GeneTargetedExons, GeneUntargetedExons
            TExons.extend(GeneTargetedExons)
            UExons.extend(GeneUntargetedExons)
            for stage in Stages:
                if stage not in tmp1:
                    tmp1[stage] = []
                if stage not in tmp2:
                    tmp2[stage] = []
                for ExonID in GeneTargetedExons:
                    tmp1[stage].extend(self.expdict[Gene][ExonID][stage])
                for ExonID in GeneUntargetedExons:
                    tmp2[stage].extend(self.expdict[Gene][ExonID][stage])
        for stage in Stages:
            Tseq.append(np.mean(tmp1[stage]))
            Useq.append(np.mean(tmp2[stage]))
        Tseq, Useq = [2 ** x for x in Tseq], [2 ** x for x in Useq]
        pre, post, bias = Bias(Tseq, Useq)
        if plot:
            fig, ax = plt.subplots()
            ax.plot(range(2,14), Tseq, color='red')
            ax.plot(range(2,14), Useq, color='blue')
            plt.show()
        return bias, TExons, UExons
        
def Bias(seq1, seq2, method=3):
    if method == 3:
        FC_pre = []
        for T, U in zip(seq1[:6], seq2[:6]):
            FC_pre.append(T/U)
        FC_post = []
        for T, U in zip(seq1[6:12], seq2[6:12]):
            FC_post.append(T/U)
        FC_pre, FC_post = np.mean(FC_pre), np.mean(FC_post)
        return FC_pre, FC_post, FC_pre / FC_post
    else:
        diff = []
        for i, stage in enumerate(Stages):
            diff.append(seq1[i]-seq2[i])
        #pre, post = sum(diff[0:6])/6, sum(diff[6:12])/6
        pre, post = sum(diff[0:6]), sum(diff[6:12])
        #return pre, post, (pre-post)/post
        return pre, post, (pre-post)

class TargetedGene:
    def __init__(self, genesymbol, TargetedExonID, UntargetedExonID, TargetedExonLength, UntargetedExonLength):
        self.genesymbol = genesymbol
        self.TargetedExons = TargetedExonID
        self.NTargeted = len(TargetedExonID)
        self.UntargetedExons = UntargetedExonID
        self.TargetedExonLength = TargetedExonLength
        self.UntargetedExonLength = UntargetedExonLength
        self.AllExons = self.TargetedExons + self.UntargetedExons
        self.AllLength = self.TargetedExonLength + self.UntargetedExonLength
        #self.AllLength = [math.log(x, 1000) for x in self.AllLength]
        #self.AllLength = [min(500, x) for x in self.AllLength]
        self.TotalLength = sum(self.AllLength) 
        self.Weights = [float(x)/self.TotalLength for x in self.AllLength]
        #self.Weights = [1 for x in self.AllLength]
        #print self.Weights
        self.cdf_vals = cdf(self.Weights)
    def Permute(self):
        res_tar, res_untar = set([]), set([])
        i = 0
        while i < self.NTargeted:
            x = random.random()
            idx = bisect.bisect(self.cdf_vals, x)
            TheOne = self.AllExons[idx]
            #tmp = self.AllExons[::]
            #random.shuffle(tmp)
            #TheOne = tmp[0]
            #if TheOne not in res_tar:
            res_tar.add(TheOne)
            i += 1 
            #else:
            #    continue
        res_untar = set(self.AllExons).difference(res_tar)
        return list(res_tar), list(res_untar)
    def reduce(self, expdict):
        res = []
        for stage in Stages:
            u = np.mean([expdict[self.genesymbol][x][stage] for x in self.TargetedExons])
            l = np.mean([expdict[self.genesymbol][x][stage] for x in self.UntargetedExons])
            #u = 2**u
            #l = 2**l
            res.append(u-l)
            #res.append(2**(u-l))
        #print res
        if np.mean(res[:6]) - np.mean(res[6:12]) < -10:
            print self.genesymbol, res 
        return np.mean(res[:6]), np.mean(res[6:12])
    def plot(self, expdict):
        fig, ax = plt.subplots(dpi=120)
        plt.title(self.genesymbol)
        Tseq, Useq, Uerr = [], [], [] 
        for stage in Stages:
            u = np.mean([expdict[self.genesymbol][x][stage] for x in self.TargetedExons])
            l = np.mean([expdict[self.genesymbol][x][stage] for x in self.UntargetedExons])
            l_err = math.sqrt(np.var([expdict[self.genesymbol][x][stage] for x in self.UntargetedExons])/len([expdict[self.genesymbol][x][stage] for x in self.UntargetedExons]))
            Uerr.append(l_err)
            #u = 2**u
            #l = 2**l
            Tseq.append(u)
            Useq.append(l)
            #res.append(2**(u-l))
        ax.plot(range(2,14), Tseq, label="Targeted", color='red')
        ax.errorbar(range(2,14), Useq, yerr = Uerr, linestyle = '--', label="Untargeted", color='grey')
        pre, post, bias = Bias(Tseq, Useq)
        print pre, post, bias
        ax.grid(True)
        ax.axvline(x=7.5, color="grey", linestyle="--")
        plt.xticks(np.arange(2,14), Stages, rotation=20)
        ax.legend(loc='upper right')
        plt.xlabel("Dev Stages")
        plt.ylabel("Expression")
        plt.show()


class TargetedGeneWithData:
    def __init__(self, genesymbol, TargetedExonID, UntargetedExonID):
        self.genesymbol = genesymbol
        self.TargetedExons = TargetedExonID
        self.NTargeted = len(TargetedExonID)
        self.UntargetedExons = UntargetedExonID
        #self.AllExons = [x for (x,y) in self.TargetedSet] + [x for (x,y) in self.UntargetedSet]
        #self.TotalLength = sum(y for (x,y) in self.TargetedSet) + sum(y for (x,y) in self.UntargetedSet)
        #self.Probs = dict(zip(self.AllExons, [float(y/self.totalLength) for (x,y) in self.TargetedExonID + self.UntargetedExonID]))
        #self.cdf_vals = cdf(weights)

    def addExp(self, ExpDict):
        for i, exonid, exonlength in enumerate(self.TargetedSet):
            self.Targetedset[i] = (exonid, exonlength, ExpDict[exonid])
        UntargetedSet = []
        for exonid in self.UntargetedSet:
            Untargetedset.append( (exonid, exonlegnth, ExpDict[exonid]) )
        self.TargetedSet = TargetedSet
        self.UntargetedSet = UntargetedSet
    def Permute(self):
        selected = set([])
        res_tar, res_untar = [], []
        i = 0
        while i < self.NTargeted:
            x = random.random()
            idx = bisect.bisect(self.cdf_vals, x)
            if x not in selected:
                selected.add(x)
                i += 1 
            else:
                continue

def cdf(weights):
    total = sum(weights)
    result = []
    cumsum = 0
    for w in weights:
        cumsum += w
        result.append(cumsum / total)
    return result
def choice(population, weights):
    assert len(population) == len(weights)
    cdf_vals = cdf(weights)
    x = random.random()
    idx = bisect.bisect(cdf_vals, x)
    return population[idx]
def GetExonProb(exon_df):
    Total_length = sum(exon_df["exon length"].values)
    Probs = [float(x)/Total_length for x in exon_df["exon length"].values]
    return exon_df.index.values, Probs

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
                return np.mean(self.non_zero), np.var(self.non_zero), np.median(self.non_zero_exps)
            except ZeroDivisionError:
                return 0, 0, 0

from sklearn import datasets, linear_model
from sklearn.metrics import mean_squared_error, r2_score
import scipy
import random
import bisect
import collections
class PredictIQFromExon:
    def __init__(self, VarFil, FamINFO, matrix):
        self.matrix = matrix
        self.df_var = pd.read_excel(VarFil)
        self.df_fam = pd.read_excel(FamINFO)
    def sameExon(self, var1, var2, gene, row_meta):
        gene_df = row_meta[row_meta["gene_symbol"]==gene]
        pos_1 = int(var1.split(":")[1])
        pos_2 = int(var2.split(":")[1])
        exon_1, exon_2 = None, None
        for i, row in gene_df.iterrows():
            start = int(row["start"]) - 2
            end = int(row["end"]) + 2
            if pos_1 >= start and pos_1 <= end:
                exon_1 = i
            if pos_2 >= start and pos_2 <= end:
                exon_2 = i
        return (exon_1 == exon_2) and (exon_1 != None) and (exon_2 != None)
    def get_exon_id(self, var, gene, row_meta):
        gene_df = row_meta[row_meta["gene_symbol"]==gene]
        pos = int(var.split(":")[1])
        for i, row in gene_df.iterrows():
            start = int(row["start"]) - 2
            end = int(row["end"]) + 2
            if pos >= start and pos <= end:
                return i
        return None
    def isLastExon(self, exon_id, row_meta):
        gene = row_meta.get_value(exon_id, "gene_symbol")
        next_gene = row_meta.get_value(exon_id+1, "gene_symbol")
        return not (gene == next_gene)
        def cdf(self, weights):
            total = sum(weights)
            result = []
            cumsum = 0
            for w in weights:
                cumsum += w
                result.append(cumsum / total)
            return result
    def choice(self, population, weights):
        assert len(population) == len(weights)
        cdf_vals = cdf(weights)
        x = random.random()
        idx = bisect.bisect(cdf_vals, x)
        return population[idx]
    def GetExonProb(self, exon_df):
        Total_length = sum(exon_df["exon length"].values)
        Probs = [float(x)/Total_length for x in exon_df["exon length"].values]
        return exon_df.index.values, Probs
    def avgseq(self, seq, start, end):
        return np.mean(seq[start: end])
    def getIQ(self, famID, which="probandNVIQ"):
        return float(self.df_fam[self.df_fam["familyId"] == famID][which].values[0])
    def regGene(self, dat):
        regr = linear_model.LinearRegression(fit_intercept=False)
        X_train = np.array([[x[4]] for x in dat])
        y_train = np.array([x[3] for x in dat])
        regr.fit(X_train, y_train)
        return regr
    def regGene2(self, X, Y):
        regr = linear_model.LinearRegression(fit_intercept=False)
        X = X.reshape(-1,1)
        regr.fit(X, Y)
        return regr
    def plotGene(self, gene, regr, dat):
        X_train = np.array([[x[4]] for x in dat])
        y_train = np.array([x[3] for x in dat])
        y_pred = regr.predict(X_train)
        #print('Coefficients:', regr.coef_)
        R, P = scipy.stats.pearsonr([x[0] for x in X_train], y_train)
        plt.title("Gene:{}, R:{}, P:{}".format(gene, R,P))
        plt.xlim(0,3)
        plt.ylim(0, max(y_train+y_pred))
        plt.scatter(X_train, y_train,  color='black')
        plt.plot(np.append(X_train, [0]), np.append(y_pred, 0), color='blue', linewidth=3)
        plt.show()
    def model(self, df_row, var2leave = None, plot=True):
        self.gene2slope = {}
        ALL_Exp = np.array([])
        ALL_Normed_IQD = np.array([])
        Genes = list(set(self.df_var["effectGene"].values))
        for gene in Genes:
            tmp_df = self.df_var[self.df_var["effectGene"]==gene]
            famids = tmp_df["familyId"].values
            varids = tmp_df["location"].values
            if var2leave != None:
                varids = [x for x in varids if x!=var2leave]
            exonids = [self.get_exon_id(x, gene, df_row) for x in varids]
            exonids = [x for x in exonids if x != None]
            #print exonids
            exonids = [x for x in exonids if not self.isLastExon(x, df_row)]
            exps = []
            for exonid in exonids:
                seq = self.matrix.iloc[exonid-1, 0:12].values
                #print seq
                exp = self.avgseq(seq, 0, 12) / 10
                exps.append(exp)
            IQs = [max(0, (100-self.getIQ(x))) for x in famids]
            #exonids = [x for x in exonids if x!=None]
            dat = zip(famids, varids, exonids, IQs, exps)
            dat = [x for x in dat if x[2]!= None]
            exps = np.array([[x[4]] for x in dat])
            IQDs = np.array([x[3] for x in dat])
            if dat == []:
                continue
            regr = self.regGene(dat)
            slope = regr.coef_[0]
            if slope == 0:
                continue
            self.gene2slope[gene] = slope
            ALL_Exp = np.append(ALL_Exp, exps)
            NormIQDs = np.array([x[3]/slope for x in dat])
            ALL_Normed_IQD = np.append(ALL_Normed_IQD, NormIQDs)
        model = linear_model.LinearRegression(fit_intercept=False)
        ALL_Exp = ALL_Exp.reshape(-1,1)
        R, P = scipy.stats.pearsonr([x[0] for x in ALL_Exp], ALL_Normed_IQD)
        if plot:
            plt.title("R:{}, P:{}".format(R,P))
            model.fit(ALL_Exp, ALL_Normed_IQD)
            plt.scatter(ALL_Exp, ALL_Normed_IQD)
            upper = max(ALL_Normed_IQD)
            plt.plot([0, upper], [0,model.predict(upper)], color='red')
            plt.show()
        else:
            print "R:{}, P:{}".format(R,P) 
    def LoadData(self, df_row, method = 'rolling_time', time_start=0, time_end=12):
        assert time_end > time_start
        self.time_start, self.time_end = time_start, time_end
        Genes = list(set(self.df_var["effectGene"].values))
        self.VarList = {}
        self.Gene2Var = {}
        self.Genes = []
        for gene in Genes:
            tmp_df = self.df_var[self.df_var["effectGene"]==gene]
            famids = tmp_df["familyId"].values
            varids = tmp_df["location"].values
            exonids = [self.get_exon_id(x, gene, df_row) for x in varids]
            exonids = [x for x in exonids if x != None]
            exonids = [x for x in exonids if not self.isLastExon(x, df_row)]
            exps = [] 
            for exonid in exonids:
                #seq = self.matrix.iloc[exonid-1, 0:12].values
                seq = self.matrix.iloc[exonid-1, 0:12].values
                if method == 'rolling_time':
                    exp = self.avgseq(seq, self.time_start, self.time_end)
                elif method == 'max':
                    exp = max(seq) 
                elif method == 'min':
                    exp = min(seq)
                exps.append(exp)
            IQs = [max(0, (100-self.getIQ(x))) for x in famids]
            dat = zip(famids, varids, exonids, IQs, exps)
            dat = [x for x in dat if x[2]!= None]
            if len(dat) < 2:
                continue
            self.Genes.append(gene)
            for a in dat:
                var = variant(a[1], a[0], a[3], gene, a[2], a[4])
                self.VarList[var.VarID] = var
                if gene in self.Gene2Var:
                    self.Gene2Var[gene].append(var)
                else:
                    self.Gene2Var[gene] = [var]
    def model2(self, plot=True, var2leave = None):
        self.gene2slope = {}
        ALL_Exp = np.array([])
        ALL_Normed_IQD = np.array([])
        for gene in self.Genes:
            Vars = self.Gene2Var[gene]
            if var2leave != None:
                Vars = [x for x in Vars if x.VarID != var2leave]
            exps = np.array([x.ExonExp for x in Vars])
            IQDs = np.array([x.ProbandIQ for x in Vars])
            regr = self.regGene2(exps, IQDs)
            slope = regr.coef_[0]
            #print gene, slope
            if slope == 0:
                continue
            self.gene2slope[gene] = slope
            ALL_Exp = np.append(ALL_Exp, exps)
            NormIQDs = np.array([x/slope for x in IQDs])
            ALL_Normed_IQD = np.append(ALL_Normed_IQD, NormIQDs)
        #print ALL_Exp, ALL_Normed_IQD
        model = linear_model.LinearRegression(fit_intercept=False)
        ALL_Exp = ALL_Exp.reshape(-1,1)
        R, P = scipy.stats.pearsonr([x[0] for x in ALL_Exp], ALL_Normed_IQD)
        if plot:
            plt.title("R:{}, P:{}, time{}-{}".format(round(R,6),round(P,6), self.time_start, self.time_end))
            model.fit(ALL_Exp, ALL_Normed_IQD)
            plt.scatter(ALL_Exp, ALL_Normed_IQD)
            upper = max(ALL_Normed_IQD)
            plt.plot([0, upper], [0,model.predict(upper)], color='red')
            plt.show()
    def Predict(self, var):
        return self.gene2slope[var.Gene] * var.ExonExp

class variant:
    def __init__(self, VarID, ProbandID, ProbandIQ, Gene, ExonID, ExonExp):
        self.VarID = VarID
        self.ProbandID = ProbandID
        self.ProbandIQ = ProbandIQ
        self.Gene = Gene
        self.ExonID = ExonID
        self.ExonExp = ExonExp
    def show(self):
        print "VarID:{} FamID:{} IQ:{} Gene:{} ExonID:{} ExonExp:{}".format(
        self.VarID, self.ProbandID, self.ProbandIQ, self.Gene, self.ExonID, self.ExonExp)
