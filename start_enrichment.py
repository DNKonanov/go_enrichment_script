from goatools.godag_plot import plot_gos, plot_results, plot_goid2goobj
import pandas as pd
from goatools.go_enrichment import GOEnrichmentStudy
from goatools.obo_parser import GODag
import argparse


parser = argparse.ArgumentParser()
parser.add_argument('--inputtype', type=str, default='coordinates', help='target input type ("coordinates" or "targetgenes", default "coordinates")')
parser.add_argument('--genome', type=str, default=None, help='strain name')
parser.add_argument('--genes', type=str, default=None, help='path to *_genes.sif file')
parser.add_argument('--b2g', type=str, default=None, help='path blast2go file')
parser.add_argument('--leftedge', type=int, default=None, help='left coordinate (if "coordinate inputtype is chosen)')
parser.add_argument('--rightedge', type=int, default=None, help='right coordinate (if "coordinate" inputtype is chosen)')
parser.add_argument('--geneslist', type=str, default=None, help='path to file with target genes list (if "taregtgenes" inputtype is chosen, genes names must be divided by end-of-line)')
parser.add_argument('--outname', type=str, default=None, help='output files names prefix')
parser.add_argument('--pvalue', type=float, default=0.05, help='p-value (default 0.05)')
args = parser.parse_args()


def enrich(genome, left_coord, inputtype, right_coord, genes_sif_file, GOs_file, geneslist, outname, pvalue):

    genome = args.genome
    genes_sif_file = args.genes
    GOs_file = args.b2g
    input_genes = pd.read_csv(genes_sif_file, sep='\t')
    input_GOs = open(GOs_file)
    obodag = GODag('gene2go.obo')

    print('Open genes..')
    background = []
    for i in range(len(input_genes.gene)):
        if str(input_genes.genome[i]) != genome:
            continue
        background.append((input_genes.start[i], input_genes.end[i], input_genes.gene[i]))

    if left_coord != None and right_coord != None:

        left_coord = args.leftedge
        right_coord = args.rightedge
        print('Create target..')
        background.sort()
        target = []
        for i in background:
            if i[0] > right_coord:
                break
            if i[1] > left_coord:
                target.append(i[2])
        
    elif geneslist is not None:
        target = [line[:-1] for line in open(geneslist, 'r')]

    background = [i[2] for i in background]
    assoc = {}

    print('Create association table..')
    i = 0
    for line in input_GOs:
        if line.startswith('Tags'):
            continue
        string = [i for i in line[:-1].split('\t') if i != '']
        
        og = string[1]
        background.append(og)
        GO = set([])
        
        for s in string:
            if ':GO:' in s:
                GOs = s.split('; ')
                for el in GOs:
                    GO.add(el[2:])
                
        assoc[og] = GO

    print('Fisher-test..')
    obj = GOEnrichmentStudy(background, assoc, obodag, method=['fdr_bh'], alpha=args.pvalue)
    goea_results_all = obj.run_study(target)
    goea_results_sig = [r for r in goea_results_all if r.get_pvalue() < args.pvalue]
    
    f_out = open(outname + '.txt', 'w')
    for g in goea_results_sig:
        f_out.write(str(g) + '\n')

    plot_results(outname + "_{NS}.png", goea_results_sig)





if args.inputtype == 'coordinates':
    if None in (args.genome, args.genes, args.b2g, args.leftedge, args.rightedge, args.outname):
        print('Invalid input! See -h')

    enrich(args.genome, args.leftedge, args.inputtype, args.rightedge, args.genes, args.b2g, args.geneslist, args.outname, args.pvalue)

elif args.inputtype == 'targetgenes':
    if None in (args.genome, args.genes, args.b2g, args.geneslist, args.outname):
        print('Invalid input! See -h')
    enrich(args.genome, args.leftedge, args.inputtype, args.rightedge, args.genes, args.b2g, args.geneslist, args.outname, args.pvalue)


else:
    print('Invalid input! See -h')

