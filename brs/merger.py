from Bio import SeqIO
import pandas as pd
import json
import os

def combine_results(directory):

    fasta_sequences = SeqIO.parse(open(f'{directory}/proteins.fasta'),'fasta')
    fasta_df = pd.DataFrame(columns=['ID', 'Sequence'])

    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)
        chrom,start,end = fasta.description.split(' ')[1].split(':')
        fasta_df = fasta_df.append({'ID': name, 'Sequence': sequence, 'chrom': chrom, 'start': start, 'end': end}, ignore_index=True)

    fasta_df = get_abricate_results(directory=directory, fasta_df=fasta_df)
    fasta_df = get_bacmet_results(directory=directory, fasta_df=fasta_df)
    fasta_df = get_platon_results(directory=directory, fasta_df=fasta_df)
    fasta_df = get_ssg_lugia_results(directory=directory, fasta_df=fasta_df)
    fasta_df = get_phispy_results(directory=directory, fasta_df=fasta_df)

    return fasta_df

def get_abricate_results(directory, fasta_df):

    df_abricate = pd.read_csv(f'{directory}/abricate.tab',sep='\t')

    abricate_hits = []

    for _,row in fasta_df.iterrows():
        best_hit = None
        best_ioc = 0
        for _,abricate_row in df_abricate.iterrows():
            if row.chrom.split('.')[0] == abricate_row.SEQUENCE:
                row_interval = set(range(int(row.start.strip('>').strip('<'))-1, int(row.end.strip('>').strip('<'))))
                row_abricate_interval = set(range(abricate_row.START-1, abricate_row.END))
                intersection = len(row_interval.intersection(row_abricate_interval))
                union = len(row_interval.union(row_abricate_interval))
                ioc = intersection / union
                if ioc > best_ioc:
                    best_hit = abricate_row.GENE
                    best_ioc = ioc
        abricate_hits.append(best_hit)
    fasta_df['abricate_genes'] = abricate_hits

    return fasta_df

def get_bacmet_results(directory, fasta_df):

    bacmet_compound = []
    bacmet_gene_name = []
    bacmet_id = []
    bacmet_organism = []
    bacmet_identity = [] 
    bacmet_e_value = [] 

    df_bacmet = pd.read_csv(f'{directory}/bacmet.csv',sep=',')
    df_bacmet['ID'] = df_bacmet['query_id'].map(lambda x: x.split(' ')[0])

    for r, row in fasta_df.iterrows():
        df_query = df_bacmet.query('ID == @row.ID')
        if df_query.shape[0] > 0:
            bacmet_compound.append(df_query.iloc[0].compound)
            bacmet_gene_name.append(df_query.iloc[0].gene_name)
            bacmet_id.append(df_query.iloc[0].bacmet_id)
            bacmet_organism.append(df_query.iloc[0].organism)
            bacmet_identity.append(df_query.iloc[0].identity)
            bacmet_e_value.append(df_query.iloc[0]['e-value'])

        else:
            bacmet_compound.append(None)
            bacmet_gene_name.append(None)
            bacmet_id.append(None)
            bacmet_organism.append(None)
            bacmet_identity.append(None)
            bacmet_e_value.append(None)

    fasta_df['bacmet_compound'] = bacmet_compound
    fasta_df['bacmet_gene'] = bacmet_gene_name
    fasta_df['bacmet_id'] = bacmet_id

    return fasta_df

def get_platon_results(directory, fasta_df):

    #passar pelo platon
    #usar o seqio para carregar as sequencias presentes no genome.plasmid.fasta
    #criar uma lista que se chamará plasmid_ids
    
    if not os.path.isfile(f'{directory}/genome.plasmid.fasta'):
        return fasta_df

    records = SeqIO.parse(f'{directory}/genome.plasmid.fasta', 'fasta')
    plasmid_ids = [record.id for record in records]

    fasta_df['platon_plasmids'] = fasta_df['chrom'].map(lambda chrom_id: chrom_id in plasmid_ids)
    return fasta_df

    #verificar se os nossos genes estão nas ilhas identificadas pelo SSG-LUGIA:

def get_ssg_lugia_results(directory, fasta_df):

    ssg_lugia_islands = []
    ssg_lugia_islands_count = 0

    if not os.path.isfile(f'{directory}/phispy/prophage_coordinates.tsv'):
        fasta_df['genomic_island'] = None
        return fasta_df
    try:
        with open(f'{directory}/ssg_lugia.json') as reader:
            ssg_lugia = json.loads(reader.read())
    except:
        return fasta_df['genomic_island']

    for chrom in ssg_lugia:
        for island in ssg_lugia[chrom]:
            ssg_lugia_islands_count = ssg_lugia_islands_count + 1
            ssg_lugia_islands.append([ssg_lugia_islands_count, chrom, int(island[0]), int(island[1])])

    gene_islands = []
            
    for r, row in fasta_df.iterrows():
        row_interval = set(range(int(row.start.strip('>').strip('<')), int(row.end.strip('>').strip('<'))))
        for island in ssg_lugia_islands:
            island_interval = set(range(island[2], island[3]))
            intersection = row_interval.intersection(island_interval)
            if row.chrom == island[1].split('.')[0] and len(intersection) > 0:
                gene_islands.append(island[0])
                break
        else:
            gene_islands.append(None)

    fasta_df['genomic_island'] = gene_islands

    return fasta_df

def get_phispy_results(directory, fasta_df):

    if not os.path.isfile(f'{directory}/phispy/prophage_coordinates.tsv'):
        fasta_df['phage'] = None
        return fasta_df
    
    try:
        df_phispy = pd.read_csv(f'{directory}/phispy/prophage_coordinates.tsv', sep='\t', header=None)
    except:
        fasta_df['phage'] = None
        return fasta_df

    phispy_islands = []
            
    for r, row in fasta_df.iterrows():
        row_interval = set(range(int(row.start.strip('>').strip('<')), int(row.end.strip('>').strip('<'))+1))
        best_phage = None
        best_phage_iou = 0
        for r_phi, row_phi in df_phispy.iterrows():
            if row_phi[1] == row.chrom:
                row_phi_interval = set(range(int(row_phi[2]), int(row_phi[7])))
                iou = len(row_interval.intersection(row_phi_interval)) / len(row_interval.union(row_phi_interval))
                if iou > best_phage_iou:
                    best_phage     = row_phi[0]
                    best_phage_iou = iou
        
        phispy_islands.append(best_phage)

    fasta_df['phage'] = phispy_islands

    return fasta_df
