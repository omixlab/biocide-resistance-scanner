from Bio import SeqIO

def run_converter_gbk_to_genome_fna(input_gbk_file, output_fna_file):

    records = SeqIO.parse(input_gbk_file, 'genbank')
    SeqIO.write(records, output_fna_file, 'fasta')    

def run_converter_gbk_to_protein_fna (input_gbk_file, output_fna_file):

    output_handle = open(output_fna_file, "w")

    for seq_record in SeqIO.parse(open(input_gbk_file), "genbank"):
        for seq_feature in seq_record.features :
            if seq_feature.type=="CDS" :
                output_handle.write(">%s %s:%s:%s\n%s\n" % (
                    seq_feature.qualifiers['locus_tag'][0],
                    seq_record.id,
                    str(seq_feature.location.start+1),
                    str(seq_feature.location.end),
                    str(seq_feature.location.extract(seq_record.seq).translate())))
                
    output_handle.close()