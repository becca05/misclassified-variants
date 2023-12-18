# Import required libraries
import os
import pandas as pd
import hail as hl


def liftover(csv_file_path, output_file_path):
    hl.init()  # Initialize Hail

    liftover_chain_file = 'grch37_to_grch38.over.chain'

    ht = hl.import_table(csv_file_path, delimiter=',', impute=True)

    grch37 = hl.get_reference('GRCh37')
    grch38 = hl.get_reference('GRCh38')
    if not grch37.has_liftover('GRCh38'):
        grch37.add_liftover(liftover_chain_file, grch38)

    standard_contigs = set(grch37.contigs)
    ht = ht.filter(hl.literal(standard_contigs).contains(ht['CHR']))

    ht = ht.annotate(locus=hl.locus(ht['CHR'], ht['BP'], reference_genome='GRCh37'))
    ht = ht.annotate(new_locus=hl.liftover(ht.locus, 'GRCh38'))

    ht = ht.filter(hl.is_defined(ht.new_locus))
    ht = ht.transmute(CHR=ht.new_locus.contig, BP=ht.new_locus.position)

    ht.export(output_file_path)

    hl.stop()


def cosmic_tsv_to_csv(input_tsv, output_csv):
    """Convert COSMIC data from TSV to CSV."""
    try:
        df = pd.read_csv(input_tsv, sep='\t')
        df.to_csv(output_csv, index=False)
        print(f"File converted successfully and saved to {output_csv}")
    except Exception as e:
        print(f"Error while converting file: {e}")


def add_ref_alt(file):
    """Add reference and alternative alleles to COSMIC data."""
    try:
        df = pd.read_csv(file, encoding='latin1')
        df['ref'] = df['GENOMIC_WT_ALLELE']
        df['alt'] = df['GENOMIC_MUT_ALLELE']
        df[['ref', 'alt']] = df[['ref', 'alt']].fillna('.')
        return df
    except Exception as e:
        print(f"Error processing file {file}: {e}")


def tsv_to_csv(file, file2, name):
    """Convert TSV file to CSV."""
    try:
        df = pd.read_csv(file, sep='\t')
        df.rename(columns={'AF': name}, inplace=True)
        df.to_csv(file2, index=False)
        print('GnomAD is successfully CSV file!')
    except Exception as e:
        print(f"Error processing file {file}: {e}")


def add_genomic_coordinate_to_gnomad(csv_file_path, output_csv_path):
    """Add genomic coordinates to gnomAD data."""
    try:
        # Use tab as the delimiter
        df = pd.read_csv(csv_file_path, delimiter='\t')

        # Function to extract chromosome number
        def extract_chromosome_number(chrom):
            if chrom == 'chrX' or chrom == 'chrY':
                return '23'
            return chrom.replace('chr', '')

        df['MUTATION_GENOME_POSITION'] = df.apply(
            lambda x: f"{extract_chromosome_number(x['CHR'])}:{x['BP']}-{x['BP']}", axis=1)

        # Save the DataFrame as a CSV file
        df.to_csv(output_csv_path, index=False)
        print(f'Genomic coordinate added and saved to {output_csv_path}')
    except Exception as e:
        print(f"Error adding genomic coordinates to gnomAD data: {e}")


def add_genomic_coordinate_to_cosmic(df, output_csv):
    """Add genomic coordinates to COSMIC data."""
    try:
        df['CHROMOSOME'] = df['CHROMOSOME'].apply(lambda x: '23' if x in ['X', 'Y'] else str(x))
        df['GENOME_START'] = df['GENOME_START'].fillna(0).astype(float).astype(int)
        coordinate_column = 'GENOME_STOP' if 'GENOME_STOP' in df.columns else 'GENOME_START'
        if coordinate_column == 'GENOME_STOP':
            df['GENOME_STOP'] = df['GENOME_STOP'].fillna(0).astype(float).astype(int)
        df['MUTATION_GENOME_POSITION'] = df['CHROMOSOME'] + ':' + df['GENOME_START'].astype(str) + '-' + df[
            coordinate_column].astype(str)
        df.to_csv(output_csv, index=False)
        print('Genomic coordinate has been added to COSMIC!')
    except Exception as e:
        print("Error adding genomic coordinates to COSMIC data:", e)


def merge_cosmic_gnomad(cosmic, file_2, file_3):
    """Merge COSMIC and gnomAD data."""
    try:
        # Check if file exists before reading
        if not os.path.exists(file_2):
            print(f"File not found: {file_2}")
            return

        gnomad = pd.read_csv(file_2)
        cosmic_df = pd.read_csv(cosmic)
        output1 = pd.merge(cosmic_df, gnomad, on=['MUTATION_GENOME_POSITION', 'ref', 'alt'], how='inner')
        output1.to_csv(file_3, index=False)
        print(f"{file_2} merged successfully.")
    except Exception as e:
        print(f"Error merging files {file_2}: {e}")


def concat_gnomad(file_paths, output_csv):
    """Concatenate all gnomAD files into one CSV."""
    try:
        dfs = []
        for file_path in file_paths:
            df = pd.read_csv(file_path)
            print(f"Columns in {file_path}: {df.columns}")
            dfs.append(df)

        concat_df = pd.concat(dfs)
        concat_df.to_csv(output_csv, index=False)
        print(f'Concatenated data saved to {output_csv}')
    except Exception as e:
        print(f"Error concatenating files: {e}")



def main(c_file, c_file2, gnomad_files, merged):
    # cosmic_tsv_to_csv(c_file, c_file2)
    # cosmic_df = add_ref_alt(c_file2)
    # add_genomic_coordinate_to_cosmic(cosmic_df, c_file2)
    # Process each gnomAD file
    # Process each gnomAD file
    for tsv_file, (csv_file, group) in gnomad_files.items():
        tsv_to_csv(tsv_file, csv_file, f'AF_{group}')
        liftover(csv_file, f'{group}_lifted.csv')
        add_genomic_coordinate_to_gnomad(f'{group}_lifted.csv', f'{group}_lifted.csv')

    lifted_files = [f"{group}_lifted.csv" for _, group in gnomad_files.values()]

    for lifted_file in lifted_files:
        merged_file = f"{os.path.splitext(lifted_file)[0]}_merged.csv"
        merge_cosmic_gnomad(c_file2, lifted_file, merged_file)
    concat_gnomad(
        [f"{group}_lifted_merged.csv" for _, group in gnomad_files.values()],
        merged)

# Example usage
gnomad_files = {
    '/Users/becca/PycharmProjects/variantmisclassifier/GnomAD-data/gnomad.genomes.r2.1.1.afr.adj.ld_scores-1.ldscore': (
    'gnomad_afr.csv', 'afr'),
    '/Users/becca/PycharmProjects/variantmisclassifier/GnomAD-data/gnomad.genomes.r2.1.1.est.adj.ld_scores-1.ldscore': (
    'gnomad_est.csv', 'est'),
    '/Users/becca/PycharmProjects/variantmisclassifier/GnomAD-data/gnomad.genomes.r2.1.1.amr.adj.ld_scores-1.ldscore': (
    'gnomad_amr.csv', 'amr'),
    '/Users/becca/PycharmProjects/variantmisclassifier/GnomAD-data/gnomad.genomes.r2.1.1.asj.adj.ld_scores-1.ldscore': (
    'gnomad_asj.csv', 'asj'),
    '/Users/becca/PycharmProjects/variantmisclassifier/GnomAD-data/gnomad.genomes.r2.1.1.eas.adj.ld_scores-1.ldscore': (
    'gnomad_eas.csv', 'eas'),
    '/Users/becca/PycharmProjects/variantmisclassifier/GnomAD-data/gnomad.genomes.r2.1.1.fin.adj.ld_scores-1.ldscore': (
    'gnomad_fin.csv', 'fin'),
    '/Users/becca/PycharmProjects/variantmisclassifier/GnomAD-data/gnomad.genomes.r2.1.1.nwe.adj.ld_scores-1.ldscore': (
    'gnomad_nwe.csv', 'nwe'),
    '/Users/becca/PycharmProjects/variantmisclassifier/GnomAD-data/gnomad.genomes.r2.1.1.seu.adj.ld_scores-1.ldscore': (
    'gnomad_seu.csv', 'seu')
}
main(
    '/Users/becca/PycharmProjects/variantmisclassifier/COSMIC-data/Cosmic_MutantCensus_Tsv_v99_GRCh38/Cosmic_MutantCensus_v99_GRCh38.tsv',
    'cosmic_final.csv', gnomad_files, 'final_merged_data.csv')
