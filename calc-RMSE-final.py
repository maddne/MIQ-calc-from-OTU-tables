import statistics
import argparse
from math import floor
import pandas as pd
import matplotlib.pyplot as plt
import csv

shotgun_reference_Z10 = {
    "Bacteria;Bacteria;Firmicutes;Bacilli;Bacillales;Bacillaceae;Bacillus;Bacillus_subtilis": 12,
    "Bacteria;Bacteria;Pseudomonadota;Gammaproteobacteria;Enterobacterales;Enterobacteriaceae;Escherichia;Escherichia-coli": 12,
    "Bacteria;Bacteria;Pseudomonadota;Gammaproteobacteria;Enterobacterales;Enterobacteriaceae;Salmonella;Salmonella-enterica": 12,
    "Bacteria;Bacteria;Firmicutes;Bacilli;Lactobacillales;Lactobacillaceae;Lactobacillus;Lactobacillus_fermentum": 12,
    "Bacteria;Bacteria;Firmicutes;Bacilli;Bacillales;Staphylococcaceae;Staphylococcus;Staphylococcus_aureus": 12,
    "Bacteria;Bacteria;Firmicutes;Bacilli;Bacillales;Listeriaceae;Listeria;Listeria_monocytogenes": 12,
    "Bacteria;Bacteria;Firmicutes;Bacilli;Lactobacillales;Enterococcaceae;Enterococcus;Enterococcus_faecalis": 12,
    "Bacteria;Bacteria;Pseudomonadota;Gammaproteobacteria;Pseudomonadales;Pseudomonadaceae;Pseudomonas;Pseudomonas-aeruginosa": 12,
    "Eukaryota;Fungi;Ascomycota;Saccharomycetes;Saccharomycetales;Saccharomycetaceae;Saccharomyces;Saccharomyces_cerevisiae": 2,
    "Eukaryota;Fungi;Basidiomycota;Tremellomycetes;Tremellales;Cryptococcaceae;Cryptococcus;Cryptococcus_neoformans": 2
}

amplicon_reference_Z21 = { ## lipsvat edni 0.08 - Salmonella enterica 0.01, E. faecalis 0.001, Methanobrevibacter smithii 0.1
    "Bacteria;Verrucomicrobiota;Verrucomicrobiae;Verrucomicrobiales;Akkermansiaceae;Akkermansia;Akkermansia-muciniphila": 0.97,
    "Bacteria;Bacteroidota;Bacteroidia;Bacteroidales;Bacteroidaceae;Bacteroides;Bacteroides-fragilis": 9.94,
    "Bacteria;Actinomycetota;Actinomycetes;Bifidobacteriales;Bifidobacteriaceae;Bifidobacterium;Bifidobacterium-adolescentis": 8.78,
    "Bacteria;Bacillota;Clostridia;Eubacteriales;Peptostreptococcaceae;Clostridioides;Clostridioides-difficile": 2.62,
    "Bacteria;Bacteria;Pseudomonadota;Gammaproteobacteria;Enterobacterales;Enterobacteriaceae;Escherichia;Escherichia-coli": 12.12,
    "Bacteria;Bacillota;Clostridia;Oscillospirales;Ruminococcaceae;Faecalibacterium;Faecalibacterium-prausnitzii": 17.63,
    "Bacteria;Fusobacteriota;Fusobacteriia;Fusobacteriales;Fusobacteriaceae;Fusobacterium;Fusobacterium-nucleatum": 7.49,
    "Bacteria;Bacillota;Bacilli;Lactobacillales;Lactobacillaceae;Lactobacillus;Lactobacillus-fermentum": 9.63,
    "Bacteria;Bacteroidota;Bacteroidia;Bacteroidales;Prevotellaceae;Prevotella;Prevotella-corporis": 4.98,
    "Bacteria;Bacillota;Clostridia;Eubacteriales;Lachnospiraceae;Roseburia;Roseburia-hominis": 9.89,
    "Bacteria;Bacillota;Negativicutes;Veillonellales;Veillonellaceae;Veillonella;Veillonella-rogosae": 15.87,
}

amplicon_reference_Z10 = {
    "Bacteria;Bacteria;Firmicutes;Bacilli;Bacillales;Bacillaceae;Bacillus;Bacillus_subtilis": 17.4,
    "Bacteria;Bacteria;Pseudomonadota;Gammaproteobacteria;Enterobacterales;Enterobacteriaceae;Escherichia;Escherichia-coli": 10.1,
    "Bacteria;Bacteria;Pseudomonadota;Gammaproteobacteria;Enterobacterales;Enterobacteriaceae;Salmonella;Salmonella-enterica": 10.4,
    "Bacteria;Bacteria;Firmicutes;Bacilli;Lactobacillales;Lactobacillaceae;Lactobacillus;Lactobacillus_fermentum": 18.4,
    "Bacteria;Bacteria;Firmicutes;Bacilli;Bacillales;Staphylococcaceae;Staphylococcus;Staphylococcus_aureus": 15.5,
    "Bacteria;Bacteria;Firmicutes;Bacilli;Bacillales;Listeriaceae;Listeria;Listeria_monocytogenes": 14.1,
    "Bacteria;Bacteria;Firmicutes;Bacilli;Lactobacillales;Enterococcaceae;Enterococcus;Enterococcus_faecalis": 9.9,
    "Bacteria;Bacteria;Pseudomonadota;Gammaproteobacteria;Pseudomonadales;Pseudomonadaceae;Pseudomonas;Pseudomonas-aeruginosa": 4.2
}

### There are 5 different E.coli strins included. Have to check whether they can be distinguished with any DB. Try Centrifuge ### 
shotgun_reference_Z21 = {
    "Bacteria;Verrucomicrobiota;Verrucomicrobiae;Verrucomicrobiales;Akkermansiaceae;Akkermansia;Akkermansia-muciniphila": 1.5,
    "Bacteria;Bacteroidota;Bacteroidia;Bacteroidales;Bacteroidaceae;Bacteroides;Bacteroides-fragilis": 14,
    "Bacteria;Actinomycetota;Actinomycetes;Bifidobacteriales;Bifidobacteriaceae;Bifidobacterium;Bifidobacterium-adolescentis": 6,
    "Bacteria;Bacillota;Clostridia;Eubacteriales;Peptostreptococcaceae;Clostridioides;Clostridioides-difficile": 1.5,
    "Bacteria;Bacteria;Pseudomonadota;Gammaproteobacteria;Enterobacterales;Enterobacteriaceae;Escherichia;Escherichia-coli": 14.1,
    "Bacteria;Bacillota;Clostridia;Oscillospirales;Ruminococcaceae;Faecalibacterium;Faecalibacterium-prausnitzii": 14,
    "Bacteria;Fusobacteriota;Fusobacteriia;Fusobacteriales;Fusobacteriaceae;Fusobacterium;Fusobacterium-nucleatum": 6,
    "Bacteria;Bacillota;Bacilli;Lactobacillales;Lactobacillaceae;Lactobacillus;Lactobacillus-fermentum": 6,
    "Bacteria;Bacteroidota;Bacteroidia;Bacteroidales;Prevotellaceae;Prevotella;Prevotella-corporis": 6,
    "Bacteria;Bacillota;Clostridia;Eubacteriales;Lachnospiraceae;Roseburia;Roseburia-hominis": 14,
    "Bacteria;Bacillota;Negativicutes;Veillonellales;Veillonellaceae;Veillonella;Veillonella-rogosae": 14,
    "Fungi;Ascomycota;Saccharomycetes;Saccharomycetales;Saccharomycetaceae;Candida;Candida-albicans": 1.5,
    "Fungi;Ascomycota;Saccharomycetes;Saccharomycetales;Saccharomycetaceae;Saccharomyces;Saccharomyces-cerevisiae": 1.4
}

def calculate_miq(sample_data, reference_dict):
    sample_rel_freq = {taxa: (count / sum(sample_data.values())) * 100 for taxa, count in sample_data.items()}
    sample_percent_of_expected = {taxa: (sample_rel_freq[taxa] / reference_dict[taxa]) * 100 for taxa in sample_data if taxa in reference_dict}

    percent_tolerance_in_standard = 15
    raw_percent_of_expected = list(sample_percent_of_expected.values())
    # print(raw_percent_of_expected)
    unadjusted_percent_errors = [100 - value for value in raw_percent_of_expected]
    adjusted_percent_errors_squared = [
        (abs(err) - percent_tolerance_in_standard) ** 2 if abs(err) > percent_tolerance_in_standard else 0
        for err in unadjusted_percent_errors
    ]
    # print(adjusted_percent_errors_squared)
    mean_deviation_squared = statistics.mean(adjusted_percent_errors_squared)
    rmse = mean_deviation_squared ** 0.5
    miq_score = 100 - rmse

    return floor(miq_score), rmse


def select_mode(type, mock):
    # reference_dict = {}
    if type == 'amplicon':
        if mock == 'zymo-10':
            reference_dict = amplicon_reference_Z10
        else:
            reference_dict = amplicon_reference_Z21
    elif type == 'shotgun':
        if mock == 'zymo-10':
            reference_dict = shotgun_reference_Z10
        else:
            reference_dict = shotgun_reference_Z21     
               
    return reference_dict


def load_sample_data(file_path):
    df = pd.read_csv(file_path, sep='\t')
    return df


def save_results_to_file(results, output_file):
    with open(output_file, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['Sample', 'MIQ Score', 'RMSE'])
        for sample, miq_score, rmse in results:
            writer.writerow([sample, miq_score, rmse])
    print(f"Results saved to {output_file}")
   
    
def plot_miq_scores(results, output_file):
    samples = [result[0] for result in results]
    miq_scores = [result[1] for result in results]

    plt.figure(figsize=(10, 6))
    plt.bar(samples, miq_scores, color='skyblue')
    plt.xlabel('Samples')
    plt.ylabel('MIQ Score')
    plt.title('MIQ Scores for Samples')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(output_file)
    print(f"Figure saved to {output_file}")
    

if __name__ == "__main__":
    import sys
    import os

    # if len(sys.argv) != 4:
    #     print("Usage: python script.py <data_type> <sample_data_file> <dna_extraction_kit>")
    #     sys.exit(1)
    parser=argparse.ArgumentParser(
        description='Calculate the Measurement Integrity Quotient (MIQ) score of a Microbial/Mock Community Standard (MCS) sample based on taxonomic table.',
            epilog="Example: python calc-RMSE-final.py --data_type amplicon --mock --input-file OTU-table.tsv --kit-name Zymo-mini-prep")
    
    parser.add_argument('--data_type', required=True, type=str, choices=['amplicon','shotgun'],
                        help="Specify the data type: 'amplicon' or 'shotgun'.")
    parser.add_argument('--mock', required=True, type=str, choices=['zymo-10', 'zymo-gut-21'],
                        help="Specify the mock variant. Currently, Zymo-D6300 with 10 microorganisms and Zymo-D6331 Gut microbiome standards are available")                    
    parser.add_argument('--input_file', required=True, type=str,
                        help="Path to the sample data tsv file. Must be formatted in advance.")
    parser.add_argument('--output_file', required=True, type=str,
                        help="Path to save the MIQ scores and RMSE values as a CSV file.")
    parser.add_argument('--plot_figure', action='store_true',
                        help="[Bool, Default: False] Generate and save a bar plot of MIQ scores.")
    
    args=parser.parse_args()
    
    data_type = args.data_type.lower()
    mock = args.mock.lower()
    reference_dict = select_mode(data_type, mock)
    df = load_sample_data(args.input_file)
    sample_names = df.columns[1:]  # Exclude the first column which is '#OTU ID'
    miq_scores = []
   
    for sample in sample_names:
        sample_data = df.set_index('#OTU ID')[sample].dropna().to_dict()
        miq_score, rmse = calculate_miq(sample_data, reference_dict)
        miq_scores.append((sample, miq_score, rmse))
        print(f"MIQ score for sample '{sample}': {miq_score}, RMSE: {rmse:.2f}")
    
     # Save results to a file
    save_results_to_file(miq_scores, args.output_file)
    
    # Generate and save a figure if requested
    if args.plot_figure:
        plot_output_file = args.output_file.replace('.csv', '.png')
        plot_miq_scores(miq_scores, plot_output_file)
