import statistics
from math import floor
import pandas as pd

shotgun_reference = {
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

amplicon_reference = {
    "Bacteria;Bacteria;Firmicutes;Bacilli;Bacillales;Bacillaceae;Bacillus;Bacillus_subtilis": 17.4,
    "Bacteria;Bacteria;Pseudomonadota;Gammaproteobacteria;Enterobacterales;Enterobacteriaceae;Escherichia;Escherichia-coli": 10.1,
    "Bacteria;Bacteria;Pseudomonadota;Gammaproteobacteria;Enterobacterales;Enterobacteriaceae;Salmonella;Salmonella-enterica": 10.4,
    "Bacteria;Bacteria;Firmicutes;Bacilli;Lactobacillales;Lactobacillaceae;Lactobacillus;Lactobacillus_fermentum": 18.4,
    "Bacteria;Bacteria;Firmicutes;Bacilli;Bacillales;Staphylococcaceae;Staphylococcus;Staphylococcus_aureus": 15.5,
    "Bacteria;Bacteria;Firmicutes;Bacilli;Bacillales;Listeriaceae;Listeria;Listeria_monocytogenes": 14.1,
    "Bacteria;Bacteria;Firmicutes;Bacilli;Lactobacillales;Enterococcaceae;Enterococcus;Enterococcus_faecalis": 9.9,
    "Bacteria;Bacteria;Pseudomonadota;Gammaproteobacteria;Pseudomonadales;Pseudomonadaceae;Pseudomonas;Pseudomonas-aeruginosa": 4.2
}

def calculate_miq(sample_data, reference_dict):
    sample_rel_freq = {taxa: (count / sum(sample_data.values())) * 100 for taxa, count in sample_data.items()}
    sample_percent_of_expected = {taxa: (sample_rel_freq[taxa] / reference_dict[taxa]) * 100 for taxa in sample_data if taxa in reference_dict}

    percent_tolerance_in_standard = 15
    raw_percent_of_expected = list(sample_percent_of_expected.values())
    unadjusted_percent_errors = [100 - value for value in raw_percent_of_expected]
    adjusted_percent_errors_squared = [
        (abs(err) - percent_tolerance_in_standard) ** 2 if abs(err) > percent_tolerance_in_standard else 0
        for err in unadjusted_percent_errors
    ]
    mean_deviation_squared = statistics.mean(adjusted_percent_errors_squared)
    rmse = mean_deviation_squared ** 0.5
    miq_score = 100 - rmse

    return floor(miq_score), rmse

def load_sample_data(file_path):
    df = pd.read_csv(file_path, sep='\t')
    return df

if __name__ == "__main__":
    import sys
    import os

    if len(sys.argv) != 4:
        print("Usage: python script.py <data_type> <sample_data_file> <dna_extraction_kit>")
        sys.exit(1)

    data_type = sys.argv[1].strip().lower()
    sample_data_file = sys.argv[2]
    name = sys.argv[3]

    if not os.path.isfile(sample_data_file):
        print(f"Error: File '{sample_data_file}' not found.")
        sys.exit(1)

    if data_type == 'amplicon':
        reference_dict = amplicon_reference
    elif data_type == 'shotgun':
        reference_dict = shotgun_reference
    else:
        raise ValueError("Invalid data type. Please enter 'amplicon' or 'shotgun'.")

    df = load_sample_data(sample_data_file)
    sample_names = df.columns[1:]  # Exclude the first column which is '#OTU ID'

    print(f"Samples name: {name}\n")

    miq_scores = []
    rmse_values = []

    for sample in sample_names:
        sample_data = df.set_index('#OTU ID')[sample].dropna().to_dict()
        miq_score, rmse = calculate_miq(sample_data, reference_dict)
        miq_scores.append((sample, miq_score, rmse))
        print(f"MIQ score for {sample}: {miq_score}, RMSE: {rmse:.2f}")

    # print("\nSummary of MIQ Scores and RMSE Values:")
    # for sample, miq_score, rmse in miq_scores:
    #     print(f"Sample: {sample}, MIQ score: {miq_score}, RMSE: {rmse:.2f}")
