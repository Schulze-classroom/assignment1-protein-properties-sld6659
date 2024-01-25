#First, import the following libraries
from Bio import SeqIO #Biopython library that works with biological seq
from pyteomics import mass, electrochem #Pyteomics library for mass calc and pI
import csv #for reading and writing CSV files
import plotly.graph_objects as go #Plotly library for creating interactive plots

#function to calc mw using pyteomics
def calculate_molecular_weight(sequence):
    return mass.calculate_mass(sequence, charge=0, average=False)

#function to calc pI using pyteomics
def calculate_isoelectric_point(sequence):
    return electrochem.pI(sequence)

#function to convert FASTA file to CSV files 
def fasta_to_csv(input_fasta, output_mw_csv, output_pi_csv):
    mw_data = []
    pi_data = []

    with open(input_fasta, 'r') as fasta_file:
        fasta_sequences = SeqIO.to_dict(SeqIO.parse(fasta_file, 'fasta'))

        # Calculate molecular weight and isoelectric point for each sequence in the FASTA file
        for seq_id, seq_record in fasta_sequences.items():
            sequence_str = str(seq_record.seq)
            weight = calculate_molecular_weight(sequence_str)
            isoelectric_point = calculate_isoelectric_point(sequence_str)

            #append data to their lists
            mw_data.append({'Protein_ID': seq_id, 'Molecular_Weight(Da)': weight})
            pi_data.append({'Protein_ID': seq_id, 'Isoelectric_Point': isoelectric_point})

    # Write to molecular weight CSV
    with open(output_mw_csv, 'w', newline='') as mw_csv_file:
        mw_csv_writer = csv.DictWriter(mw_csv_file, fieldnames=['Protein_ID', 'Molecular_Weight(Da)'])
        mw_csv_writer.writeheader()
        mw_csv_writer.writerows(mw_data)

    # Write to isoelectric point CSV
    with open(output_pi_csv, 'w', newline='') as pi_csv_file:
        pi_csv_writer = csv.DictWriter(pi_csv_file, fieldnames=['Protein_ID', 'Isoelectric_Point'])
        pi_csv_writer.writeheader()
        pi_csv_writer.writerows(pi_data)

#function to create a box plot using Plotly and export it as a jpg
def plot_boxplot(csv_file, property_name, export_file):
    data = []

    with open(csv_file, 'r') as csv_file:
        csv_reader = csv.DictReader(csv_file)
        for row in csv_reader:
            data.append(float(row[property_name]))

#creating boxplot
    fig = go.Figure()

    fig.add_trace(
        go.Box(
            y=data,
            name=property_name,
            boxpoints='all',
            jitter=0.3, #shifting points horizontally to avoid overlap
            pointpos=-1.8 #points will be positioned to the left of the boxplot, at a distance of 1.8 time the width of the box
        )
    )

#updates the layout
    fig.update_layout(
        title=f'Boxplot of {property_name}',
        yaxis_title=property_name
    )

    # Export the graph as a JPG file
    fig.write_image(export_file)

#outputs
if __name__ == "__main__":
    input_fasta_file = 'N_meningitidis.fasta' #put the name of your FASTA file here
    mw_csv_file = 'output_mw.csv'
    pi_csv_file = 'output_pi.csv'
    mw_boxplot_file = 'MW_boxplot.jpg'
    pi_boxplot_file = 'pI_boxplot.jpg'

    # Calculate and save molecular weights and isoelectric points to CSVs
    fasta_to_csv(input_fasta_file, mw_csv_file, pi_csv_file)
    print(f"Molecular weights and isoelectric points calculated and saved to CSVs.")

    # Plot a boxplot of molecular weights and export as JPG
    plot_boxplot(mw_csv_file, 'Molecular_Weight(Da)', mw_boxplot_file)
    print(f"Boxplot exported as {mw_boxplot_file}")

    # Plot a boxplot of molecular weights and export as JPG
    plot_boxplot(pi_csv_file, 'Isoelectric_Point', pi_boxplot_file)
    print(f"Boxplot exported as {pi_boxplot_file}")
