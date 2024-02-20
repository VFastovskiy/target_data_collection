import tabula as tb
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns
import time
import requests
import json
import xml.etree.ElementTree as ET
import csv


def read_csv_column_unique_vals(csv_path, column_name):

    output_file_path = csv_path.replace('.csv', '.out')

    df = pd.read_csv(csv_path)
    entries = df[column_name].unique().tolist()
    print('len of the list is ', len(entries))

    with open(output_file_path, 'w') as file:
        # Write each item to the file with a space separator
        file.write(' '.join(map(str, entries)))

    # Print a message indicating success
    print(f'The list has been written to {output_file_path}')




def write_to_csv(protein_dict, output_file):
    # Define the CSV file header
    fieldnames = ['UniProt_ID', 'Organism', 'Gene_Name', 'Synonyms', 'Other_Accessions', 'PDB_ID', 'Method', 'Resolution', 'Chains']

    with open(output_file, mode='w', newline='') as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=fieldnames)

        # Write the header
        writer.writeheader()

        for uniprot_id, data in protein_dict.items():
            # Flatten the data for writing to CSV
            for pdb_info in data['pdb_ids']:
                writer.writerow({
                    'UniProt_ID': uniprot_id,
                    'Organism': data['ncbi_taxonomy_id'],
                    'Other_Accessions': ', '.join(data['other_accessions']),
                    'Gene_Name': data['gene_name'],
                    'Synonyms': ', '.join(data['synonyms']),
                    'PDB_ID': pdb_info.get('id', ''),
                    'Method': pdb_info.get('method', ''),
                    'Resolution': pdb_info.get('resolution', ''),
                    'Chains': pdb_info.get('chains', '')
                })


def parse_xml(xml_file):
    # Define the namespace mapping
    ns = {'uniprot': 'http://uniprot.org/uniprot'}

    tree = ET.parse(xml_file)
    root = tree.getroot()

    protein_dict = {}

    for entry in root.findall('.//uniprot:entry', namespaces=ns):
        uniprot_id = entry.find('uniprot:accession', namespaces=ns).text

        # Find all other accessions
        other_accessions = [acc.text for acc in entry.findall('.//uniprot:accession', namespaces=ns) if acc.text != uniprot_id]

        # Find NCBI Taxonomy information
        ncbi_taxonomy_ref = entry.find('.//uniprot:dbReference[@type="NCBI Taxonomy"]', namespaces=ns)
        ncbi_taxonomy_id = ncbi_taxonomy_ref.get('id') if ncbi_taxonomy_ref is not None else None

        # Find all PDB references within the entry
        pdb_references = entry.findall('.//uniprot:dbReference[@type="PDB"]', namespaces=ns)

        # Extract PDB IDs and associated information
        pdb_ids = []
        for pdb_ref in pdb_references:
            pdb_id = pdb_ref.get('id')

            # Extract method, resolution, and chains
            method = pdb_ref.find('.//uniprot:property[@type="method"][@value]', namespaces=ns)
            resolution = pdb_ref.find('.//uniprot:property[@type="resolution"][@value]', namespaces=ns)
            chains = pdb_ref.find('.//uniprot:property[@type="chains"][@value]', namespaces=ns)

            pdb_info = {
                'id': pdb_id,
                'method': method.get('value') if method is not None else None,
                'resolution': resolution.get('value') if resolution is not None else None,
                'chains': chains.get('value') if chains is not None else None
            }

            pdb_ids.append(pdb_info)

        # Extract gene names and synonyms
        gene_name_element = entry.find('.//uniprot:gene/uniprot:name[@type="primary"]', namespaces=ns)
        gene_name = gene_name_element.text if gene_name_element is not None else None

        # Extract synonyms
        synonyms = [name.text for name in entry.findall('.//uniprot:gene/uniprot:name[@type="synonym"]', namespaces=ns)]

        protein_dict[uniprot_id] = {
            'pdb_ids': pdb_ids,
            'gene_name': gene_name,
            'synonyms': synonyms,
            'other_accessions': other_accessions,
            'ncbi_taxonomy_id': ncbi_taxonomy_id
        }

    # Print the results
    for uniprot_id, info in protein_dict.items():
        print(f'Uniprot ID: {uniprot_id}')
        print(f'taxonomy_id: {info["ncbi_taxonomy_id"]}')
        print(f'Other Accessions: {info["other_accessions"]}')
        print(f'Gene Name: {info["gene_name"]}')
        print(f'Synonyms: {info["synonyms"]}')
        print(f'PDB IDs:')
        for pdb_info in info["pdb_ids"]:
            print(
                f'  - ID: {pdb_info["id"]}, Method: {pdb_info["method"]}, Resolution: {pdb_info["resolution"]}, Chains: {pdb_info["chains"]}')
    print('\n')
    print(f'Number of Proteins with proper mapping: {len(protein_dict)}')
    print('\n')

    # Write dictionary into a csv
    # Extract the base name of the input file (excluding path and extension)
    base_name = os.path.splitext(os.path.basename(xml_file))[0]
    output_file_path = f'./data/{base_name}_mapped.csv'

    write_to_csv(protein_dict, output_file_path)





def chembl_mapping(csv_file, tsv_file):

    df_csv = pd.read_csv(csv_file)

    # Load the TSV file into a DataFrame with columns 'uniprot_id' and 'chembl_id'
    df_mapping = pd.read_csv(tsv_file, delimiter='\t', names=['UniProt_ID', 'ChEMBL_ID'])
    df_mapping = df_mapping.drop(0)

    df_merged = pd.merge(df_csv, df_mapping, on='UniProt_ID', how='right')
    fieldnames = ['UniProt_ID', 'ChEMBL_ID', 'PDB_ID', 'Gene_Name', 'Synonyms', 'Other_Accessions', 'Organism', 'Method', 'Resolution', 'Chains']
    df_merged = df_merged[fieldnames]

    # Extract the base name of the input file (excluding path and extension)
    base_name = os.path.splitext(os.path.basename(csv_file))[0]
    output_file_path = f'./data/{base_name}_chemblID_mapped.csv'

    df_merged.to_csv(output_file_path, index=False)






def main():

    current_dir = os.path.dirname(__file__)
    #pyridones_scaffold_csv = os.path.join(current_dir, 'data/pyridones_scaffold.csv')
    #diarylamines_scaffold_csv = os.path.join(current_dir, 'data/diarylamines_scaffold_output.csv')

    pyridones_scaffold_xml = os.path.join(current_dir, '../data/data_mapping/pyridones_pdbID_to_uniprotID.xml')
    diarylamines_scaffold_xml = os.path.join(current_dir, '../data/data_mapping/diarylamines_pdbID_to_uniprotID.xml')

    #pyridones_scaffold_csv = os.path.join(current_dir, 'data/pyridones_uniprotID_mapped.csv')
    #diarylamines_scaffold_csv = os.path.join(current_dir, 'data/diarylamines_uniprotID_mapped.csv')

    pyridones_scaffold_tsv = os.path.join(current_dir, '../data/data_mapping/pyridones_uniprotID_to_chemblID.tsv')
    pyridones_scaffold_csv = os.path.join(current_dir, '../data/pyridones_pdbID_to_uniprotID_mapped.csv')

    diarylamines_scaffold_tsv = os.path.join(current_dir, '../data/data_mapping/diarylamines_uniprotID_to_chemblID.tsv')
    diarylamines_scaffold_csv = os.path.join(current_dir, '../data/diarylamines_pdbID_to_uniprotID_mapped.csv')



    #read_csv(pyridones_scaffold_csv)
    #read_csv(diarylamines_scaffold_csv)


    #reshape_csv(pyridones_scaffold_csv)
    #reshape_csv(diarylamines_scaffold_csv)

    #get_uniprot_id_for_proteins(pyridones_scaffold_csv)

    #result = get_uniprot_id('Abl')
    #print(result)

    #get_uniprot_pdb_mapping(pyridones_scaffold_json, pyridones_scaffold_csv)


    #parse_xml(pyridones_scaffold_xml)
    #parse_xml(diarylamines_scaffold_xml)

    #column_name = 'Uniprot_ID'

    #read_csv_column_unique_vals(pyridones_scaffold_csv, column_name)
    #read_csv_column_unique_vals(diarylamines_scaffold_csv, column_name)

    chembl_mapping(pyridones_scaffold_csv, pyridones_scaffold_tsv)
    chembl_mapping(diarylamines_scaffold_csv, diarylamines_scaffold_tsv)


if __name__ == "__main__":
    main()











