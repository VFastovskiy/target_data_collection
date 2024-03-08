from parsing import parsing
from visualization import visualization
from uniprot_mapper_request import uniprot_mapper
import matplotlib.pyplot as plt
import seaborn as sns

import pandas as pd
import numpy as np
import os
import glob
import requests
import xml.etree.ElementTree as ET
import csv
import math
from pathlib import Path
from zipfile import ZipFile
from tempfile import TemporaryDirectory
from rdkit.Chem import PandasTools
from chembl_webresource_client.new_client import new_client
from tqdm.auto import tqdm
import time
import re
import json
import zlib
from urllib.parse import urlparse, parse_qs, urlencode
from requests.adapters import HTTPAdapter, Retry
import matplotlib.ticker as ticker


POLLING_INTERVAL = 3
API_URL = "https://rest.uniprot.org"


retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
session = requests.Session()
session.mount("https://", HTTPAdapter(max_retries=retries))



def print_boxplot(column, name):
    data_to_plot = pd.to_numeric(column, errors='coerce')

    fig, ax = plt.subplots(figsize=(10, 7))
    bp = ax.boxplot(data_to_plot, patch_artist=True, notch=True, vert=False)
    colors = ['#0000FF', '#00FF00', '#FFFF00', '#FF00FF']

    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)

    for whisker in bp['whiskers']:
        whisker.set(color='#8B008B', linewidth=1.5, linestyle=":")

    for cap in bp['caps']:
        cap.set(color='#8B008B', linewidth=2)

    for median in bp['medians']:
        median.set(color='red', linewidth=3)

    for flier in bp['fliers']:
        flier.set(marker='.', color='#e7298a', alpha=0.5)

    ax.set_yticklabels(['pKi', 'data_2', 'data_3', 'data_4'])
    plt.title(f'Box Plot - {name}')

    # Set x-axis ticks and labels (adjust the frequency as needed)
    x_ticks = np.arange(np.floor(data_to_plot.min()), np.ceil(data_to_plot.max()) + 0.2, step=0.5)
    ax.set_xticks(x_ticks)
    ax.set_xticklabels([f'{tick:.1f}' for tick in x_ticks])

    # Remove top and right axes ticks
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()

    plt.savefig(f'../results/step3_data_collection/target_data/{name}_boxplot.svg', bbox_inches='tight')




def print_histogram(column, name):
    x_ticks_step = 0.2
    data_to_plot = pd.to_numeric(column, errors='coerce')

    # Create a histogram or kernel density plot
    plt.figure(figsize=(8, 6))
    sns.distplot(data_to_plot, bins=60, kde=True, hist_kws={'edgecolor': 'black'})

    # Control x ticks
    x_ticks = np.arange(data_to_plot.min(), data_to_plot.max() + 1, x_ticks_step)
    plt.gca().xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: '{:.1f}'.format(x)))
    plt.xticks(x_ticks, fontsize=3)

    plt.title(f'Distribution of pKi')
    plt.xlabel(name)
    plt.ylabel('Density')

    plt.savefig(f'../results/step3_data_collection/target_data/{name}_histogram.svg', bbox_inches='tight')






def get_data_from_binding_db(uniprot_id, affinity_cutoff=None):

    filters_list = ['IC50', 'Ki']

    #1. Make a request to BindingDB
    api_url = "https://bindingdb.org/axis2/services/BDBService/getLigandsByUniprot"
    params = {"uniprot": uniprot_id, "response": "application/xml"}
    if affinity_cutoff is not None:
        params["IC50cutoff"] = affinity_cutoff
    response = requests.get(api_url, params=params)

    vals = {uniprot_id: {'IC50': [], 'Ki': []}}

    if response.status_code == 200:

        #xml_file_path = f'../results/step3_data_collection/bindingdb_xml_response_{uniprot_id}.xml'
        #with open(xml_file_path, 'w', encoding='utf-8') as xml_file:
        #    xml_file.write(response.text)

        for filter in filters_list:
            print(filter)
            data_rows = parse_bindingdb_response(response.text, filter)

            if data_rows:
                header = ["Monomer_ID", "SMILES", "Affinity_Type", "Affinity_Unit", "Affinity_Value"] #"Affinity_Unit"
                output_file_path = f"../results/step3_data_collection/bindingdb_{uniprot_id}_{filter}.csv"
                write_to_csv(data_rows, output_file_path, header)
                print(f"Data written to {output_file_path}")

                input_file_path = output_file_path
                df = pd.read_csv(input_file_path)

                # Convert 'Affinity_Value' column to strings
                # df['Affinity_Value'] = df['Affinity_Value'].astype(str)
                # Convert 'Affinity_Value' column to strings
                #df['Affinity_Value'] = df['Affinity_Value'].astype(str)

                # Filter rows where the 'Affinity_Value' column does not start with '>' or '<'
                #filtered_df = df.loc[~df['Affinity_Value'].str.match('[<>]')]

                #df['Affinity_Value'] = df['Affinity_Value'].astype(float)

                # Use .loc to avoid SettingWithCopyWarning
                #filtered_df = filtered_df.sort_values(by='Affinity_Value', ascending=False)

                #vals[uniprot_id][filter].append(filtered_df['Affinity_Value'].count())
                vals[uniprot_id][filter] = int(df['Affinity_Value'].count())

                # Save the filtered DataFrame to a new CSV file
                output_file_path = f"../results/step3_data_collection/bindingdb_{uniprot_id}_{filter}_clean.csv"
                df.to_csv(output_file_path, index=False)

                #name = f'{uniprot_id}_{filter}'
                #print_boxplot(filtered_df['Affinity_Value'], name)
                #print_histogram(filtered_df['Affinity_Value'], name)

                # Collect SMILES strings for a STD
                #input_file_path = output_file_path
                #smiles_column = filtered_df['SMILES']

                # Save SMILES to a .smi file with each SMILES on a new line
                #with open(f"../results/step3_data_collection/ligands_data_bindingdb_{uniprot_id}_{filter}_clean.smi", "w") as smi_file:
                #    for smiles in smiles_column:
                #        smi_file.write(f"{smiles}\n")

            else:
                print("No data to write.")
        print(f'inside f\n')
        print(vals)
        return vals
    else:
        print(f"Error for UNIPROT ID {uniprot_id}: {response.status_code} - {response.text}")
        return vals






def get_data_from_chembl(uniprot_id, affinity_cutoff=None):
    targets_api = new_client.target
    compounds_api = new_client.molecule
    bioactivities_api = new_client.activity

    targets = targets_api.get(
        target_components__accession=uniprot_id,
        organism="Homo sapiens",
        target_type="SINGLE PROTEIN"
    ).only(
        "target_chembl_id", "organism", "pref_name", "target_type"
    )

    targets_df = pd.DataFrame.from_records(targets)
    print(targets_df)

    if not targets_df.empty:
        target = targets_df.iloc[0]
        chembl_id = target["target_chembl_id"]
        print(f"ChEMBL ID: {chembl_id}")

        bioactivities = bioactivities_api.filter(
            target_chembl_id=chembl_id, type="Ki", relation="=", assay_type="B"
        ).only(
            "activity_id",
            "assay_chembl_id",
            "assay_description",
            "assay_type",
            "molecule_chembl_id",
            "type",
            "standard_units",
            "relation",
            "standard_value",
            "target_chembl_id",
            "target_organism",
        )
        print(f"Length and type of bioactivities object: {len(bioactivities)}, {type(bioactivities)}")

        print(f"Length and type of first element: {len(bioactivities[0])}, {type(bioactivities[0])}")
        print(bioactivities[0])

        bioactivities_df = pd.DataFrame.from_dict(bioactivities)
        print(f"DataFrame shape: {bioactivities_df.shape}")
        print(bioactivities_df["units"].unique())
        print(bioactivities_df.dtypes)
        bioactivities_df = bioactivities_df.astype({"standard_value": "float64"})
        print(bioactivities_df.dtypes)

        bioactivities_df.dropna(axis=0, how="any", inplace=True)
        print(f"DataFrame shape: {bioactivities_df.shape}")

        print(f"Units in downloaded data: {bioactivities_df['standard_units'].unique()}")
        print(
            f"Number of non-nM entries:\
            {bioactivities_df[bioactivities_df['standard_units'] != 'nM'].shape[0]}"
        )

        bioactivities_df = bioactivities_df[bioactivities_df["standard_units"] == "nM"]
        print(f"Units after filtering: {bioactivities_df['standard_units'].unique()}")

        bioactivities_df.drop_duplicates("molecule_chembl_id", keep="first", inplace=True)
        print(f"DataFrame shape: {bioactivities_df.shape}")

        bioactivities_df.reset_index(drop=True, inplace=True)
        print(bioactivities_df)
        #bioactivities_df.head()



    else:
        print("No target information found.")





def parse_bindingdb_response(xml_content, filter):
    root = ET.fromstring(xml_content)

    data_rows = []
    for hit in root.findall('.//bdb:affinities', namespaces={'bdb': 'http://ws.bindingdb.org/xsd'}):
        monomer_id = hit.find('./bdb:monomerid', namespaces={'bdb': 'http://ws.bindingdb.org/xsd'}).text
        smiles = hit.find('./bdb:smiles', namespaces={'bdb': 'http://ws.bindingdb.org/xsd'}).text
        affinity_type = hit.find('./bdb:affinity_type', namespaces={'bdb': 'http://ws.bindingdb.org/xsd'}).text
        #affinity_unit = hit.find('./bdb:affinity_unit', namespaces={'bdb': 'http://ws.bindingdb.org/xsd'}).text
        affinity_value = hit.find('./bdb:affinity', namespaces={'bdb': 'http://ws.bindingdb.org/xsd'}).text.strip()

        # Extracting affinity unit from the affinity_type element
        #affinity_unit_element = hit.find('./bdb:affinity_type', namespaces={'bdb': 'http://ws.bindingdb.org/xsd'})

        affinity_unit = 'nM'

        # Only consider entries with filterVal as affinity type
        if affinity_type == filter:
            data_rows.append([monomer_id, smiles, affinity_type, affinity_unit, affinity_value]) # affinity_unit

    return data_rows






def write_to_csv(data, output_file, header):

    with open(output_file, 'w', newline='') as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(header)
        writer.writerows(data)




def parse_pdb_to_uniprot_json(from_db, to_db, input_file):
    # Converting response to a CSV file PDB_ID : UniProt_ID
    with open(input_file, 'r') as file:
        json_data = json.load(file)

    data = []

    for result in json_data['results']:
        pdb_id = result['from']
        uniprot_id = result['to']['primaryAccession']
        organism = result['to']['organism']['scientificName'] if result['to']['organism']['taxonId'] == 9606 else None
        data.append({'PDB_ID': pdb_id, 'UNIPROT_ID': uniprot_id, 'ORGANISM': organism})

    df = pd.DataFrame(data)

    # 1) Remove Primary Accessions which not Homo Sapiens
    # 2) Keep only experimentally characterized proteins (starts from P)
    df = df[(df['ORGANISM'] == 'Homo sapiens') & (df['UNIPROT_ID'].str.startswith('P'))]

    output_file_path = f'../results/step2_mapping/{from_db}_to_{to_db}_mapping_results_clean.csv'
    df.to_csv(output_file_path, index=False)

    return output_file_path




def parse_uniprot_to_chembl_json(from_bd, to_bd, json_response_input):
    with open(json_response_input, 'r') as file:
        json_data = json.load(file)

    data = [{"uniprot_id": entry["from"], "chembl_id": entry["to"]} for entry in json_data["results"]]
    df = pd.DataFrame(data)
    output_file_path = f'../results/step2_mapping/{from_bd}_to_{to_bd}_mapping_results_clean.csv'
    df.to_csv(output_file_path, index=False)


def print_collected_data_report(collected_data_report):
    with open(collected_data_report, 'r', encoding='utf-8') as json_file:
        data_dictionary = json.load(json_file)

    # Data for some proteins was not absent, and we should change [] to a 0
    for protein_id, values in data_dictionary.items():
        for key in ['IC50', 'Ki']:
            if isinstance(values[key], list) and not values[key]:
                data_dictionary[protein_id][key] = 0

    protein_ids = list(data_dictionary.keys())
    ic50_points = [entry['IC50'] for entry in data_dictionary.values()]
    ki_points = [entry['Ki'] for entry in data_dictionary.values()]

    sorted_indices = sorted(range(len(ki_points)), key=lambda k: ki_points[k], reverse=True)
    protein_ids_sorted = [protein_ids[i] for i in sorted_indices]
    ic50_points_sorted = [ic50_points[i] for i in sorted_indices]
    ki_points_sorted = [ki_points[i] for i in sorted_indices]

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

    # Plotting the bar chart for all data points
    bar_width = 0.1
    index = range(len(protein_ids_sorted))

    bar1 = ax1.bar(index, ic50_points_sorted, bar_width, label='IC50', color='blue', edgecolor='black', linewidth=0.1)
    bar2 = ax1.bar([i + bar_width for i in index], ki_points_sorted, bar_width, label='Ki', color='red',
                   edgecolor='black', linewidth=0.1)

    ax1.set_ylabel('# of Data Points', fontsize=12)
    ax1.set_title('# of Data Points per Protein', fontsize=12)
    ax1.set_xticks([i + bar_width / 2 for i in index])
    ax1.set_xticklabels(protein_ids_sorted, rotation=45, ha='right', fontsize=4)
    ax1.legend(fontsize=5)

    # Add value labels on top of Ki bars
    for bar, values in zip([bar2], [ki_points_sorted]):
        for i, value in enumerate(values):
            ax1.text(bar[i].get_x() + bar[i].get_width() / 2, value + 0.1, str(value), ha='center', va='bottom',
                     fontsize=3)

    # Plotting the bar chart for Ki > 100
    ki_threshold = 100
    ki_indices_gt_100 = [i for i, ki in enumerate(ki_points_sorted) if ki > ki_threshold]
    ki_protein_ids_gt_100 = [protein_ids_sorted[i] for i in ki_indices_gt_100]
    ki_points_gt_100 = [ki_points_sorted[i] for i in ki_indices_gt_100]

    bar3 = ax2.bar(ki_indices_gt_100, ki_points_gt_100, bar_width, label='Ki', color='red', edgecolor='black',
                   linewidth=0.1)

    ax2.set_xlabel('UniProt_ID', fontsize=12)
    ax2.set_ylabel('# of Data Points per Protein (Ki > 100)', fontsize=12)
    ax2.set_xticks([i for i in ki_indices_gt_100])
    ax2.set_xticklabels(ki_protein_ids_gt_100, rotation=45, ha='right', fontsize=4)
    ax2.legend(fontsize=5)

    # Add value labels on top of Ki bars (Ki > 100)
    for bar, values in zip([bar3], [ki_points_gt_100]):
        for i, value in enumerate(values):
            ax2.text(bar[i].get_x() + bar[i].get_width() / 2, value + 0.1, str(value), ha='center', va='bottom',
                     fontsize=3)

    # Set a different scale for the y-axis in the second subplot
    ax2.set_ylim(0, max(ki_points_gt_100) + 400)

    # Adjust x-axis limits for the second subplot
    #ax2.set_xlim(-1, len(ki_indices_gt_100))

    # Adjust layout and save the plot
    plt.tight_layout()
    plt.savefig('../results/step3_data_collection/collected_data_raw_with_subplots.png', dpi=300)


def main():

    # Get the directory of the current script and construct the relative paths
    current_dir = os.path.dirname(__file__)
    pdf_path = os.path.join(current_dir, '../data/data.pdf')
    csv_directory = os.path.join(current_dir, '../results/step1_pdf_parsing/no_mapping_csvs')
    joined_csv = os.path.join(current_dir, '../results/step1_pdf_parsing/no_mapping_csvs/joined/concatenated_scaffolds.csv')
    joined_csv_directory = os.path.join(current_dir, '../results/step1_pdf_parsing/no_mapping_csvs/joined')

    # Flags for controlling
    parsing_flag = False
    vz_flag = False
    mapping_flag = False
    data_collection_flag = True
    conformer_generation = False


    # Step 1. Parsing: pdf -> csv; cvs files creating at ../results/step1_pdf_parsing/no_mapping_csvs
    # Visualisation of csv: unique and common vals for a list of columns
    parsing(parsing_flag, pdf_path)

    if vz_flag:
        visualization(True, True, joined_csv_directory)
        visualization(True, False, csv_directory)


    # Step 2. Mapping: PDB ID -> UniProt ID + BindingDB ID + Chembl ID
    # The goal is to collect binding data from a few resources by corresponding ids
    if mapping_flag:

        # 2.1. PDB_ID -> UniProt_ID
        # Retrieve a JSON formatted response from Uniprot Mapper
        json_response_path = uniprot_mapper('PDB_ID', 'PDB', 'UniProtKB', joined_csv)

        # Converte response to a CSV file (columns: PDB_ID, UniProt_ID, Organism)
        pdb_to_uniprot_mapping = parse_pdb_to_uniprot_json('PDB', 'UniProtKB', json_response_path)


        # 2.2. UniProt_ID -> ChEMBL_ID
        # Retrieve a JSON formatted response from Uniprot Mapper
        json_response_path = uniprot_mapper('UNIPROT_ID', 'UniProtKB_AC-ID', 'ChEMBL', pdb_to_uniprot_mapping)

        # Converte response to a CSV file (columns: UniProt_ID, ChEMBL_ID)
        parse_uniprot_to_chembl_json('UniProtKB_AC-ID', 'ChEMBL', json_response_path)



    # Step 3. Data collection for a chosen target CDK2 (Protein_Name) P24941 (Uniprot_ID)
    # 3.1. bindingDB request
    if data_collection_flag:

        ######### PART 1. Define the target. -> P00918 with 6K Ki. ###################

        #uniprot_id = 'P24941'
        #chembl_id = 'CHEMBL301'
        #mapped_csv = os.path.join(current_dir, '../results/step2_mapping/PDB_to_UniProtKB_mapping_results_clean.csv')


        #df = pd.read_csv(mapped_csv)
        #protein_lst = df['UNIPROT_ID'].unique().tolist()
        #print(protein_lst)
        #print(len(protein_lst)) # 115 proteins

        #protein_lst_slice = protein_lst[0:9]
        #print(protein_lst_slice)


        #affinity_cutoff = None
        #proteins_data_points = {}

        #for protein_id in protein_lst:
        #    print(protein_id + '\n')
        #    data_points = get_data_from_binding_db(protein_id, affinity_cutoff)
        #    proteins_data_points.update(data_points)

        output_file_path = '../results/step3_data_collection/collected_data.json'

        # Writing the dictionary to a JSON file with indentation for alignment
        #with open(output_file_path, 'w', encoding='utf-8') as json_file:
        #    json.dump(proteins_data_points, json_file, indent=4)

        print_collected_data_report(output_file_path)


        #    print(proteins_data_points)
        #get_data_from_chembl(uniprot_id, affinity_cutoff)













if __name__ == "__main__":
    main()














