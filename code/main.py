from parsing import parsing
from visualization import visualization
from uniprot_mapper_request import uniprot_mapper

import pandas as pd
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


POLLING_INTERVAL = 3
API_URL = "https://rest.uniprot.org"


retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
session = requests.Session()
session.mount("https://", HTTPAdapter(max_retries=retries))








def get_data_from_binding_db(uniprot_id, affinity_cutoff=None):

    filters_list = ['IC50', 'Ki']

    #1. Make a request to BindingDB
    api_url = "https://bindingdb.org/axis2/services/BDBService/getLigandsByUniprot"
    params = {"uniprot": uniprot_id, "response": "application/xml"}
    if affinity_cutoff is not None:
        params["IC50cutoff"] = affinity_cutoff
    response = requests.get(api_url, params=params)


    if response.status_code == 200:

        xml_file_path = f'../results/step3_data_collection/bindingdb_xml_response_{uniprot_id}.xml'
        with open(xml_file_path, 'w', encoding='utf-8') as xml_file:
            xml_file.write(response.text)

        for filter in filters_list:
            data_rows = parse_bindingdb_response(response.text, filter)

            if data_rows:
                header = ["Monomer_ID", "SMILES", "Affinity_Type", "Affinity_Unit", "Affinity_Value"] #"Affinity_Unit"
                output_file_path = f"../results/step3_data_collection/bindingdb_{uniprot_id}_{filter}.csv"
                write_to_csv(data_rows, output_file_path, header)
                print(f"Data written to {output_file_path}")

                input_file_path = output_file_path
                df = pd.read_csv(input_file_path)

                # Filter rows where the last column contains '>' or '<'
                filtered_df = df.loc[~df['Affinity_Value'].str.contains('[<>]')]

                # If 'Ki' values are stored as strings, convert them to numeric for sorting
                filtered_df['Affinity_Value'] = pd.to_numeric(filtered_df['Affinity_Value'], errors='coerce')

                # Use .loc to avoid SettingWithCopyWarning
                filtered_df = filtered_df.sort_values(by='Affinity_Value', ascending=False)

                # Save the filtered DataFrame to a new CSV file
                output_file_path = f"../results/step3_data_collection/bindingdb_{uniprot_id}_{filter}_clean.csv"
                filtered_df.to_csv(output_file_path, index=False)

                # Collect SMILES strings for a STD
                input_file_path = output_file_path
                smiles_column = filtered_df['SMILES']

                # Save SMILES to a .smi file with each SMILES on a new line
                with open(f"../results/step3_data_collection/ligands_data_bindingdb_{uniprot_id}_{filter}_clean.smi", "w") as smi_file:
                    for smiles in smiles_column:
                        smi_file.write(f"{smiles}\n")
            else:
                print("No data to write.")
    else:
        print(f"Error for UNIPROT ID {uniprot_id}: {response.status_code} - {response.text}")






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

    if not targets_df.empty:
        target = targets_df.iloc[0]
        chembl_id = target["target_chembl_id"]
        print(f"ChEMBL ID: {chembl_id}")

        bioactivities = bioactivities_api.filter(
            target_chembl_id=chembl_id, type__in=["IC50", "Ki"], relation="=", assay_type="B"
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
        json_response_path = uniprot_mapper('UNIPROT_ID', 'UniProtKB_AC-ID', 'ChEMBL', pdb_to_uniprot_mapping)

        # Converte response to a CSV file (columns: UniProt_ID, ChEMBL_ID)
        parse_uniprot_to_chembl_json('UniProtKB_AC-ID', 'ChEMBL', json_response_path)



    # Step 3. Data collection for a chosen target CDK2 (Protein_Name) P24941 (Uniprot_ID)
    # 3.1. bindingDB request
    if data_collection_flag:
        uniprot_id = 'P24941'
        chembl_id = 'CHEMBL301'
        affinity_cutoff = None
        get_data_from_binding_db(uniprot_id, affinity_cutoff)
        #get_data_from_chembl(uniprot_id, affinity_cutoff)






if __name__ == "__main__":
    main()














