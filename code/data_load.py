from parsing import parsing
from visualization import visualization




import tabula as tb
import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm
import glob
import plotly.express as px
import plotly.graph_objects as go
import plotly.subplots as sp
from plotly.subplots import make_subplots
import plotly
import math
import networkx as nx
from matplotlib import cm
import numpy as np
import matplotlib.gridspec as gridspec
import mplcursors
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



def get_data_from_binding_db(uniprot_id, affinity_cutoff=None):

    filters_list = ['IC50', 'Ki']

    api_url = "https://bindingdb.org/axis2/services/BDBService/getLigandsByUniprot"
    params = {"uniprot": uniprot_id, "response": "application/xml"}
    if affinity_cutoff is not None:
        params["IC50cutoff"] = affinity_cutoff
    response = requests.get(api_url, params=params)


    if response.status_code == 200:
        for filter in filters_list:
            data_rows = parse_bindingdb_response(response.text, filter)

            if data_rows:
                header = ["Monomer_ID", "SMILES", "Affinity_Type", "Affinity_Value"]
                output_file_path = f"../results/step2_data_collection/ligands_data_bindingdb_{uniprot_id}_{filter}.csv"
                write_to_csv(data_rows, output_file_path, header)
                print(f"Data written to {output_file_path}")
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
        affinity_value = hit.find('./bdb:affinity', namespaces={'bdb': 'http://ws.bindingdb.org/xsd'}).text.strip()

        # Only consider entries with filterVal as affinity type
        if affinity_type == filter:
            data_rows.append([monomer_id, smiles, affinity_type, affinity_value])

    return data_rows




def write_to_csv(data, output_file, header):

    with open(output_file, 'w', newline='') as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(header)
        writer.writerows(data)










def main():

    # Get the directory of the current script and construct the relative paths
    current_dir = os.path.dirname(__file__)
    pdf_path = os.path.join(current_dir, '../data/data.pdf')
    csv_directory = os.path.join(current_dir, '../results/step1_pdf_parsing/no_mapping_csvs')
    joined_csv_directory = os.path.join(current_dir, '../results/step1_pdf_parsing/no_mapping_csvs/joined')

    # Flags for controlling
    parsing_flag = False
    vz_flag = True
    mapping_flag = False
    data_collection_flag = False

    # Step 1. Parsing: pdf -> csv; cvs files creating at ../results/step1_pdf_parsing/no_mapping_csvs
    parsing(parsing_flag, pdf_path)

    # Step 2. Visualisation of csv: unique and common vals for a list of columns
    if vz_flag:
        visualization(True, True, joined_csv_directory)
        visualization(True, False, csv_directory)

    # Step 3. Mapping: PDB ID -> UniProt ID + BindingDB ID + Chembl ID
    # The goal is to collect binding data from a few resources

    # 3.1. Protein Name -> UniProt ID

    if mapping_flag:
        csv_path = glob.glob(os.path.join(csv_directory, '*.csv'))
        for csv in csv_path:
            df_name = f'{os.path.splitext(os.path.basename(csv_path))[0]}'
            df = pd.read_csv(csv_path)


    # Step 4. Data collection for a chosen target CDK2 (Protein_Name) P24941 (Uniprot_ID)
    # 4.1. bindingDB request
    if data_collection_flag:
        uniprot_id = 'P24941'
        chembl_id = 'CHEMBL301'
        affinity_cutoff = None
        #get_data_from_binding_db(uniprot_id, affinity_cutoff)
        get_data_from_chembl(uniprot_id, affinity_cutoff)




if __name__ == "__main__":
    main()


