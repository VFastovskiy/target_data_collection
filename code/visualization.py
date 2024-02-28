import glob
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import plotting




def vz_pdb_id_distribution(csv_paths, sf_concatination_flag):

    for csv in csv_paths:
        df = pd.read_csv(csv)
        df_name = os.path.splitext(os.path.basename(csv))[0]

        # Group by 'Protein_Name' and count the number of unique 'PDB_ID'
        pdb_counts = df.groupby('Protein_Name')['PDB_ID'].nunique().reset_index()

        # Sort the DataFrame by the number of PDB IDs in descending order
        pdb_counts = pdb_counts.sort_values(by='PDB_ID', ascending=False)

        # Select the top 40 proteins
        top_proteins = pdb_counts.head(40)

        # Create a bar plot
        plt.figure(figsize=(12, 6))
        ax = sns.barplot(x='Protein_Name', y='PDB_ID', data=top_proteins, palette='viridis')
        plt.xlabel('Protein Name')
        plt.ylabel('Number of Unique PDB IDs')
        plt.title(f'Distribution of PDB IDs per Protein (Top 40) {df_name}')
        plt.xticks(rotation=45, ha='right')

        # Set x-axis and y-axis ticks with step size 1
        plt.xticks(range(len(top_proteins)), top_proteins['Protein_Name'])
        plt.yticks(range(0, top_proteins['PDB_ID'].max() + 1, 10))

        # Add annotations to each bar with exact number of PDB IDs
        for i, v in enumerate(top_proteins['PDB_ID']):
            ax.text(i, v + 0.1, str(v), ha='center', va='bottom', fontsize=8)

        if sf_concatination_flag:
            output_path = f'../results/step1_pdf_parsing/plots/joined/{df_name}_pdb_id_distribution.svg'
        else:
            output_path = f'../results/step1_pdf_parsing/plots/{df_name}_pdb_id_distribution.svg'

        plt.savefig(output_path)



def csv_to_unique_vals_dict(csv_paths, columns_list):
    df_list = []
    unique_values_dict = {}

    for csv_path in csv_paths:
        df_name = f'{os.path.splitext(os.path.basename(csv_path))[0]}'
        df = pd.read_csv(csv_path)
        df_list.append(df)

        # Get unique values for the specified columns and present it as a {'column_key': num_of_unique_val}
        unique_values = df[columns_list].apply(lambda x: len(x.unique())).to_dict()
        unique_values_dict[df_name] = unique_values

    return df_list, unique_values_dict




def csv_to_common_vals_dict(df_list, columns_list):
    common_unique_values = {}

    for df in df_list:
        for column in df[columns_list].columns:
            unique_values_column = set(df[column].unique())

            if column in common_unique_values:
                common_unique_values[column].intersection_update(unique_values_column)
            else:
                common_unique_values[column] = unique_values_column.copy()

    return common_unique_values








def vz_comparison_of_csvs(csv_paths, columns_list):

    df_list, unique_values_dict = csv_to_unique_vals_dict(csv_paths, columns_list)
    common_unique_values = csv_to_common_vals_dict(df_list, columns_list)

    plot_common_vals(common_unique_values)
    plot_unique_vals(unique_values_dict, columns_list)







def visualization(vz_flag, sf_concatination_flag, csv_directory):
    csv_paths = glob.glob(os.path.join(csv_directory, '*.csv'))
    columns_list = ['Protein_Name', 'PDB_ID', 'Chemical_ID']

    if vz_flag & sf_concatination_flag:
        vz_pdb_id_distribution(csv_paths, True)

    elif vz_flag and not sf_concatination_flag:
        # Plotting reports and writing DataFrames to CSV
        vz_comparison_of_csvs(csv_paths, columns_list)
        #vz_pdb_id_distribution(csv_paths)
        #vz_tree_diagram(csv_paths, columns_list)
        #vz_intersections(csv_paths, columns_list)

    else:
        print('Visualization was skipped or done previously.')