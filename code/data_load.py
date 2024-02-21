import tabula as tb
import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm
import glob
import plotly.express as px
import plotly.graph_objects as go
import plotly
import math



class NamedDataFrame(pd.DataFrame):
    """
    A custom DataFrame with an associated name attribute.

    Parameters:
        *args, **kwargs: Passed to the pandas DataFrame constructor.

    Attributes:
        name (str): A name associated with the DataFrame.
    """

    def __init__(self, *args, **kwargs):
        """
        Initialize a NamedDataFrame.

        Parameters:
            *args, **kwargs: Passed to the pandas DataFrame constructor.
        """
        super().__init__(*args, **kwargs)
        self.name = None




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




def plot_common_vals(common_unique_values):
    # Plotting a bar chart with specified colors
    fig, ax = plt.subplots()
    colors = sns.dark_palette("navy", n_colors=len(common_unique_values))

    # Iterate over the common_unique_values dictionary
    for i, (column_name, values_list) in enumerate(common_unique_values.items()):
        bar_height = len(values_list)

        # Plot each bar
        ax.bar(i, bar_height, color=colors[i], label=column_name)

        # Annotate each bar with the exact number and list of elements (inside the bar)
        annotation_text = f'{bar_height}\n'
        annotation_text += '\n'.join(values_list)

        ax.text(i, bar_height / 2, annotation_text,
                ha='center', va='center', color='white')

    # Set x-axis labels and title
    ax.set_xticks(range(len(common_unique_values)))
    ax.set_xticklabels(common_unique_values.keys())
    #plt.title('Comparison of Common Unique Values')
    #plt.xlabel('Columns')
    plt.ylabel('Number of Common Unique Values')

    plt.savefig('../results/step1_pdf_parsing/plots/common_vals.svg')




def plot_unique_vals(unique_values_dict, columns_list):
    # Convert the dictionary to a DataFrame for plotting
    unique_values_df = pd.DataFrame(unique_values_dict)

    # Plotting a bar chart with specified colors
    colors = sns.dark_palette("navy", n_colors=len(columns_list))
    ax = unique_values_df.plot(kind='bar', rot=0, color=colors)
    #plt.title('Comparison of Unique Values')
    plt.ylabel('Number of Unique Values')

    # Annotate each bar with a number of elements
    for p in ax.patches:
        ax.annotate(f'{int(p.get_height())}', (p.get_x() + p.get_width() / 2., p.get_height()),
                    ha='center', va='center', xytext=(0, 5), textcoords='offset points', clip_on=False)

    plt.savefig('../results/step1_pdf_parsing/plots/unique_vals_with_comparison.svg')





def vz_comparison_of_csvs(csv_paths, columns_list):

    df_list, unique_values_dict = csv_to_unique_vals_dict(csv_paths, columns_list)
    common_unique_values = csv_to_common_vals_dict(df_list, columns_list)

    plot_common_vals(common_unique_values)
    plot_unique_vals(unique_values_dict, columns_list)




def rename_dataframe_columns(dataframe, original_column_names, new_column_names):
    """
    Rename columns in a DataFrame based on specified old and new column names.

    Parameters:
        dataframe (pd.DataFrame): The DataFrame to be modified.
        old_column_names (list): A list of old column names to be replaced.
        new_column_names (list): A list of new column names to replace the old ones.

    Returns:
        pd.DataFrame: A DataFrame with updated column names.
    """
    header_names_list = ['Chemical', 'PDB', 'Protein']

    # Check if at least one element of old_column_names is present in header_names_list
    contains_common_element = any(elem in original_column_names for elem in header_names_list)

    if contains_common_element:
        return dataframe.drop(index=0).rename(columns=dict(zip(original_column_names, new_column_names))).reset_index(drop=True)

    else:
        # Rename columns using new_row_dict
        dataframe = dataframe.rename(columns=dict(zip(original_column_names, new_column_names)))

        # Create a new row with data from renamed row
        new_row = dict(zip(new_column_names, original_column_names))

        # Add a new row at the first place
        return pd.concat([pd.DataFrame(new_row, index=[0]), dataframe], ignore_index=True)




def read_pdf_and_create_dataframe(page_from, page_till, pdf_path):
    """
    Read PDF pages and create a DataFrame by concatenating data from each page.

    Parameters:
        page_from (int): The starting page number.
        page_till (int): The ending page number + 1 (exclusive).
        pdf_path (str): The path to the PDF file.

    Returns:
        pd.DataFrame or None: The concatenated DataFrame or None if no data is found.
    """
    all_dataframes = []

    for page_number in tqdm(range(page_from, page_till), desc='Page parsing', total=page_till-page_from):
        try:
            dataframe_list = tb.read_pdf(pdf_path, pages=str(page_number))
        except Exception as e:
            print(f"Error reading PDF page {page_number}: {e}")
            continue

        if dataframe_list and len(dataframe_list[0].columns) in [4, 5]:
            original_column_names = dataframe_list[0].columns
            common_new_names = ['Chemical_ID', 'PDB_ID', 'Protein_Name']

            if len(dataframe_list[0].columns) == 4:
                new_column_names = common_new_names + ['Torsion_Angle']
            elif len(dataframe_list[0].columns) == 5:
                new_column_names = common_new_names + ['Torsion_Angle_1', 'Torsion_Angle_2']
            else:
                print('Table does not match the expected format')
                continue

            renamed_dataframe = rename_dataframe_columns(dataframe_list[0], original_column_names, new_column_names)
            all_dataframes.append(renamed_dataframe)
        else:
            print(f'No data found in PDF page {page_number}.')

    if all_dataframes:

        # Concatenate all DataFrames vertically
        return pd.concat(all_dataframes, axis=0).reset_index(drop=True)

    else:
        print('No data found in the specified pages.')
        return None




def parse_pdf(pdf_path):

    # Setting page numbers of the source PDF for each DataFrame to be formed: [start_page, end_page+1]
    df_shapes = {'pyridones_scaffold': [3, 6], 'diarylamines_scaffold': [6, 32]}
    df_names = list(df_shapes.keys())
    named_df_list = []

    # Creating named DataFrames of set shapes
    for key, name in tqdm(zip(df_shapes.keys(), df_names), desc="Processing PDFs", total=len(df_shapes)):
        df = read_pdf_and_create_dataframe(df_shapes[key][0], df_shapes[key][1], pdf_path)
        if df is not None and not df.empty:
            named_df = NamedDataFrame(df)
            named_df.name = name
            named_df_list.append(named_df)
        else:
            print(f'Error during reading the .pdf file into a DataFrame: {name}')

    return named_df_list




def save_dataframe_to_csv(dataframe, file_name, directory='../results/step1_pdf_parsing/no_mapping_csvs'):
    """
    Save a DataFrame to a CSV file.

    Parameters:
        dataframe (pd.DataFrame): The DataFrame to be saved.
        file_name (str): The name of the CSV file (without extension).
        directory (str): The directory where the CSV file will be saved.

    Returns:
        bool: True if the operation is successful, False otherwise.
    """
    try:
        # Ensure the directory exists, create it if not
        os.makedirs(directory, exist_ok=True)

        # Construct the full path for the CSV file
        csv_file_path = os.path.join(directory, f'{file_name}.csv')

        # Write the DataFrame to the CSV file
        dataframe.to_csv(csv_file_path, index=False)
        return True
    except Exception as e:
        print(f'Error saving DataFrame to a CSV: {e}')
        return False






def parsing(parsing_flag, pdf_path):

    if parsing_flag:
        named_df_list = parse_pdf(pdf_path)

        if named_df_list:
            for df in named_df_list:
                if not save_dataframe_to_csv(df, str(df.name)):
                    parsing_flag = False
                    print(f'Failed to save DataFrame {df.name}.')
                    break
        else:
            print('Error during reading the .pdf file occurred.')
    else:
        print('Parsing was skipped or done previously')




def vz_tree_diagram(csv_paths, columns_list):
    for csv in csv_paths:
        df = pd.read_csv(csv)
        df_name = os.path.splitext(os.path.basename(csv))[0]

        # Group by 'Protein_Name' and count the number of unique 'PDB_ID'
        pdb_counts = df.groupby('Protein_Name')['PDB_ID'].nunique().reset_index()

        # Sort the DataFrame by the number of PDB IDs in descending order
        pdb_counts = pdb_counts.sort_values(by='PDB_ID', ascending=False)

        # Create a new DataFrame with protein names in sorted order
        protein_names_df = pdb_counts[['Protein_Name']]

        if len(protein_names_df) < 50:
            # Create a single sunburst diagram
            top_proteins = pdb_counts.head(len(protein_names_df))
            fig = px.sunburst(df, path=columns_list)

            # Add annotations for levels of hierarchy
            annotations = []
            rev_column_list = columns_list[::-1]
            for level, column in enumerate(rev_column_list):
                num_unique_values = len(df[column].unique())
                annotation = {
                    'text': f'<b>{column}</b>',
                    'x': 0.5,
                    'y': 1.05 - 0.17 * level,  # Adjust the vertical position
                    'xref': 'paper',
                    'yref': 'paper',
                    'showarrow': False,
                    'font': {'size': 10}
                }
                annotations.append(annotation)

            fig.update_layout(annotations=annotations)

            output_path = f'../results/step1_pdf_parsing/plots/{df_name}_tree_diagram.svg'
            plotly.io.write_image(fig, output_path, format='svg')
        else:
            # Split DataFrame into chunks of 4 protein_names and create sunburst diagram for each chunk
            chunks = [protein_names_df[i:i+4] for i in range(0, len(protein_names_df), 4)]

            for i, chunk_proteins in enumerate(chunks, 1):
                # Plot sunburst diagram for each chunk
                chunk_df = df[df['Protein_Name'].isin(chunk_proteins['Protein_Name'])]
                fig = px.sunburst(chunk_df, path=columns_list)

                # Add annotations for levels of hierarchy
                annotations = []
                rev_column_list = columns_list[::-1]
                for level, column in enumerate(rev_column_list):
                    num_unique_values = len(df[column].unique())
                    annotation = {
                        'text': f'<b>{column}</b>',
                        'x': 0.5,
                        'y': 1.05 - 0.17 * level,  # Adjust the vertical position
                        'xref': 'paper',
                        'yref': 'paper',
                        'showarrow': False,
                        'font': {'size': 10}
                    }
                    annotations.append(annotation)

                fig.update_layout(annotations=annotations)

                output_path = f'../results/step1_pdf_parsing/plots/{df_name}_tree_diagram_part_{i}.svg'
                plotly.io.write_image(fig, output_path, format='svg')




def vz_pdb_id_distribution(csv_paths):

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

        output_path = f'../results/step1_pdf_parsing/plots/{df_name}_pdb_id_distribution.svg'
        plt.savefig(output_path)

def vz_intersections():




def visualization(vz_flag, csv_directory):
    if vz_flag:
        csv_paths = glob.glob(os.path.join(csv_directory, '*.csv'))

        # Plotting reports and writing DataFrames to CSV
        columns_list = ['Protein_Name', 'PDB_ID', 'Chemical_ID']
        #vz_comparison_of_csvs(csv_paths, columns_list)
        #vz_pdb_id_distribution(csv_paths)
        #vz_tree_diagram(csv_paths, columns_list)
        vz_intersections()
    else:
        print('Visualization was skipped of done previously.')




def main():

    # Get the directory of the current script and construct the relative paths
    current_dir = os.path.dirname(__file__)
    pdf_path = os.path.join(current_dir, '../data/data.pdf')
    csv_directory = os.path.join(current_dir, '../results/step1_pdf_parsing/no_mapping_csvs')

    # Flags for controlling
    vz_flag = True
    parsing_flag = False
    mapping_flag = False

    # Step 1. Parsing: pdf -> csv; cvs files creating at ../results/step1_pdf_parsing/no_mapping_csvs
    parsing(parsing_flag, pdf_path)

    # Step 2. Visualisation of csv: unique and common vals for a list of columns
    visualization(vz_flag, csv_directory)

    # Step 3. Mapping: PDB ID -> UniProt ID
    if mapping_flag:
        csv_path = glob.glob(os.path.join(csv_directory, '*.csv'))
        for csv in csv_path:
            df_name = f'{os.path.splitext(os.path.basename(csv_path))[0]}'
            df = pd.read_csv(csv_path)






if __name__ == "__main__":
    main()