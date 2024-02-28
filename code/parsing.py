import pandas as pd
import os
import tqdm
import tabula as tb



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