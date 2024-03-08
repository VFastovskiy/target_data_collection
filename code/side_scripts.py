# Assuming you have a column named 'column_to_compare' in each DataFrame
unique_values = set(named_df[column_to_compare].unique())
named_df_dict[f'unique_values_{df_name}'] = unique_values

for df_number in range(len(named_df_list)):
    df_name = f'{named_df_list[df_number].name}'
    named_df_dict[df_name] = named_df_list[df_number]

# Get unique values for each CSV as lists
unique_values1 = set(df1[column_to_compare].unique())
unique_values2 = set(df2[column_to_compare].unique())

# Find common values
common_values = unique_values1.intersection(unique_values2)

# Choose colors from Seaborn's dark_palette
colors = sns.dark_palette("navy", n_colors=3)

# Calculate the number of unique proteins
num_unique1 = len(unique_values1)
num_unique2 = len(unique_values2)
num_common = len(common_values)

# Create a bar plot
plt.bar(['Pyridones', 'Diarylamines', 'Common'], [num_unique1, num_unique2, num_common], color=colors)

# Add exact number of unique proteins on top of each bar
for i, val in enumerate([num_unique1, num_unique2]):
    plt.text(i, val + 0.1, str(val), ha='center', va='bottom', fontsize=10)

# Add annotations for common proteins
common_proteins_annotation = "\n".join(common_values)
plt.annotate(common_proteins_annotation, xy=(2, num_common), ha='center', va='bottom', fontsize=10)

# Add labels and title
plt.xlabel('CSV Files')
plt.ylabel('Number of Proteins')
plt.title(f'Difference in Unique Values for {column_to_compare}')


































def visualize_dataframe(df, columns_list):
    """
    Visualize a DataFrame by creating a bar plot showing the number of unique values in each column.

    Parameters:
        df (pd.DataFrame): The DataFrame to be visualized.
        columns_list (list): A list of column names to be included in the plot.

    Returns:
        None
    """
    subset_df = df[columns_list]

    # Create a bar plot showing the number of unique values in each column
    plt.figure(figsize=(12, 6))
    sns.set(style="whitegrid")
    sns.set_palette("dark")

    sns.barplot(x=subset_df.nunique(), y=subset_df.columns)

    # Add plot labels and title
    plt.title(f'Number of Unique Proteins in CSV ({df.name})')
    plt.xlabel('Number of Unique Entries')
    plt.ylabel('Columns')

    # Save the plot as a PNG file
    plot_filename = f'unique_values_plot_{df.name}.png'
    plt.savefig(plot_filename, bbox_inches='tight')
    print(f'Visualization saved to {plot_filename}')
































    common_unique_values = {}

    for csv_path in csv_paths:
        df_name = f'{os.path.splitext(os.path.basename(csv_path))[0]}'
        df = pd.read_csv(csv_path)
        for column in df.columns:
            unique_values_df1 = set(df1[column].unique())
            unique_values_df2 = set(df2[column].unique())
            common_values = unique_values_df1.intersection(unique_values_df2)
            common_unique_values[column] = list(common_values)




    i = 0
    j = 0
    keys_list = list(unique_values_dict.keys())

    while j < len(keys_list) and i <= (len(unique_values_dict) - 1):
        for column in unique_values_dict[keys_list[j]]:
            vals_df_i = (unique_values_dict[keys_list[j]])[column]
            vals_df_ii = (unique_values_dict[keys_list[j+1]])[column]
            print(vals_df_i, vals_df_ii)
        i+=1
        j+=1





    # Convert the dictionary to a DataFrame for plotting
    unique_values_df = pd.DataFrame(unique_values_dict)
    common_unique_values_df = pd.DataFrame(common_unique_values)

    # Plotting a bar chart with specified colors
    colors = sns.dark_palette("navy", n_colors=len(columns_list))
    ax = unique_values_df.plot(kind='bar', rot=0, color=colors)
    plt.title('Comparison of Unique Values')
    plt.xlabel('Columns')
    plt.ylabel('Number of Unique Values')

    # Annotate each bar with a number of elements
    for p in ax.patches:
        ax.annotate(f'{int(p.get_height())}', (p.get_x() + p.get_width() / 2., p.get_height()),
                    ha='center', va='center', xytext=(0, 10), textcoords='offset points')





    # Convert the dictionary to a DataFrame for plotting
    unique_values_df = pd.DataFrame(unique_values_dict)

    # Plotting a bar chart with specified colors
    colors = sns.dark_palette("navy", n_colors=len(columns_list))
    ax = unique_values_df.plot(kind='bar', rot=0, color=colors)
    plt.title('Comparison of Unique Values')
    plt.xlabel('Columns')
    plt.ylabel('Number of Unique Values')

    # Annotate each bar with a number of elements
    for p in ax.patches:
        ax.annotate(f'{int(p.get_height())}', (p.get_x() + p.get_width() / 2., p.get_height()),
                    ha='center', va='center', xytext=(0, 10), textcoords='offset points')







def vz_tree_diagram(csv_paths, columns_list):
    for csv in csv_paths:
        df = pd.read_csv(csv)
        protein_names = df['Protein_Name'].unique().tolist()
        print(len(protein_names))

        if len(protein_names) < 50:
            # Create a single sunburst diagram
            fig = px.sunburst(df, path=columns_list)
            df_name = os.path.splitext(os.path.basename(csv))[0]
            output_path = f'../results/step1_pdf_parsing/plots/{df_name}_tree_diagram.svg'
            plotly.io.write_image(fig, output_path, format='svg')
        else:
            # Split DataFrame by protein names and create sunburst diagram for each chunk
            protein_name_chunks = split_dataframe_by_protein_names(df, chunk_size=50)
            for i, chunk in enumerate(protein_name_chunks, 1):
                fig = px.sunburst(chunk, path=columns_list)
                df_name = os.path.splitext(os.path.basename(csv))[0]
                output_path = f'../results/step1_pdf_parsing/plots/{df_name}_tree_diagram_part_{i}.svg'
                plotly.io.write_image(fig, output_path, format='svg')




def vz_graph(df, columns_list, df_name):
    # Assuming your filtered_df has columns: 'Protein_Name', 'PDB_ID', 'Chemical_ID'
    G = nx.Graph()

    # Add nodes and edges to the graph
    for _, row in df.iterrows():
        protein_name = row[columns_list[0]]
        pdb_id = row[columns_list[1]]
        chemical_id = row[columns_list[2]]

        G.add_node(protein_name, level=1)  # Assign level 1 to protein_name
        G.add_node(pdb_id, level=2)        # Assign level 2 to pdb_id
        G.add_edge(protein_name, pdb_id)

        G.add_node(chemical_id, level=3)   # Assign level 3 to chemical_id
        G.add_edge(pdb_id, chemical_id)

    # STEP 1
    # Draw the graph with different colors for each level
    levels = nx.get_node_attributes(G, 'level')
    colors = [levels[node] for node in G.nodes]

    # Create a list of edge colors based on the levels of the connected nodes
    edge_colors = [levels[u] for u, v in G.edges]

    # Adjust the scale parameter to control the length of edges
    pos = nx.fruchterman_reingold_layout(G, scale=10.0, k=0.11, seed=127)

    nx.draw(G, pos, with_labels=True, font_size=4, node_color=colors, cmap=plt.cm.spring,
            node_size=60, font_color='black', font_weight='bold', width=0.5, edge_color=edge_colors)

    # Save the graph as an SVG file
    output_path = f'../results/step1_pdf_parsing/plots/{df_name}_common_vals_graph.svg'
    plt.savefig(output_path, format='svg')










def vz_graph(df, columns_list, df_name, common_values):
    # Assuming your filtered_df has columns: 'Protein_Name', 'PDB_ID', 'Chemical_ID'
    G = nx.Graph()

    # Add nodes and edges to the graph
    for _, row in df.iterrows():
        protein_name = row[columns_list[0]]
        pdb_id = row[columns_list[1]]
        chemical_id = row[columns_list[2]]

        G.add_node(protein_name, level=1)  # Assign level 1 to protein_name
        G.add_node(pdb_id, level=2)        # Assign level 2 to pdb_id
        G.add_edge(protein_name, pdb_id)

        G.add_node(chemical_id, level=3)   # Assign level 3 to chemical_id
        G.add_edge(pdb_id, chemical_id)

    # STEP 1
    # Draw the graph with different colors for each level
    levels = nx.get_node_attributes(G, 'level')
    colors = [levels[node] for node in G.nodes]

    # Create a list of node colors based on common_values
    node_colors = ['red' if node in common_values else 'grey' for node in G.nodes]

    # Create a list of edge colors based on the levels of the connected nodes
    edge_colors = [levels[u] for u, v in G.edges]

    # Adjust the scale parameter to control the length of edges
    pos = nx.fruchterman_reingold_layout(G, scale=8.0, k=0.10, seed=127)

    nx.draw(G, pos, with_labels=False, font_size=3, node_color=node_colors, cmap=plt.cm.spring,
            node_size=10, font_color='black', font_weight='bold', width=1, edge_color=edge_colors)

    # Draw node labels above the nodes
    labels = {node: node for node in G.nodes}
    nx.draw_networkx_labels(G, pos, labels, font_size=4, font_color='black', font_weight='bold', verticalalignment='bottom')

    # Save the graph as an SVG file
    output_path = f'../results/step1_pdf_parsing/plots/{df_name}_common_vals_graph.svg'
    plt.savefig(output_path, format='svg')




    # STEP 2. Making a list of Subgraphs
    # Identify connected components (subgraphs)
    subgraphs = list(nx.connected_components(G))

    # Draw each connected component separately
    for i, subgraph_nodes in enumerate(subgraphs):
        subgraph = G.subgraph(subgraph_nodes)

        # Assign distinct colors to edges based on their levels
        edge_colors = [levels[u] for u, v in subgraph.edges]

        # Create a list of node colors based on common_values
        node_colors = ['red' if node in common_values else 'grey' for node in G.nodes]

        # Save each subgraph as an SVG file
        output_path = f'../results/step1_pdf_parsing/plots/{df_name}_subgraph_{i + 1}.svg'
        pos_subgraph = nx.fruchterman_reingold_layout(subgraph, scale=5.0, k=0.10, seed=127)
        nx.draw(subgraph, pos_subgraph, with_labels=True, font_size=4, node_size=20, font_color='black',
                font_weight='bold', width=0.5, edge_color=edge_colors)
        plt.savefig(output_path, format='svg')












   # STEP 2
    # Identify connected components (subgraphs)
    # STEP 2
    # Identify connected components (subgraphs)
    subgraphs = list(nx.connected_components(G))

    # Draw each connected component separately
    for i, subgraph_nodes in enumerate(subgraphs):
        subgraph = G.subgraph(subgraph_nodes)

        levels = nx.get_node_attributes(subgraph, 'level')
        colors = [levels[node] for node in subgraph.nodes]

        # Create a list of node colors based on common_values
        node_colors = ['red' if node in common_values else 'grey' for node in subgraph.nodes]

        # Assign distinct colors to edges based on their levels
        edge_colors = [levels[u] for u, v in subgraph.edges]

        # Save each subgraph as an SVG file
        output_path = f'../results/step1_pdf_parsing/plots/{df_name}_subgraph_{i + 1}.svg'
        plt.figure()

        # Adjust the scale parameter to control the length of edges
        pos = nx.fruchterman_reingold_layout(subgraph, scale=0.5, k=0.05, seed=131)

        nx.draw(subgraph, pos, with_labels=False, font_size=3, node_color=node_colors, cmap=plt.cm.spring,
                node_size=10, font_color='black', font_weight='bold', width=1, edge_color=edge_colors)

        # Draw node labels above the nodes
        labels = {node: node for node in subgraph.nodes}
        nx.draw_networkx_labels(subgraph, pos, labels, font_size=3, font_color='black', font_weight='bold',
                                verticalalignment='bottom')

        plt.savefig(output_path, format='svg')












        #print(proteins_data_points)

        proteins_data_points = {'P00519': {'IC50': 2181, 'Ki': 641},
                                'P00352': {'IC50': 224, 'Ki': 11},
                                'P48730': {'IC50': 1529, 'Ki': 3},
                                'P00533': {'IC50': 9192, 'Ki': 209},
                                'P08246': {'IC50': 1718, 'Ki': 776},
                                'P29320': {'IC50': 9, 'Ki': 0},
                                'P00742': {'IC50': 1821, 'Ki': 2556},
                                'P04035': {'IC50': 168, 'Ki': 9},
                                'P01116': {'IC50': 3117, 'Ki': 6}}













        # 2.3. Handle duplicates: mean val for monomer_id with more than 1 entry
        grouped_df = filtered_df.groupby('Monomer_ID')['pKi'].mean().reset_index()
        clean_df = pd.merge(filtered_df, grouped_df, on='Monomer_ID', how='inner', suffixes=('', '_mean'))
        clean_df.drop('pKi', axis=1, inplace=True)
        clean_df.rename(columns={'pKi_mean': 'pKi'}, inplace=True)
        clean_df = clean_df.drop_duplicates(subset=['Monomer_ID'], ignore_index=True)

        print("Duplicated Rows:")
        print(filtered_df[filtered_df.duplicated(subset=['Monomer_ID'])])
        num_duplicates = len(filtered_df) - len(clean_df)
        print(f"Number of duplicate rows: {num_duplicates}")
        num_dropped_rows = len(filtered_df) - len(clean_df)
        percentage_dropped_out = (num_dropped_rows / total_rows) * 100
        print(
            f"Percentage of filtered out rows: {percentage_dropped_out:.2f}% which is {num_dropped_rows} rows in total")


        output_file_path = '../results/step3_data_collection/target_data/bindingdb_P00918_Ki_dropped_dupl.csv'
        clean_df.to_csv(output_file_path)







        # 2.4. Save SMILES to a .smi file with each SMILES on a new line for Standardization
        smiles_column = clean_df['SMILES']

        with open("../results/step3_data_collection/target_data/bindingdb_P00918_clean.smi", "w") as smi_file:
            for i, smiles in enumerate(smiles_column):
                smi_file.write(smiles)
                if i < len(smiles_column) - 1:
                    smi_file.write("\n")

        output_file_path = '../results/step3_data_collection/target_data/bindingdb_P00918_Ki_clean.csv'
        clean_df.to_csv(output_file_path, index=False)


        # 2.4. Standardization of smiles
        # indigo_std_v2_static_Ubuntu22
        # input: bindingdb_P00918_clean.smi
        # output: bindingdb_P00918_clean_std.smi

        # 2.5 filter out non-valid smiles

        input_smiles_file_path = "../results/step3_data_collection/target_data/bindingdb_P00918_clean_std.smi"
        output_csv_file_path = "../results/step3_data_collection/target_data/bindingdb_P00918_clean_std.csv"


        # 2.5.1. Read SMILES file and assign internal IDs to drop non-valid smiles on 2.5.3.
        smiles_list = []
        to_drop_ids = []
        with open(input_smiles_file_path, "r") as smi_file:
            for i, line in enumerate(smi_file, start=1):
                smiles = line.strip()
                if smiles == 'O':
                    to_drop_ids.append(int(i))
                smiles_list.append({"Internal_ID": i, "SMILES": smiles})

        invalid_smiles = "../results/step3_data_collection/target_data/bindingdb_P00918_invalid_smiles_ids.out"

        with open(invalid_smiles, 'w') as file:
            for value in to_drop_ids:
                file.write(f'{value}\n')


        # 2.5.2. Write data to CSV file
        fieldnames = ["Internal_ID", "SMILES"]
        with open(output_csv_file_path, "w", newline="", encoding="utf-8") as csv_file:
            writer = csv.DictWriter(csv_file, fieldnames=fieldnames)

            # Write header
            writer.writeheader()

            # Write rows
            writer.writerows(smiles_list)


        # 2.5.3. filter out non-standard smiles
        input_csv_file_path = "../results/step3_data_collection/target_data/bindingdb_P00918_Ki_clean.csv"
        output_csv_file_path = "../results/step3_data_collection/target_data/bindingdb_P00918_Ki_clean_std_final.csv"

        df = pd.read_csv(input_csv_file_path)

        # Remove rows with specified IDs
        df = df[~df['Internal_ID'].isin(to_drop_ids)]
        df.to_csv(output_csv_file_path)










def print_boxplot(column, name):

    # Convert the column of strings to numeric values
    data_to_plot = pd.to_numeric(column, errors='coerce')
    print(data_to_plot)

    # Create a horizontal box plot with seaborn for better styling
    plt.figure(figsize=(8, 6))  # Adjust the figure size as needed
    sns.boxplot(y=data_to_plot, orient='h', width=0.5, color='skyblue', fliersize=5)

    plt.title(f'Box Plot {name}')
    plt.xlabel(name)

    plt.savefig(f'../results/step3_data_collection/target_data/{name}_horizontal_boxplot.svg',
                bbox_inches='tight')










   ######### PART 2. Clean and Report the target data ###################

        report_dict = {'P00918': {}}



        # 2.1. Filter out >< values
        input_file_path = '../results/step3_data_collection/target_data/bindingdb_P00918_Ki.csv'
        df = pd.read_csv(input_file_path)

        total_rows = len(df)
        condition = ~df['Affinity_Value'].str.match('[<>]')
        filtered_df = df[condition]
        filtered_rows = len(filtered_df)
        filtered_out_rows = total_rows - filtered_rows
        percentage_filtered_out = (filtered_out_rows / total_rows) * 100
        rounded = round(percentage_filtered_out, 3)
        report_dict['P00918']['inaccurate_values_num'] = filtered_out_rows
        report_dict['P00918']['inaccurate_values_perc'] = rounded
        print(f"Percentage of filtered out rows: {percentage_filtered_out:.2f}% which is {filtered_out_rows} rows in total.")


        # 2.2. Ki -> pKi: nM->M then np.log10 then round vals
        filtered_df = filtered_df.astype({"Affinity_Value": "float64"})
        filtered_df['Affinity_Value'] = (-np.log10(filtered_df['Affinity_Value'] * 1e-9)).round(3)
        filtered_df.drop(['Affinity_Unit', 'Affinity_Type'], axis=1, inplace=True)
        filtered_df.rename(columns={'Affinity_Value': 'pKi'}, inplace=True)


        # 2.3. Handle duplicates: mean val for monomer_id with more than 1 entry
        grouped_df = filtered_df.groupby('Monomer_ID')['pKi'].mean().round(3).reset_index()
        clean_df = pd.merge(filtered_df, grouped_df, on='Monomer_ID', how='inner', suffixes=('', '_mean'))
        clean_df.drop('pKi', axis=1, inplace=True)
        clean_df.rename(columns={'pKi_mean': 'pKi'}, inplace=True)
        clean_df = clean_df.drop_duplicates(subset=['Monomer_ID'], ignore_index=True)
        clean_df.sort_values(by='pKi', ascending=False, inplace=True)

        # Add internal_id to delete non-valid rows on step 2.5.3.
        clean_df['Internal_ID'] = list(range(1, len(clean_df) + 1))
        column_order = ['Monomer_ID', 'Internal_ID', 'SMILES', 'pKi']
        clean_df = clean_df[column_order]

        print("Duplicated Rows:")
        print(filtered_df[filtered_df.duplicated(subset=['Monomer_ID'])])
        num_duplicates = len(filtered_df) - len(clean_df)
        print(f"Number of duplicated rows: {num_duplicates}")
        num_dropped_rows = len(filtered_df) - len(clean_df)
        percentage_dropped_out = (num_dropped_rows / total_rows) * 100
        rounded = round(percentage_dropped_out, 3)
        report_dict['P00918']['duplicated_rows_num'] = num_dropped_rows
        report_dict['P00918']['duplicated_rows_perc'] = rounded
        print(
            f"Percentage of filtered out rows: {percentage_dropped_out:.2f}% which is {num_dropped_rows} rows in total")

        output_file_path = '../results/step3_data_collection/target_data/bindingdb_P00918_Ki_dropped_dupl.csv'
        clean_df.to_csv(output_file_path, index=False)


        # 2.4. Save SMILES to a .smi file with each SMILES on a new line for Standardization
        smiles_column = clean_df['SMILES']
        smiles_to_std = '../results/step3_data_collection/target_data/bindingdb_P00918_clean.smi'

        with open(smiles_to_std, "w") as smi_file:
            for i, smiles in enumerate(smiles_column):
                smi_file.write(smiles)
                if i < len(smiles_column) - 1:
                    smi_file.write("\n")



        # 2.5. Standardization of smiles
        # indigo_std_v2_static_Ubuntu22
        # input: bindingdb_P00918_clean.smi
        # output: bindingdb_P00918_clean_STDRTSED.smi


        # 2.6 filter out non-valid smiles

        # 2.6.1. Read SMILES file and assign internal IDs to drop non-valid smiles on 2.6.2.
        input_smiles_file_path = "../results/step3_data_collection/target_data/bindingdb_P00918_clean_STDRTSED.smi"

        smiles_list = []
        to_drop_ids = []
        with open(input_smiles_file_path, "r") as smi_file:
            for i, line in enumerate(smi_file, start=1):
                smiles = line.strip()
                if smiles == 'O':
                    to_drop_ids.append(int(i))
                smiles_list.append({"Internal_ID": i, "SMILES": smiles})

        invalid_smiles_num = len(to_drop_ids)
        invalid_smiles_perc = ((invalid_smiles_num / total_rows) * 100)
        rounded_invalid_smiles_perc = round(invalid_smiles_perc, 3)

        report_dict['P00918']['invalid_smiles_num'] = invalid_smiles_num
        report_dict['P00918']['invalid_smiles_perc'] = rounded_invalid_smiles_perc



        # Collect invalid ids as a part of data preparation report
        invalid_smiles = "../results/step3_data_collection/target_data/bindingdb_P00918_invalid_smiles_ids.out"

        with open(invalid_smiles, 'w') as file:
            for value in to_drop_ids:
                file.write(f'{value}\n')


        # 2.6.2. drop non-valid smiles using invalid_smiles_id from to_drop_ids
        output_file_path = '../results/step3_data_collection/target_data/bindingdb_P00918_Ki_dropped_dupl_and_inalid_smiles.csv'

        clean_df = clean_df[~clean_df['Internal_ID'].isin(to_drop_ids)]
        clean_df['SMILES'] = clean_df['SMILES'].astype(str)
        clean_df.to_csv(output_file_path, index=False)


        # 2.6.3. Target Data Report
        dtypes = clean_df.dtypes
        report_dict['P00918']['dtypes'] = dtypes.to_dict()
        dtypes_str = {key: str(value) for key, value in dtypes.items()}
        report_dict['P00918']['dtypes'] = dtypes_str

        pKi_statistics = clean_df['pKi'].describe()
        report_dict['P00918']['pKi_statistics'] = {
            'mean': (pKi_statistics['mean']).round(3),
            'std': (pKi_statistics['std']).round(3),
            '25%': pKi_statistics['25%'],
            '50%': pKi_statistics['50%'],
            '75%': pKi_statistics['75%'],
            'max': pKi_statistics['max']
        }

        loss = filtered_out_rows + num_dropped_rows + invalid_smiles_num
        loss_perc = ((filtered_out_rows + num_dropped_rows + invalid_smiles_num) / total_rows) * 100
        loss_perc_round = round(loss_perc, 3)

        report_dict['P00918']['dropped_points_num'] = loss
        report_dict['P00918']['dropped_points_perc'] = loss_perc_round
        report_dict['P00918']['raw_amount_of_points'] = total_rows
        report_dict['P00918']['resulting_amount_of_points'] = len(clean_df)

        output_path = '../results/step3_data_collection/target_data/bindingdb_P00918_data_report.json'

        # Convert the dictionary to JSON with indentation
        with open(output_path, 'w') as json_file:
            json.dump(report_dict, json_file, indent=4)



        print_boxplot(clean_df['pKi'], 'P00918')
        print_histogram(clean_df['pKi'], 'P00918')
