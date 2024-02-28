import glob
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from plotting import plot_common_vals, plot_unique_vals
import plotly.express as px
import networkx as nx
import plotly
import matplotlib.gridspec as gridspec
import numpy as np




def vz_pdb_id_distribution(csv_paths, sf_concatination_flag, columns_list):

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
            df_list, unique_values_dict = csv_to_unique_vals_dict(csv_paths, columns_list)
            plot_unique_vals(unique_values_dict, columns_list, sf_concatination_flag)
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
    plot_unique_vals(unique_values_dict, columns_list, False)







def vz_intersections(csv_paths, columns_list):
    df_list = []
    common_values = set()

    # Check if any value in a row is in the common values list
    def check_common_values(row):
        return any(val in row.values for val in common_values)

    for csv in csv_paths:
        df = pd.read_csv(csv)
        df = df[columns_list]
        df_list.append(df)
        df_name = os.path.splitext(os.path.basename(csv))[0]
        #print(df_name)

    common_vals_dict = csv_to_common_vals_dict(df_list, columns_list)
    for values_list in common_vals_dict.values():
        common_values.update(values_list)

    #print(common_values)

    for i, df in enumerate(df_list):
        name_for_all_dataset = f'all_dataset_{i}'
        # Apply the check_common_values function to filter rows
        filtered_df = df[df.apply(check_common_values, axis=1)]
        vz_graph(df, columns_list, name_for_all_dataset, common_values)
        vz_graph(filtered_df, columns_list, str(i), common_values)

        # Create a dictionary to map specific elements to colors
        color_discrete_map = {}
        for val in common_values:
            color_discrete_map[val] = 'yellow'
        #print(color_discrete_map)

        fig = px.sunburst(filtered_df, path=columns_list, color_discrete_map=color_discrete_map)

        # Add annotations for levels of hierarchy
        annotations = set_annotation(columns_list, filtered_df)
        fig.update_layout(annotations=annotations)

        output_path = f'../results/step1_pdf_parsing/plots/{i}_tree_diagram_common.svg'
        plotly.io.write_image(fig, output_path, format='svg')










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

    plt.figure()
    pos = nx.fruchterman_reingold_layout(G, scale=8.0, k=0.10, seed=131)

    nx.draw(G, pos, with_labels=False, font_size=3, node_color=node_colors, cmap=plt.cm.spring,
            node_size=10, font_color='black', font_weight='bold', width=1, edge_color=edge_colors)


    # Draw node labels above the nodes for level 1 nodes
    labels_level_1 = {node: node for node, level in levels.items() if level == 1}
    node_label_color_level_1 = 'black'
    nx.draw_networkx_labels(G, pos, labels_level_1, font_size=5, font_color=node_label_color_level_1,
                            font_weight='bold', verticalalignment='bottom')

    # Draw node labels above the nodes for level 2 nodes
    labels_level_2 = {node: node for node, level in levels.items() if level == 2}
    node_label_colors_level_2 = 'green'
    nx.draw_networkx_labels(G, pos, labels_level_2, font_size=5, font_color=node_label_colors_level_2,
                            font_weight='bold', verticalalignment='bottom')

    # Draw node labels above the nodes for level 3 nodes
    labels_level_3 = {node: node for node, level in levels.items() if level == 3}
    node_label_colors_level_3 = 'purple'
    nx.draw_networkx_labels(G, pos, labels_level_3, font_size=5, font_color= node_label_colors_level_3,
                            font_weight='bold', verticalalignment='bottom')

    # Save the graph as an SVG file
    output_path = f'../results/step1_pdf_parsing/plots/{df_name}_common_vals_graph.svg'
    plt.savefig(output_path, format='svg')



    # STEP 2
    # Identify connected components (subgraphs)
    subgraphs = list(nx.connected_components(G))

    # Define the layout for the main figure with graph and subplots
    fig = plt.figure(figsize=(12, 8))
    gs = gridspec.GridSpec(2, 2, width_ratios=[1, 1], height_ratios=[1, 0.5])

    # Draw each connected component separately
    for i, subgraph_nodes in enumerate(subgraphs):
        subgraph = G.subgraph(subgraph_nodes)

        levels = nx.get_node_attributes(subgraph, 'level')
        colors = [levels[node] for node in subgraph.nodes]

        # Create a list of node colors based on common_values
        node_colors = ['red' if node in common_values else 'grey' for node in subgraph.nodes]

        # Assign distinct colors to edges based on their levels
        # edge_colors = [levels[u] for u, v in subgraph.edges]

        # Define the position for the subplots
        ax0 = plt.subplot(gs[0, 0])
        ax1 = plt.subplot(gs[0, 1])
        ax2 = plt.subplot(gs[1, :])

        # Save each subgraph as an SVG file
        output_path = f'../results/step1_pdf_parsing/plots/{df_name}_subgraph_{i + 1}.svg'

        # Adjust the scale parameter to control the length of edges
        pos = nx.fruchterman_reingold_layout(subgraph, scale=0.5, k=0.05, seed=131)

        # Draw subgraph on the main figure
        nx.draw(subgraph, pos, with_labels=False, font_size=14, node_color=node_colors, cmap=plt.cm.spring,
                node_size=10, font_color='black', font_weight='bold', width=0.5, edge_color='black', ax=ax0)

        # Draw node labels above the nodes for each level
        for level_value in set(levels.values()):
            labels_level = {node: node for node, level in levels.items() if level == level_value}
            node_label_color_level = 'black' if level_value == 1 else 'green' if level_value == 2 else 'purple'
            nx.draw_networkx_labels(G, pos, labels_level, font_size=8 if level_value == 1 else 5,
                                    font_color=node_label_color_level, font_weight='bold', verticalalignment='bottom',
                                    ax=ax0)

        # Plot degree rank plot
        degree_sequence = sorted((d for n, d in subgraph.degree()), reverse=True)
        ax1.plot(degree_sequence, "b-", marker="o")
        ax1.set_title("Degree Rank Plot")
        ax1.set_ylabel("Degree")
        ax1.set_xlabel("Rank")

        # Plot degree histogram
        ax2.bar(*np.unique(degree_sequence, return_counts=True))
        ax2.set_title("Degree Histogram")
        ax2.set_xlabel("Degree")
        ax2.set_ylabel("# of Nodes")

        # Save the main figure with graph and subplots as an SVG file
        plt.savefig(output_path, format='svg')

        # Clear the figure for the next iteration
        plt.clf()




def vz_tree_diagram(csv_paths, columns_list):
    df_list = []

    vz_flag = True

    if vz_flag:
        for csv in csv_paths:
            df = pd.read_csv(csv)
            df_list.append(df)
            df_name = os.path.splitext(os.path.basename(csv))[0]

            # Group by 'Protein_Name' and count the number of unique 'PDB_ID'
            pdb_counts = df.groupby('Protein_Name')['PDB_ID'].nunique().reset_index()

            # Sort the DataFrame by the number of PDB IDs in descending order
            pdb_counts = pdb_counts.sort_values(by='PDB_ID', ascending=False)

            # Create a new DataFrame with protein names in sorted order
            protein_names_df = pdb_counts[['Protein_Name']]


            if len(protein_names_df) < 50:
                # Create a single sunburst diagram
                build_tree_without_split(pdb_counts, df, protein_names_df, columns_list, df_name)

            else:
                # Split DataFrame into chunks of 4 protein_names and create sunburst diagram for each chunk
                build_tree_with_split(protein_names_df, df, columns_list, df_name)





def build_tree_without_split(pdb_counts, df, protein_names_df, columns_list, df_name):
    # Create a single sunburst diagram
    top_proteins = pdb_counts.head(len(protein_names_df))
    fig = px.sunburst(df, path=columns_list)

    # Add annotations for levels of hierarchy
    annotations = set_annotation(columns_list, df)
    fig.update_layout(annotations=annotations)

    output_path = f'../results/step1_pdf_parsing/plots/{df_name}_tree_diagram.svg'
    plotly.io.write_image(fig, output_path, format='svg')



def build_tree_with_split(protein_names_df, df, columns_list, df_name):
    chunks = [protein_names_df[i:i + 4] for i in range(0, len(protein_names_df), 4)]

    for i, chunk_proteins in enumerate(chunks, 1):
        # Plot sunburst diagram for each chunk
        chunk_df = df[df['Protein_Name'].isin(chunk_proteins['Protein_Name'])]
        fig = px.sunburst(chunk_df, path=columns_list)

        # Add annotations for levels of hierarchy
        annotations = set_annotation(columns_list, df)
        fig.update_layout(annotations=annotations)

        output_path = f'../results/step1_pdf_parsing/plots/{df_name}_tree_diagram_part_{i}.svg'
        plotly.io.write_image(fig, output_path, format='svg')



def build_tree_for_common_elements(common_elements_dict, df_name, columns_list, df_list):
    for df in df_list:
        filtered_df = df[columns_list]

        # Lists of values for each column
        protein_name_list = common_elements_dict['Protein_Name']
        pdb_id_list = common_elements_dict['PDB_ID']
        chemical_id_list = common_elements_dict['Chemical_ID']

        # Initialize a mask with True values
        mask = pd.Series([True] * len(filtered_df), index=filtered_df.index)

        # Apply filtering for each condition
        for column, values in common_elements_dict.items():
            mask &= filtered_df[column].isin(values)

        # Filter the DataFrame based on the combined mask
        plot_df = filtered_df[mask]

        fig = px.sunburst(plot_df, path=columns_list)

        # Add annotations for levels of hierarchy
        annotations = set_annotation(columns_list, plot_df)
        fig.update_layout(annotations=annotations)

        output_path = f'../results/step1_pdf_parsing/plots/{df_name}_tree_diagram.svg'
        plotly.io.write_image(fig, output_path, format='svg')
        #print(plot_df)





def set_annotation(columns_list, df):
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

    return annotations





def visualization(vz_flag, sf_concatination_flag, csv_directory):
    csv_paths = glob.glob(os.path.join(csv_directory, '*.csv'))
    columns_list = ['Protein_Name', 'PDB_ID', 'Chemical_ID']

    if vz_flag & sf_concatination_flag:
        vz_pdb_id_distribution(csv_paths, True, columns_list)

    elif vz_flag and not sf_concatination_flag:
        # Plotting reports and writing DataFrames to CSV
        vz_comparison_of_csvs(csv_paths, columns_list)
        vz_pdb_id_distribution(csv_paths, False, columns_list)
        vz_tree_diagram(csv_paths, columns_list)
        vz_intersections(csv_paths, columns_list)

    else:
        print('Visualization was skipped or done previously.')