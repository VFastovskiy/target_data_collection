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