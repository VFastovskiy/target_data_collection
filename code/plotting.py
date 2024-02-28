import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd



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
    # plt.title('Comparison of Common Unique Values')
    # plt.xlabel('Columns')
    plt.ylabel('Number of Common Unique Values')

    plt.savefig('../results/step1_pdf_parsing/plots/common_vals.svg')




def plot_unique_vals(unique_values_dict, columns_list, sf_concatination_flag):
    # Convert the dictionary to a DataFrame for plotting
    unique_values_df = pd.DataFrame(unique_values_dict)

    # Plotting a bar chart with specified colors
    colors = sns.dark_palette("navy", n_colors=len(columns_list))
    ax = unique_values_df.plot(kind='bar', rot=0, color=colors)
    # plt.title('Comparison of Unique Values')
    plt.ylabel('Number of Unique Values')

    # Annotate each bar with a number of elements
    for p in ax.patches:
        ax.annotate(f'{int(p.get_height())}', (p.get_x() + p.get_width() / 2., p.get_height()),
                    ha='center', va='center', xytext=(0, 5), textcoords='offset points', clip_on=False)

    if sf_concatination_flag:
        plt.savefig('../results/step1_pdf_parsing/plots/joined/unique_vals_joined_scaffolds.svg')
    else:
        plt.savefig('../results/step1_pdf_parsing/plots/unique_vals_with_comparison.svg')

