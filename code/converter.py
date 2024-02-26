import nbformat

def script_to_notebook(script_file, notebook_file):
    # Read the Python script
    with open(script_file, 'r') as f:
        script_content = f.read()

    # Create a new Jupyter Notebook
    notebook_content = nbformat.v4.new_notebook()

    # Add a code cell with the script content
    notebook_content.cells.append(nbformat.v4.new_code_cell(script_content))

    # Save the notebook to a file
    with open(notebook_file, 'w') as f:
        nbformat.write(notebook_content, f)



# Replace 'your_script.py' and 'output_notebook.ipynb' with your actual filenames
script_to_notebook('data_load.py', 'data_load.ipynb')
