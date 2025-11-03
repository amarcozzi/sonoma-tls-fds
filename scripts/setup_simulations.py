import sys
from pathlib import Path
from string import Template
from tqdm import tqdm

def setup_simulations():
    """
    Sets up simulation directories and input files based on .bdf files
    using the pathlib module.
    """
    try:
        # Define paths relative to the script's location
        script_path = Path(__file__).resolve()
        scripts_dir = script_path.parent
        root_dir = scripts_dir.parent
        aux_files_dir = root_dir / 'Auxiliary_Files'
        simulations_dir = root_dir / 'simulations'
        template_path = root_dir / 'template.fds'
    except NameError:
        # Handle case where __file__ is not defined (e.g., interactive interpreter)
        print("Error: This script is intended to be run as a file.", file=sys.stderr)
        return

    # Create the main simulations directory if it doesn't exist
    simulations_dir.mkdir(exist_ok=True)

    # Find all .bdf files in the Auxiliary_Files directory
    bdf_files = list(aux_files_dir.glob('*.bdf'))

    if not bdf_files:
        print(f"No .bdf files found in '{aux_files_dir}'")
        return

    identifiers = []

    # Read the template file
    if not template_path.is_file():
        print(f"Error: template.fds not found at '{template_path}'", file=sys.stderr)
        return
        
    template_content = template_path.read_text()
    fds_template = Template(template_content)

    print("Setting up simulation directories...")
    # Loop through each .bdf file and create the corresponding simulation setup
    for bdf_path in tqdm(bdf_files, desc="Processing files"):
        bdf_filename = bdf_path.name
        
        # Extract the identifier (e.g., 'c4_p26') from the filename
        identifier = '_'.join(bdf_filename.split('_')[:2])
        identifiers.append(identifier)

        # Create the specific simulation directory
        sim_path = simulations_dir / identifier
        sim_path.mkdir(exist_ok=True)

        # Substitute placeholders in the template
        new_fds_content = fds_template.substitute(
            identifier=identifier,
            bdf_file=bdf_filename
        )

        # Write the new input.fds file
        output_fds_path = sim_path / 'input.fds'
        output_fds_path.write_text(new_fds_content)

    print(f"\nSuccessfully created {len(identifiers)} simulation cases.")
    
    # Format identifiers for bash array
    # Using sorted() to ensure a consistent order
    bash_array_string = ' '.join(f'"{i}"' for i in sorted(identifiers))
    print("\nIdentifiers for bash script array:")
    print(f"({bash_array_string})")


if __name__ == "__main__":
    setup_simulations()