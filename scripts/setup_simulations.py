import sys
from pathlib import Path
from string import Template
from tqdm import tqdm

FOLIAGE_SUFFIX = '_foliage.bdf'
WOOD_SUFFIX = '_wood.bdf'


def find_pairs(aux_files_dir):
    """
    Maps each identifier (e.g. 'c4_p26') to its (foliage, wood) .bdf filenames.
    Exits with an error if a foliage file has no matching wood file or vice versa.
    """
    # Extract the identifier (e.g. 'c4_p26') from the filename
    def identifier_of(path):
        return '_'.join(path.name.split('_')[:2])

    foliage = {identifier_of(p): p.name for p in aux_files_dir.glob('*' + FOLIAGE_SUFFIX)}
    wood = {identifier_of(p): p.name for p in aux_files_dir.glob('*' + WOOD_SUFFIX)}

    unpaired = sorted(set(foliage) ^ set(wood))
    if unpaired:
        for identifier in unpaired:
            if identifier in foliage:
                print(f"Error: '{foliage[identifier]}' has no matching *{WOOD_SUFFIX}", file=sys.stderr)
            else:
                print(f"Error: '{wood[identifier]}' has no matching *{FOLIAGE_SUFFIX}", file=sys.stderr)
        sys.exit(1)

    return {identifier: (foliage[identifier], wood[identifier]) for identifier in foliage}


def setup_simulations():
    """
    Sets up simulation directories, input files, and an identifiers list based on
    paired foliage/wood .bdf files
    """
    try:
        # Define paths relative to the script's location
        script_path = Path(__file__).resolve()
        scripts_dir = script_path.parent
        root_dir = scripts_dir.parent
        aux_files_dir = root_dir / 'Auxiliary_Files'
        simulations_dir = root_dir / 'simulations'
        template_path = root_dir / 'template.fds'
        identifiers_path = root_dir / 'identifiers.txt'
    except NameError:
        # Handle case where __file__ is not defined (e.g., interactive interpreter)
        print("Error: This script is intended to be run as a file.", file=sys.stderr)
        return

    # Pair each *_foliage.bdf with its *_wood.bdf in the Auxiliary_Files directory
    pairs = find_pairs(aux_files_dir)

    if not pairs:
        print(f"No *{FOLIAGE_SUFFIX} / *{WOOD_SUFFIX} pairs found in '{aux_files_dir}'")
        return

    # Read the template file
    if not template_path.is_file():
        print(f"Error: template.fds not found at '{template_path}'", file=sys.stderr)
        return

    template_content = template_path.read_text()
    fds_template = Template(template_content)

    # Create the main simulations directory if it doesn't exist
    simulations_dir.mkdir(exist_ok=True)

    print("Setting up simulation directories...")
    # Loop through each foliage/wood pair and create the corresponding simulation setup
    for identifier, (foliage_bdf, wood_bdf) in tqdm(sorted(pairs.items()), desc="Processing pairs"):
        # Create the specific simulation directory
        sim_path = simulations_dir / identifier
        sim_path.mkdir(exist_ok=True)

        # Substitute placeholders in the template
        new_fds_content = fds_template.substitute(
            identifier=identifier,
            foliage_bdf=foliage_bdf,
            wood_bdf=wood_bdf,
        )

        # Write the new input.fds file
        output_fds_path = sim_path / 'input.fds'
        output_fds_path.write_text(new_fds_content)

    # Write the identifiers to a file for the Slurm script
    sorted_identifiers = sorted(pairs)
    identifiers_path.write_text('\n'.join(sorted_identifiers) + '\n')

    print(f"\nSuccessfully created {len(pairs)} simulation cases.")
    print(f"A list of identifiers has been saved to: {identifiers_path}")


if __name__ == "__main__":
    setup_simulations()
