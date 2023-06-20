import pandas as pd
from pathlib import Path
import xarray as xr


def get_example_data_file_path(filename, data_dir='example_data'):
    # Path(__file__).resolve() returns the absolute path of the current module.
    # __file__ is a special variable path of the current file (here
    # data_io.py). We can use it as base path to construct other paths
    # that should end up correct on other machines or when the package is
    # installed
    module_directory = Path(__file__).resolve().parent
    data_path = Path(module_directory, data_dir, filename)
    return data_path


def load_data(data_file):
    if str(data_file).endswith('.txt'):
        input_file = pd.read_csv(data_file, sep=' ')
        output_data = input_file.to_xarray()
    elif str(data_file).endswith('.nc'):
        output_data = xr.open_dataset(data_file)
    else:
        raise TypeError('file must be a text or NetCDF file')
    return output_data
