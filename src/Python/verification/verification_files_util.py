import pandas as pd

reference_directory = 'verification_reference'

def generate_csv(df, filename: str):
    with open(filename, 'w') as f:
        df.to_csv(f, index=False, line_terminator='\n')

reference_file_names = [
    'impact_data.csv',
    'angles_data.csv',
    'dispersion_data.csv',
    'postpenetration_data.csv',
    'trajectory_data.csv'
]

def generate_reference_csvs(*args):
    if len(args) == len(reference_file_names):
        for i, dataframe in enumerate(args):
            generate_csv(dataframe, f'{reference_directory}/{reference_file_names[i]}')
    else:
        raise TypeError(f'Missing arguments, expected {len(reference_file_names)}; found {len(args)}')

def read_reference_csvs():
    return [pd.read_csv(f'{reference_directory}/{file}') for file in reference_file_names]