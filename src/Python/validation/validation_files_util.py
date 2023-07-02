import pandas as pd

reference_directory = 'verification_reference'

def generate_csv(df, filename: str):
    with open(filename, 'w') as f:
        df.to_csv(f, index=False, line_terminator='\n')

def generate_reference_csvs(impact_dataframe, angles_dataframe, dispersion_dataframe, post_penetration_dataframe):
    generate_csv(impact_dataframe, f'{reference_directory}/impact_data.csv')
    generate_csv(angles_dataframe, f'{reference_directory}/angles_data.csv')
    generate_csv(dispersion_dataframe, f'{reference_directory}/dispersion_data.csv')
    generate_csv(post_penetration_dataframe, f'{reference_directory}/postpenetration_data.csv')

def read_reference_csvs():
    return (
        pd.read_csv(f'{reference_directory}/impact_data.csv'),
        pd.read_csv(f'{reference_directory}/angles_data.csv'),
        pd.read_csv(f'{reference_directory}/dispersion_data.csv'),
        pd.read_csv(f'{reference_directory}/postpenetration_data.csv')
    )