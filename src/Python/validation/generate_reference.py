import generate_data
from validation_files_util import generate_reference_csvs

impact_dataframe, angles_dataframe, dispersion_dataframe, post_penetration_dataframe, trajectory_dataframe = generate_data.get_dataframes()
generate_reference_csvs(impact_dataframe, angles_dataframe, dispersion_dataframe, post_penetration_dataframe, trajectory_dataframe)