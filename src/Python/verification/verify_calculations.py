import pandas as pd

import generate_data
from verification_files_util import read_reference_csvs

impact_dataframe, angles_dataframe, dispersion_dataframe, post_penetration_dataframe, trajectory_dataframe = generate_data.get_dataframes()

impact_dataframe_reference, angles_dataframe_reference, dispersion_dataframe_reference, post_penetration_dataframe_reference, trajectory_dataframe_reference = read_reference_csvs()

def check_columns(test_dataframe, reference_dataframe):
    if len(test_dataframe.columns) != len(reference_dataframe.columns):
        print('Column Length Difference Found')
        return False
    reference_dataframe_columns_set = set(reference_dataframe.columns)

    name_check = True
    for name in test_dataframe.columns:
        if not name in reference_dataframe_columns_set:
            print(F'Column Label {name} not found')
            name_check = False

    return name_check 

def check_difference(test_dataframe, reference_dataframe):
    if not check_columns(test_dataframe, reference_dataframe):
        return (None, None)
    
    difference_dataframe = {}
    percent_difference_dataframe = {}
    for column_name, column_data in test_dataframe.items():
        difference = column_data - reference_dataframe[column_name]
        percent_difference = difference / reference_dataframe[column_name] * 100
        difference_dataframe[column_name] = difference
        percent_difference_dataframe[column_name] = percent_difference
        print(f'{column_name} Abs Diff: {difference.abs().max()} % Diff: {percent_difference.abs().max()}')
    
    return (pd.DataFrame(difference_dataframe), pd.DataFrame(percent_difference_dataframe))

def validate_difference(difference_dataframe):
    validation_dataframe = {}
    all_valid = True
    for column_name, column_data in difference_dataframe.items():
        max_difference = column_data.abs().max()
        valid = max_difference < 1e-6
        if not valid and all_valid:
            all_valid = False
        validation_dataframe[column_name] = (
            max_difference,
            valid
        )
    if all_valid:
        print("All Valid")
    else:
        print("Validation Failed")
    return pd.DataFrame(validation_dataframe)

def validate_dataframe(test_dataframe, reference_dataframe):
    difference, percent_difference = check_difference(test_dataframe, reference_dataframe)
    print(validate_difference(difference))

validate_dataframe(impact_dataframe, impact_dataframe_reference)
validate_dataframe(angles_dataframe, angles_dataframe_reference)
validate_dataframe(dispersion_dataframe, dispersion_dataframe_reference)
validate_dataframe(post_penetration_dataframe, post_penetration_dataframe_reference)
validate_dataframe(trajectory_dataframe, trajectory_dataframe_reference)

