from cgi import test
import numpy as np
import pandas as pd

import wows_shell

reference_directory = 'verification_reference'

def generate_csv(df, filename: str):
    with open(filename, 'w') as f:
        df.to_csv(f, index=False, line_terminator='\n')

def generate_reference_csvs(impact_dataframe, angles_dataframe, dispersion_dataframe, post_penetration_dataframe):
    generate_csv(impact_dataframe, f'{reference_directory}/impact_data.csv')
    generate_csv(angles_dataframe, f'{reference_directory}/angles_data.csv')
    generate_csv(dispersion_dataframe, f'{reference_directory}/dispersion_data.csv')
    generate_csv(post_penetration_dataframe, f'{reference_directory}/postpenetration_data.csv')

def generate_ordered_enum_list(enum_dict):
    ordered = [None] * len(enum_dict)
    for name, data in enum_dict.items():
        ordered[int(data)] = name
    
    return ordered

def read_reference_csvs():
    return (
        pd.read_csv(f'{reference_directory}/impact_data.csv'),
        pd.read_csv(f'{reference_directory}/angles_data.csv'),
        pd.read_csv(f'{reference_directory}/dispersion_data.csv'),
        pd.read_csv(f'{reference_directory}/postpenetration_data.csv')
    )

test_shell = wows_shell.shell(
    wows_shell.shellParams(.460, 780, .292, 1460, 2574, 6, .033, 76, 45, 60, 0),
    wows_shell.dispersionParams(10, 2.8, 1000, 5000, 0.5, 0.2, 0.6, 0.8, 26630, 2.1),
    "Yamato")

calculator = wows_shell.shellCalc()
calculator.setMax(25)
calculator.setMin(0)
calculator.setPrecision(0.1)
calculator.setX0(0)
calculator.setY0(0)
calculator.setDtMin(0.02)

calculator.calcImpactForwardEuler(test_shell, addTraj=True)
calculator.calcAngles(test_shell, 410, -15)
calculator.calcDispersion(test_shell, int(wows_shell.verticalTypes.horizontal))
angles = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60]
calculator.calcPostPen(test_shell, 70, 0, angles, True, False)

impact_data = test_shell.getImpact()
#print(impact_data.shape)

impact_dataframe = pd.DataFrame(impact_data.transpose(), 
    columns=generate_ordered_enum_list(wows_shell.impactIndices.__members__))

#print(impact_dataframe)

angles_data = test_shell.getAngles()
#print(angles_data.shape)

angles_dataframe = pd.DataFrame(angles_data.transpose(),
    columns=generate_ordered_enum_list(wows_shell.angleIndices.__members__))

#print(angles_dataframe)

dispersion_data = test_shell.getDispersion()
#print(dispersion_data.shape)

dispersion_dataframe = pd.DataFrame(dispersion_data.transpose(),
    columns=generate_ordered_enum_list(wows_shell.dispersionIndices.__members__))

#print(dispersion_dataframe)

post_penetration_data = test_shell.getPostPen()
post_penetration_data_shape = post_penetration_data.shape
#print(post_penetration_data_shape)
#print(post_penetration_data)

post_penetration_dataframe = pd.DataFrame(
    post_penetration_data
        .reshape((post_penetration_data_shape[0], post_penetration_data_shape[1] * post_penetration_data_shape[2]))
        .transpose(),
    columns=generate_ordered_enum_list(wows_shell.postPenIndices.__members__))
#print(post_penetration_dataframe)

#generate_reference_csvs(impact_dataframe, angles_dataframe, dispersion_dataframe, post_penetration_dataframe)

#Validation Phase

impact_dataframe_reference, angles_dataframe_reference, dispersion_dataframe_reference, post_penetration_dataframe_reference = read_reference_csvs()

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

