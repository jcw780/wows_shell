import pandas as pd

import wows_shell

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

def generate_ordered_enum_list(enum_dict):
    ordered = [None] * len(enum_dict)
    for name, data in enum_dict.items():
        ordered[int(data)] = name
    
    return ordered

def get_dataframes():
    impact_data = test_shell.getImpact()
    angles_data = test_shell.getAngles()
    dispersion_data = test_shell.getDispersion()
    post_penetration_data = test_shell.getPostPen()
    trajectory_data = test_shell.getTrajectory(250)

    impact_dataframe = pd.DataFrame(impact_data.transpose(), 
        columns=generate_ordered_enum_list(wows_shell.impactIndices.__members__))

    angles_dataframe = pd.DataFrame(angles_data.transpose(),
        columns=generate_ordered_enum_list(wows_shell.angleIndices.__members__))

    dispersion_dataframe = pd.DataFrame(dispersion_data.transpose(),
        columns=generate_ordered_enum_list(wows_shell.dispersionIndices.__members__))

    post_penetration_data = test_shell.getPostPen()
    post_penetration_data_shape = post_penetration_data.shape

    post_penetration_dataframe = pd.DataFrame(
        post_penetration_data
            .reshape((post_penetration_data_shape[0], post_penetration_data_shape[1] * post_penetration_data_shape[2]))
            .transpose(),
        columns=generate_ordered_enum_list(wows_shell.postPenIndices.__members__))
    
    trajectory_dataframe = pd.DataFrame(
        trajectory_data.transpose(),
        columns=['x', 'y', 'y compressed']
    )
    
    return (impact_dataframe, angles_dataframe, dispersion_dataframe, post_penetration_dataframe, trajectory_dataframe)