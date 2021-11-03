let Module, calculator;
const chartIds = [
    'chart0', 'chart1', 'chart2', 'chart3', 'chart4', 'chart5', 'chart6', 'chart7', 'chart8', 'chart9'
];
let charts = chartIds.map(v => {
    let ctx = document.getElementById(v).getContext('2d'); 
    return new Chart(
        ctx, {
            type: 'scatter',
            data: {
                datasets: [{}]
            },
            options: {
                scales: {
                    xAxes: [{
                        type: 'linear',
                        position: 'bottom'
                    }],
                    yAxes: [{
                        id: 'A',
                        type: 'linear',
                        position: 'left'
                    },{
                        id: 'B',
                        type: 'linear',
                        position: 'right'
                    }],
                }
            }
        }
    );
});

let angles = [0, 10, 20];
let armor = 70;
let inclination = -20;

const testFunction = () => {
    const shellDataStruct = [
        [{
            caliber: .460, v0: 780, cD: .292, mass: 1460, krupp: 2574,
            normalization: 6, fuseTime: .033, threshold: 76, ricochet0: 45,
            ricochet1: 60, nonAP: 0
        }, {
            idealRadius: 10, minRadius: 2.8, idealDistance: 1000,
            taperDistance: 5000, delim: 0.5, zeroRadius: 0.2,
            delimRadius: 0.6, maxRadius: 0.8, maxDistance: 26630, sigma: 2.1
        }],
    ];
    
    const shells = shellDataStruct.map(data => {
        const sp = new Module.shellParams(data[0]);
        const dp = new Module.dispersionParams(data[1]);
        return new Module.shell(sp, dp, '');
    })

    shells.forEach(shell => {
        calculator.calcImpact(shell);
        calculator.calcDispersion(shell, Module.verticalTypes.horizontal.value);
        calculator.calcAngles(shell, armor, inclination);
        calculator.calcPostPen(shell, armor, inclination, angles, false, true);        
        let datasets = [];
        datasets.push(
            {
                label: 'Range',
                yAxisID: 'A',
                data: Module.getImpactSizedPointArray(shell,
                    [
                        Module.calcIndices.impact.value,
                        Module.impactIndices.launchA.value
                    ], 
                    [
                        Module.calcIndices.impact.value,
                        Module.impactIndices.distance.value
                    ]
                ) 
            },
        );

        charts[0].data.datasets = datasets;
        charts[0].update();

        datasets = [];
        datasets.push(
            {
                label: 'Impact Velocity',
                yAxisID: 'A',
                data: Module.getImpactSizedPointArray(shell,
                    [
                        Module.calcIndices.impact.value,
                        Module.impactIndices.distance.value
                    ], 
                    [
                        Module.calcIndices.impact.value,
                        Module.impactIndices.impactV.value
                    ]
                ) 
            },{
                label: 'Time to Target',
                yAxisID: 'B',
                data: Module.getImpactSizedPointArray(shell,
                    [
                        Module.calcIndices.impact.value,
                        Module.impactIndices.distance.value
                    ], 
                    [
                        Module.calcIndices.impact.value,
                        Module.impactIndices.tToTarget.value
                    ]
                ) 
            },{
                label: 'Time to Target Adjusted',
                yAxisID: 'B',
                data: Module.getImpactSizedPointArray(shell,
                    [
                        Module.calcIndices.impact.value,
                        Module.impactIndices.distance.value
                    ], 
                    [
                        Module.calcIndices.impact.value,
                        Module.impactIndices.tToTargetA.value
                    ]
                ) 
            }
        );
        charts[1].data.datasets = datasets;
        charts[1].update();

        datasets = [];
        datasets.push(
            {
                label: 'Raw Penetration',
                yAxisID: 'A',
                data: Module.getImpactSizedPointArray(shell,
                    [
                        Module.calcIndices.impact.value,
                        Module.impactIndices.distance.value
                    ], 
                    [
                        Module.calcIndices.impact.value,
                        Module.impactIndices.rawPen.value
                    ]
                ) 
            },{
                label: 'Effective Penetration',
                yAxisID: 'A',
                data: Module.getImpactSizedPointArray(shell,
                    [
                        Module.calcIndices.impact.value,
                        Module.impactIndices.distance.value
                    ], 
                    [
                        Module.calcIndices.impact.value,
                        Module.impactIndices.ePenH.value
                    ]
                ) 
            },{
                label: 'Effective Penetration Normalized',
                yAxisID: 'A',
                data: Module.getImpactSizedPointArray(shell,
                    [
                        Module.calcIndices.impact.value,
                        Module.impactIndices.distance.value
                    ], 
                    [
                        Module.calcIndices.impact.value,
                        Module.impactIndices.ePenHN.value
                    ]
                ) 
            },{
                label: 'Impact Angle Belt',
                yAxisID: 'B',
                data: Module.getImpactSizedPointArray(shell,
                    [
                        Module.calcIndices.impact.value,
                        Module.impactIndices.distance.value
                    ], 
                    [
                        Module.calcIndices.impact.value,
                        Module.impactIndices.impactAHD.value
                    ]
                ) 
            }
        );
        charts[2].data.datasets = datasets;
        charts[2].update();

        datasets = [];
        datasets.push(
            {
                label: 'Effective Penetration Deck',
                yAxisID: 'A',
                data: Module.getImpactSizedPointArray(shell,
                    [
                        Module.calcIndices.impact.value,
                        Module.impactIndices.distance.value
                    ], 
                    [
                        Module.calcIndices.impact.value,
                        Module.impactIndices.ePenD.value
                    ]
                ) 
            },{
                label: 'Effective Penetration Deck Normalized',
                yAxisID: 'A',
                data: Module.getImpactSizedPointArray(shell,
                    [
                        Module.calcIndices.impact.value,
                        Module.impactIndices.distance.value
                    ], 
                    [
                        Module.calcIndices.impact.value,
                        Module.impactIndices.ePenDN.value
                    ]
                ) 
            },{
                label: 'Impact Angle Deck',
                yAxisID: 'B',
                data: Module.getImpactSizedPointArray(shell,
                    [
                        Module.calcIndices.impact.value,
                        Module.impactIndices.distance.value
                    ], 
                    [
                        Module.calcIndices.impact.value,
                        Module.impactIndices.impactADD.value
                    ]
                ) 
            }
        );
        charts[3].data.datasets = datasets;
        charts[3].update();

        datasets = [];
        datasets.push(
            {
                label: 'Max Horizontal Dispersion',
                yAxisID: 'A',
                data: Module.getImpactSizedPointArray(shell,
                    [
                        Module.calcIndices.impact.value,
                        Module.impactIndices.distance.value
                    ], 
                    [
                        Module.calcIndices.dispersion.value,
                        Module.dispersionIndices.maxHorizontal.value
                    ]
                ) 
            },{
                label: 'Std Dev Horizontal Dispersion',
                yAxisID: 'A',
                data: Module.getImpactSizedPointArray(shell,
                    [
                        Module.calcIndices.impact.value,
                        Module.impactIndices.distance.value
                    ], 
                    [
                        Module.calcIndices.dispersion.value,
                        Module.dispersionIndices.standardHorizontal.value
                    ]
                ) 
            },{
                label: 'Half Horizontal Dispersion',
                yAxisID: 'A',
                data: Module.getImpactSizedPointArray(shell,
                    [
                        Module.calcIndices.impact.value,
                        Module.impactIndices.distance.value
                    ], 
                    [
                        Module.calcIndices.dispersion.value,
                        Module.dispersionIndices.halfHorizontal.value
                    ]
                ) 
            }
        );
        charts[4].data.datasets = datasets;
        charts[4].update();

        datasets = [];
        datasets.push(
            {
                label: 'Max Vertical Dispersion',
                yAxisID: 'A',
                data: Module.getImpactSizedPointArray(shell,
                    [
                        Module.calcIndices.impact.value,
                        Module.impactIndices.distance.value
                    ], 
                    [
                        Module.calcIndices.dispersion.value,
                        Module.dispersionIndices.maxVertical.value
                    ]
                ) 
            },{
                label: 'Std Dev Vertical Dispersion',
                yAxisID: 'A',
                data: Module.getImpactSizedPointArray(shell,
                    [
                        Module.calcIndices.impact.value,
                        Module.impactIndices.distance.value
                    ], 
                    [
                        Module.calcIndices.dispersion.value,
                        Module.dispersionIndices.standardVertical.value
                    ]
                ) 
            },{
                label: 'Half Horizontal Dispersion',
                yAxisID: 'A',
                data: Module.getImpactSizedPointArray(shell,
                    [
                        Module.calcIndices.impact.value,
                        Module.impactIndices.distance.value
                    ], 
                    [
                        Module.calcIndices.dispersion.value,
                        Module.dispersionIndices.halfVertical.value
                    ]
                ) 
            }
        );
        charts[5].data.datasets = datasets;
        charts[5].update();

        datasets = [];
        datasets.push(
            {
                label: 'Max Dispersion Area',
                yAxisID: 'A',
                data: Module.getImpactSizedPointArray(shell,
                    [
                        Module.calcIndices.impact.value,
                        Module.impactIndices.distance.value
                    ], 
                    [
                        Module.calcIndices.dispersion.value,
                        Module.dispersionIndices.maxArea.value
                    ]
                ) 
            },{
                label: 'Std Dev Dispersion Area',
                yAxisID: 'A',
                data: Module.getImpactSizedPointArray(shell,
                    [
                        Module.calcIndices.impact.value,
                        Module.impactIndices.distance.value
                    ], 
                    [
                        Module.calcIndices.dispersion.value,
                        Module.dispersionIndices.standardArea.value
                    ]
                ) 
            },{
                label: 'Half Dispersion Area',
                yAxisID: 'A',
                data: Module.getImpactSizedPointArray(shell,
                    [
                        Module.calcIndices.impact.value,
                        Module.impactIndices.distance.value
                    ], 
                    [
                        Module.calcIndices.dispersion.value,
                        Module.dispersionIndices.halfArea.value
                    ]
                ) 
            }
        );
        charts[6].data.datasets = datasets;
        charts[6].update();

        datasets = [];
        datasets.push(
            {
                label: 'Ricochet 0',
                yAxisID: 'A',
                data: Module.getImpactSizedPointArray(shell,
                    [
                        Module.calcIndices.impact.value,
                        Module.impactIndices.distance.value
                    ], 
                    [
                        Module.calcIndices.angle.value,
                        Module.angleIndices.ra0D.value
                    ]
                ) 
            },{
                label: 'Ricochet 1',
                yAxisID: 'A',
                data: Module.getImpactSizedPointArray(shell,
                    [
                        Module.calcIndices.impact.value,
                        Module.impactIndices.distance.value
                    ], 
                    [
                        Module.calcIndices.angle.value,
                        Module.angleIndices.ra1D.value
                    ]
                ) 
            },{
                label: 'Minimum Angle for Fusing',
                yAxisID: 'A',
                data: Module.getImpactSizedPointArray(shell,
                    [
                        Module.calcIndices.impact.value,
                        Module.impactIndices.distance.value
                    ], 
                    [
                        Module.calcIndices.angle.value,
                        Module.angleIndices.fuseD.value
                    ]
                ) 
            },{
                label: 'Maximum Angle for Perforation',
                yAxisID: 'A',
                data: Module.getImpactSizedPointArray(shell,
                    [
                        Module.calcIndices.impact.value,
                        Module.impactIndices.distance.value
                    ], 
                    [
                        Module.calcIndices.angle.value,
                        Module.angleIndices.armorD.value
                    ]
                ) 
            }
        );
        charts[7].data.datasets = datasets;
        charts[7].update();

        datasets = [];
        datasets.push(
            {
                label: 'Ricochet 0',
                yAxisID: 'A',
                data: Module.getImpactSizedPointArray(shell,
                    [
                        Module.calcIndices.impact.value,
                        Module.impactIndices.distance.value
                    ], 
                    [
                        Module.calcIndices.post.value,
                        Module.postPenIndices.x.value,
                        2
                    ]
                ) 
            },{
                label: 'Ricochet 1',
                yAxisID: 'A',
                data: Module.getImpactSizedPointArray(shell,
                    [
                        Module.calcIndices.impact.value,
                        Module.impactIndices.distance.value
                    ], 
                    [
                        Module.calcIndices.post.value,
                        Module.postPenIndices.y.value,
                        2
                    ]
                ) 
            },{
                label: 'Minimum Angle for Fusing',
                yAxisID: 'A',
                data: Module.getImpactSizedPointArray(shell,
                    [
                        Module.calcIndices.impact.value,
                        Module.impactIndices.distance.value
                    ], 
                    [
                        Module.calcIndices.post.value,
                        Module.postPenIndices.z.value,
                        2
                    ]
                ) 
            },{
                label: 'Maximum Angle for Perforation',
                yAxisID: 'A',
                data: Module.getImpactSizedPointArray(shell,
                    [
                        Module.calcIndices.impact.value,
                        Module.impactIndices.distance.value
                    ], 
                    [
                        Module.calcIndices.post.value,
                        Module.postPenIndices.xwf.value,
                        2
                    ]
                ) 
            }
        );
        charts[8].data.datasets = datasets;
        charts[8].update();
        
        let trajectories = Module.getTrajectoryPointArrays(shell, shell.getImpactSize()-1);
        console.log(trajectories);
        datasets = [];
        datasets.push(
            {
                label: 'Trajectory Original',
                yAxisID: 'A',
                data:  trajectories['original']
            },
            {
                label: 'Trajectory Compressed',
                yAxisID: 'A',
                data:  trajectories['compressed']
            }
        )

        charts[9].data.datasets = datasets;
        charts[9].update();
    })
}

ShellWasm().then(M => {
    Module = M;
    calculator = new Module.shellCalc();
    testFunction();
})