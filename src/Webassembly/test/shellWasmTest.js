console.log('Started');
let Module;

function instantiateModule(M) {
    console.log('Ready', M);
    Module = M;
    runFunc();
}


if (false) {
//if (typeof wasmFeatureDetect !== 'undefined') {
    console.log('wasmFeatureDetect Working');
    wasmFeatureDetect.threads().then(threadsSupported => {
        wasmFeatureDetect.simd().then(simdSupported => {
            if (threadsSupported) {
                console.log('threaded');
                if (simdSupported) {
                    console.log('vectorized');
                    ShellWasmTV().then(instantiateModule);
                } else {
                    console.log('unvectorized');
                    ShellWasmT().then(instantiateModule)
                }
            } else {
                console.log('unthreaded');
                if (simdSupported) {
                    console.log('vectorized');
                    ShellWasmV().then(instantiateModule);
                } else {
                    console.log('unvectorized');
                    ShellWasm().then(instantiateModule)
                }
            }
        });
    });
} else {
    console.log('wasFeatureDetect Not Working');
    console.log('unvectorized');
    ShellWasm().then(instantiateModule)
}

const shellData = [
    [.460, 780, .292, 1460, 2574, 6, .033, 76, 45, 60, 0],
    [.220, 995, .2549, 207, 2574, 7, .022, 37, 55, 65, 0]
];
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
const angles = [10];

const runFunc = () => {
    const shells = [];
    const calc = new Module.shellCalc(); 
    //WARNING: DO NOT DELETE IF MULTI-THREADED, WILL DEADLOCK
    calc.setMax(30);
    calc.setMin(0);
    calc.setPrecision(0.1);
    calc.setDtMin(0.01);

    for (const data of shellDataStruct) {
        const sp = new Module.shellParams();
        sp.setValues(data[0]);
        const dp = new Module.dispersionParams(data[1]);
        dp.setValues(data[1]);
        shells.push(new Module.shell(sp, dp, ''));
    }

    for (const shell of shells) {
        calc.calcImpact(shell);
        calc.calcAngles(shell, 70, -20);
        calc.calcDispersion(shell, Module.verticalTypes.horizontal.value);
        calc.calcPostPen(shell, 70, -20, angles, true, false);

        const hash = Module.generateShellHash(shell)
        console.log(hash, typeof hash);
        console.log(shell.maxDist());
    }

    const testPointArrays = (shell) => {
        console.log(Module.getImpactSizedPointArray(shell, 
            [Module.calcIndices.impact.value, Module.impactIndices.distance.value], 
            [Module.calcIndices.impact.value, Module.impactIndices.rawPen.value], 
        ));
        
        console.log(Module.getImpactSizedPointArray(shell, 
            [Module.calcIndices.impact.value, Module.impactIndices.distance.value], 
            [Module.calcIndices.angle.value, Module.angleIndices.ra0D.value, 0], 
        ));

        console.log(Module.getImpactSizedPointArray(shell, 
            [Module.calcIndices.impact.value, Module.impactIndices.distance.value], 
            [Module.calcIndices.post.value, Module.postPenIndices.x.value, 0], 
        ));

        console.log(Module.getImpactSizedPointArrayFuseStatus(shell, 
            [Module.calcIndices.impact.value, Module.impactIndices.distance.value], 
            [Module.calcIndices.post.value, Module.postPenIndices.x.value, 0], 
            0, true
        ));

        console.log(Module.getImpactSizedPointArrayFuseStatus(shell, 
            [Module.calcIndices.impact.value, Module.impactIndices.distance.value], 
            [Module.calcIndices.post.value, Module.postPenIndices.x.value, 0], 
            0, false
        ));

        console.log(Module.getImpactSizedPointArray(shell,
            [Module.calcIndices.impact.value, Module.impactIndices.distance.value],
            [Module.calcIndices.post.value, Module.postPenIndices.x.value, 0]
        ));

        console.log(Module.getImpactSizedPointArray(shell,
            [Module.calcIndices.impact.value, Module.impactIndices.distance.value],
            [Module.calcIndices.dispersion.value, Module.dispersionIndices.maxHorizontal.value]
        ));
        
        console.log(Module.getImpactSizedPointArray(shell,
            [Module.calcIndices.impact.value, Module.impactIndices.distance.value], 
            [Module.calcIndices.impact.value, Module.impactIndices.launchA.value],
        ));
        console.log(Module.getTrajectoryPointArrays(shell, 50));
    }

    for (const shell of shells) {
        //testPoints(shell);
        testPointArrays(shell);
    }
}