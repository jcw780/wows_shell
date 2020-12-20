console.log('Started');
let Module;

function testFunction() {
    //testFull();
    //runFunc();
    runFuncSplit();
}

function instantiateModule(M) {
    console.log('Ready', M);
    Module = M;
    testFunction();
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


function runFunc() {
    console.log('Entered');
    const numShells = shellData.length;
    const instance = new Module.shellCombined(numShells);
    for (const [i, shell] of shellData.entries()) {
        instance.setValues(...shell, i);
    }

    //console.log(instance);
    instance.calcImpact();
    //instance.printImpact();
    console.log("ran");
    instance.calcAngles(70, -20);
    instance.calcPostPen(70, -20, angles, true, false);

    console.log("ran");
    //instance.printPostPen();

    const testPoints = (instance, num) => {
        for (let i = 0; i < instance.getImpactSize(); i++) {
            console.log(i, instance.getAnglePoint(num, i, Module.angleIndices.distance.value),
                instance.getAnglePoint(num, i, Module.angleIndices.ra0.value)
            );
        }

        for (let i = 0; i < instance.getImpactSize(); i++) {
            //console.log(i);
            console.log(i, instance.getPostPenPoint(num, i, Module.postPenIndices.distance.value, 0),
                instance.getPostPenPoint(num, i, Module.postPenIndices.x.value, 0),
                instance.getPostPenPoint(num, i, Module.postPenIndices.xwf.value, 0)
            );
        }
    }

    const testPointArrays = (instance, num) => {
        console.log(instance.getImpactPointArray(
            num,
            Module.impactIndices.distance.value,
            Module.impactIndices.rawPen.value
        ));

        console.log(instance.getAnglePointArray(
            num,
            Module.angleIndices.distance.value,
            Module.angleIndices.ra0D.value
        ));

        console.log(instance.getPostPenPointArray(
            num, 0,
            Module.postPenIndices.distance.value,
            Module.postPenIndices.x.value
        ));

        console.log(instance.getPostPenPointArrayFuseStatus(
            num, true, 0,
            Module.postPenIndices.distance.value,
            Module.postPenIndices.x.value
        ));

        console.log(instance.getPostPenPointArrayFuseStatus(
            num, false, 0,
            Module.postPenIndices.distance.value,
            Module.postPenIndices.x.value
        ));
    }

    for (let num = 0; num < numShells; num++) {
        testPoints(instance, num);
        //testPointArrays(instance, num);
    }
    instance.delete();
}
const runFuncSplit = () => {
    const shells = [];
    const calc = new Module.shellCalc();
    calc.setMax(30);
    calc.setMin(0);
    calc.setPrecision(0.1);
    calc.setDtMin(0.01);

    /*for (const data of shellData) {
        shells.push(new Module.shell(...data, ''));
        console.log(Module.generateHash(...data));
    }*/

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
        calc.calcDispersion(shell);
        calc.calcPostPen(shell, 70, -20, angles, true, false);

        const hash = Module.generateShellHash(shell)
        console.log(hash, typeof hash);
        console.log(shell.maxDist());
    }

    const testPoints = (shell) => {
        let impactSize = shell.getImpactSize();
        let postpenSize = shell.getPostPenSize();
        for (let i = 0; i < impactSize; i++) {
            console.log(shell.getImpactPoint(i, Module.impactIndices.distance.value),
                shell.getImpactPoint(i, Module.impactIndices.rawPen.value));
        }
        console.log('done');

        for (let i = 0; i < impactSize; i++) {
            console.log(shell.getAnglePoint(i, Module.angleIndices.distance.value),
                shell.getAnglePoint(i, Module.angleIndices.ra0D.value));
        }
        console.log('done');

        for (let i = 0; i < postpenSize; i++) {
            console.log(shell.getPostPenPoint(i, Module.postPenIndices.distance.value, 0),
                shell.getPostPenPoint(i, Module.postPenIndices.x.value, 0),
                shell.getPostPenPoint(i, Module.postPenIndices.xwf.value, 0));
        }
        console.log('done');
    }

    const testPointArrays = (shell) => {
        console.log(shell.getImpactPointArray(
            Module.impactIndices.distance.value,
            Module.impactIndices.rawPen.value
        ));

        console.log(shell.getAnglePointArray(
            Module.angleIndices.distance.value,
            Module.angleIndices.ra0D.value
        ));

        console.log(shell.getPostPenPointArray(
            0,
            Module.postPenIndices.distance.value,
            Module.postPenIndices.x.value
        ));

        console.log(shell.getPostPenPointArrayFuseStatus(
            true, 0,
            Module.postPenIndices.distance.value,
            Module.postPenIndices.x.value
        ));

        console.log(shell.getPostPenPointArrayFuseStatus(
            false, 0,
            Module.postPenIndices.distance.value,
            Module.postPenIndices.x.value
        ));

        console.log(Module.getImpactSizedPointArray(shell,
            [Module.calcIndices.impact.value, Module.impactIndices.distance.value],
            [Module.calcIndices.post.value, Module.postPenIndices.x.value, 0]
        ));

        console.log(Module.getImpactSizedPointArray(shell,
            [Module.calcIndices.impact.value, Module.impactIndices.distance.value],
            [Module.calcIndices.dispersion.value, Module.dispersionIndices.maxHorizontal.value]
        ));
    }

    for (const shell of shells) {
        //testPoints(shell);
        testPointArrays(shell);
    }
}