
#define _USE_MATH_DEFINES
#include <math.h>
#include <omp.h>
#include "constant.h"

#include <iostream>
#include <string>
#include <vector>
#include <stdlib.h>
#include <cstring>

using namespace std;

class shell{
    private:
    double v0; // = 780;
    double caliber; // = .460;
    double krupp; // = 2574;
    double mass; // = 1460;
    double normalization; // = 6;
    double cD; // = .292;
    string name;

    double max = 10;
    double min = 0;
    double precision = .01;
    double x0 = 0, y0 = 0;
    double dt = .01;

    public:
    double k, cw_2, pPPC, normalizationR;
    private:
    vector<double> oneVector;
    vector<double> temp;

    /*standardOut1 - indicies refer to multiples of size
    [0:1] x [1:2] y [2:3] t [3:4] v_x [4:5] v_y
    */

    vector<double> standardOut0;

    /*trajectories output
    [0           ]trajx 0        [1           ]trajy 1
    ...
    [size * 2 - 2]trajx size - 1 [size * 2 - 1]trajy size - 1
    */
    vector<vector<double>> trajectories;

    /*standard output - Indicies refer to multiples of size [Num Colums of 12
    [0 : 1)-distance,          [1 : 2)-launch angle,      [2 : 3)-impact angle - R
    [3 : 4)-impact angle - D,  [4 : 5)-impact velocity,   [5 : 6)-raw pen,           
    [6 : 7)-effectivepenH,     [7 : 8)-effectivepenH w/N, [8 : 9)-IA_d,              
    [9 :10)-effectivepenD,     [10:11)-effectivepenD w/N, [11:12)-ttt,          
    [12:13)-ttta
    */

    public:
    vector<double> standardData;
    unsigned int size;

    double calcNormalizationR(double angle){ //Input in radians
        if(fabs(angle) > normalizationR){
            return fabs(angle) - normalizationR;
        }else{
            return 0;
        }
    }

    template <typename T> inline constexpr
    int signum(T x, std::false_type is_signed) {
        return T(0) < x;
    }

    template <typename T> inline constexpr
    int signum(T x, std::true_type is_signed) {
        return (T(0) < x) - (x < T(0));
    }

    template <typename T> inline constexpr
    int signum(T x) {
        return signum(x, std::is_signed<T>());
    }

    private:

    void preProcessVariables(){
        k = 0.5 * cD * pow((caliber/2),2) * M_PI / mass;
        cw_2 = 100+1000/3*caliber;
        pPPC = C * krupp/2400 * pow(mass,0.55) / pow((caliber*1000),0.65);
        normalizationR = normalization / 180 * M_PI;
    }

    void setArrays(){
        size = (unsigned int) (max - min) / precision;
        oneVector.resize(size);
        memset(oneVector.data(), 1, sizeof(double) * size);

        trajectories.resize(2 * size);
        standardData.resize(13 * size);

        standardOut0.resize(size * 4);
        temp.resize(size * 7);
    }
    
    void preProcessAV(){
        double angle, angleR;
        #pragma omp parallel for
        for(int i=0; i < size; i++){
            angle = i * precision + min;
            angleR = angle * M_PI / 180;
            standardData[i + size] = angle;
            standardOut0[i + 3 * size] = cos(angleR) * v0;
            standardOut0[i + 4 * size] = sin(angleR) * v0;
        }
    }

    void postProcessStandard(){
        //Copies first section - [0]x->[0]distance
        copy(standardOut0.begin(), standardOut0.begin() + size, standardData.begin());
        //Copies 11th section - [2]t->[10]ttt
        //copy(standardOut0.begin() + 2 * size, standardOut0.begin() + 3 * size, standardData.begin() + 11 * size);

        double iAR, iADR, iV, rP;
        #pragma omp parallel for private(iAR, iADR, iV, rP)
        for(int i=0; i<size; i++){
            //Calculate [2]IA , [7]IA_D
            iAR = atan(standardOut0[i+size*3]/standardOut0[i+size*2]);
            standardData[i+size*2] = iAR;
            standardData[i+size*3] = iAR / M_PI * 180;
            iADR = M_PI / 2 + iAR;
            standardData[i+size*8] = iADR / M_PI * 180;

            //Calculate [3]iV,  [4]rP
            iV = sqrt(pow(standardOut0[i+size*3],2) + pow(standardOut0[i+size*2],2));
            standardData[i+size*4] = iV;
            rP = pow(iV, 1.1) * pPPC;
            standardData[i+size*5] = rP;
            //Calculate [5]EPH  [8]EPV
            standardData[i+size*6] = cos(iAR) * rP;
            standardData[i+size*9] = cos(iADR) * rP;

            standardData[i+size*7] = cos(calcNormalizationR(iAR)) * rP;
            standardData[i+size*10] = cos(calcNormalizationR(iADR)) * rP;

            standardData[i+size*12] = standardData[i+size*11] / 3.1;
        }

    }

    void singleTraj(unsigned int i){
        //printf("%d ", i);
        //printf("%f %f %f %f\n", x, y, v_x, v_y);
        double T, p, rho, t, x, y, v_x, v_y;
        v_x = standardOut0[i+size*3];
        v_y = standardOut0[i+size*4];
        x = x0;
        y = y0;
        t = 0;
        trajectories[2*i  ].push_back(x);
        trajectories[2*i+1].push_back(y);
        //printf("%f ",dt);
        while(y >= 0){
            x = x + dt*v_x;
            y = y + dt*v_y;
            //printf("%f ", x);
            T = t0 - L*y;
            p = p0*pow((1-L*y/t0),(g*M/(R*L)));
            rho = p*M/(R*T);

            v_x = v_x - dt*k*rho*(cw_1*v_x*v_x+cw_2*v_x);
            v_y = v_y - dt*(g - k*rho*(cw_1*v_y*v_y+cw_2*fabs(v_y))*signum(v_y));
            t = t + dt;
            trajectories[2*i  ].push_back(x);
            trajectories[2*i+1].push_back(y);
        }
        standardOut0[i       ] = x;
        standardOut0[i+size  ] = y;
        //standardOut0[i+size*2] = t;
        standardData[i+size*11] = t;
        standardOut0[i+size*2] = v_x;
        standardOut0[i+size*3] = v_y;
    }


    public:
    shell(){
    }

    void assignShellValues(double v0, double caliber, double krupp, double mass,
    double normalization, double cD, string name){
        this->v0 = v0;
        this->caliber = caliber;
        this->krupp = krupp;
        this->mass = mass;
        this->normalization = normalization;
        this->cD = cD;
        this->name = name;
        preProcessVariables();
    }

    shell(double v0, double caliber, double krupp, double mass,
    double normalization, double cD, string name){
        assignShellValues(v0, caliber, krupp, mass, normalization, cD, name);
    }

    void editTestParameters(double max, double min, double precision, double x0, double y0, double dt){
        if(!max){
            this->max = max;
        }
        if(!min){
            this->min = max;
        }
        if(!precision){
            this->precision = max;
        }
        if(!x0){
            this->x0 = max;
        }
        if(!y0){
            this->y0 = max;
        }
        if(!dt){
            this->dt = max;
        }
    }

    void calculateStandard(){
        setArrays();
        preProcessAV();

        //printf("%d\n", size);
        omp_set_num_threads(6);
        #pragma omp parallel for
        for(int i=0; i<size; i++){
            singleTraj(i);
            //printf("i %d \n", i);
        }
        
        postProcessStandard();
    }

    void printStandardOut(){
        for(int i=0; i<size; i++){
            printf("%f %f %f %f %f\n", standardOut0[i], standardOut0[i+size], 
            standardOut0[i+size*2], standardOut0[i+size*3], standardOut0[i+size*4]);
        }
    }
    void printStandardData(){
        for(int i=0; i<size; i++){
            for(int j=0; j<13; j++){
                printf("%f ", standardData[i + size * j]);
            }
            printf("\n");
        }
    }
    void printTrajectory(unsigned int target){
        printf("Index:[%d] X Y\n", target);
        for(std::vector<double>::size_type i = 0; i != trajectories[target*2].size(); i++) {
            printf("%f %f\n", trajectories[target*2][i], trajectories[target*2+1][i]);
        }

    }

};

class postPen: public shell{
    private:
    double threshold, fuseTime, dtf = 0.0001;
    double xf0 = 0, yf0 = 0;
    bool completed = false;
    vector<double> velocities;
    
    bool includeNormalization = true;
    bool nChangeTrajectory = true;
    
    void initializeArrs(){
        postPenSize = size * angles.size();
        postPenData.resize(6 * postPenSize);
        velocities.resize(4 * postPenSize);

        #pragma omp parallel for
        for(int i=0; i<angles.size(); i++){
            printf("%f\n", angles[i]);
            for(int j=0; j<size; j++){
                printf("%f\n", standardData[j]);
            }
            fill(postPenData.begin() + i * size, postPenData.begin() + (i + 1) * size, angles[i]);
            copy(standardData.begin(), standardData.begin() + size, postPenData.begin() + postPenSize + i * size);
        }
    }

    void calcVelocities(double thickness){
        unsigned int distIndex, anglesIndex;
        double hAngle, vAngle, cAngle, nCAngle, aAngle, pPV, ePenetration, hFAngle, vFAngle;
        #pragma omp parallel for private(distIndex, anglesIndex, hAngle, vAngle, cAngle, nCAngle, aAngle, pPV, ePenetration, hFAngle, vFAngle)
        for(int i=0; i < size * angles.size(); i++){
            distIndex = i % size;
            anglesIndex = i / size;

            hAngle = angles[anglesIndex] /180*M_PI;
            vAngle = standardData[distIndex+size*2];
            cAngle = acos(cos(hAngle) * cos(vAngle));

            if(includeNormalization){
                nCAngle = calcNormalizationR(cAngle);
            }else{
                nCAngle = cAngle;
            }
                        
            ePenetration = standardData[distIndex+size*5]*cos(nCAngle);

            if(ePenetration > thickness){
                pPV = (1-exp(1-ePenetration/thickness)) * standardData[distIndex+size*4];
            }else{
                pPV = 0;
            }

            aAngle = nCAngle / cAngle;
            hFAngle = hAngle * aAngle;
            vFAngle = vAngle * aAngle;

            velocities[i              ] = pPV * cos(vFAngle) * cos(vFAngle);
            velocities[i+  postPenSize] = pPV * sin(vFAngle);
            velocities[i+2*postPenSize] = pPV * cos(vFAngle) * sin(hFAngle);
            velocities[i+3*postPenSize] = thickness / cos(nCAngle);
        }
    }

    void printVelocities(){
        printf("Velocities\n");
        for(int i=0; i<postPenSize; i++){
            for(int j=0; j<4; j++){
                printf("%f ", velocities[i + j * postPenSize]);
            }
            printf("\n");
        }
    }

    void postPenTraj(unsigned int i){
        //printf("%d ", i);
        //printf("%f %f %f %f\n", x, y, v_x, v_y);
        double T, p, rho, t, dt, x, y, z, v_x, v_y, v_z;
        v_x = velocities[i           ];
        v_z = velocities[i + postPenSize * 2];
        v_y = velocities[i + postPenSize    ];
        //printf("%d %f %f %f\n", i, v_x, v_y, v_z);
        x = xf0;
        y = yf0;
        z = xf0;
        t = 0;
        //printf("%f ",dt);
        if(v_x > 0){
            while(t < fuseTime){
                x = x + dtf*v_x;
                z = z + dtf*v_z;
                y = y + dtf*v_y;
                if(i == 100){
                    printf("%f %f %f %f\n", x, y, z, t);
                }
                //printf("%f ", x);
                T = t0 - L*y;
                p = p0*pow((1-L*y/t0),(g*M/(R*L)));
                rho = p*M/(R*T);

                v_x = v_x - dtf*k*rho*(cw_1*v_x*v_x+cw_2*v_x);
                v_z = v_z - dtf*k*rho*(cw_1*v_z*v_z+cw_2*v_z);
                v_y = v_y - dtf*(g - k*rho*(cw_1*v_y*v_y+cw_2*fabs(v_y))*signum(v_y));
                t = t + dtf;
            }
            postPenData[i + postPenSize * 2] = x;
            postPenData[i + postPenSize * 3] = y;
            postPenData[i + postPenSize * 4] = z;
            if(velocities[i + postPenSize * 3] > threshold){
                postPenData[i + postPenSize * 5] = x;
            }else{
                postPenData[i + postPenSize * 5] = -1;
            }
        }else{
            postPenData[i + postPenSize * 2] = 0;
            postPenData[i + postPenSize * 3] = 0;
            postPenData[i + postPenSize * 4] = 0;
            postPenData[i + postPenSize * 5] = 0;
        }
    }

    public:
    unsigned int postPenSize;
    vector<double> angles;
    vector<double> postPenData;

    postPen(double v0, double caliber, double krupp, double mass,
    double normalization, double cD, string name, double threshold, double fuseTime){
        this->fuseTime = fuseTime;
        assignShellValues(v0, caliber, krupp, mass, normalization, cD, name);
        if(!threshold){
            threshold = caliber / 6;   
        }
    }

    void calculateStandard(){
        completed = true;
        shell::calculateStandard();
        printf("Completed\n");
    }

    void setAngles(vector<double> angles){
        this->angles = angles;

    }

    void calculatePostPen(double thickness){
        initializeArrs();
        calcVelocities(thickness);

        omp_set_num_threads(6);
        #pragma omp parallel for
        for(int i=0; i<postPenSize; i++){
            postPenTraj(i);
        }
    }

    void printPostPen(){
        for(int i=0; i<postPenSize; i++){
            for(int j=0; j<6; j++){
                printf("%f ", postPenData[i + j * postPenSize]);
            }
            printf("\n");
        }
    }

};

int main(){
    //shell test(780, .460, 2574, 1460, 6, .033, .292, 76, "Yamato");
    //shell test(780, .460, 2574, 1460, 6, .292, "Yamato");
    postPen test(780, .460, 2574, 1460, 6, .292, "Yamato", 76, .033);
    test.calculateStandard();
    test.printStandardData();
    vector<double> angle;
    angle.push_back(0);
    angle.push_back(10);
    angle.push_back(20);
    test.setAngles(angle);
    test.calculatePostPen(100);
    //test.printTrajectory(2499);

}