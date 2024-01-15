#include <iostream>
#include <fstream>
#include <time.h>
#include<bits/stdc++.h>

#include "doublependulum.h"

#include "VtkData.h"
#include "VtkData.cpp"
using namespace std;

// Calculation of the angle 1
double DoublePendulum::calc1(double angle1, double a1, double h)
{
    double k1 = angle1;                                                 // t0
    double k2 = angle1 + k1 * 0.5 * h;                                  // t0 + h/2
    double k3 = angle1 + k2 * 0.5 * h;                                  // t0 + h/2
    double k4 = angle1 + k3 * h;                                        // t0 + h
    return angle1 + ( h * a1 * (k1+ 2*k2 + 2*k3 + k4) );
}

// Calculation of the angle 2
double DoublePendulum::calc2(double angle2, double a1, double h)
{
    double k1 = angle2;
    double k2 = angle2 + k1 * 0.5 * h;
    double k3 = angle2 + k2 * 0.5 * h;
    double k4 = angle2 + k3 * h;
    return angle2 + ( h * a1 * (k1+ 2*k2 + 2*k3 + k4) );

}

// Calculation of the angular velocity 1
double DoublePendulum::func3(double gravity, double mass1, double mass2, double length1, double length2, double angle1, double angle2,
                             double angVel1, double angVel2 ,double u1, double u2, double u3, double u4, double u5,double u7, double a1, double h)
{
    double M11 = (mass1+mass2) * length1;
    double M12 = (mass2 * length2) * u4;
    double M21 = length1 * u4;
    double M22 = length2;

    double F1 = ( (-1) * mass2 * length2 * angVel2 * angVel2 * u3) - (mass1 + mass2) * gravity * u1;
    double F2 = (length1 * angVel1 * angVel1 * u3 - (gravity * u7));


    double angAcc1 = 1.0 / (M11*M22 - M12*M21) * (M22*F1 - M12*F2);

    double k1 = angAcc1;                                                 // t0
    double k2 = angAcc1 + k1 * 0.5 * h;                                  // t0 + h/2
    double k3 = angAcc1 + k2 * 0.5 * h;                                  // t0 + h/2
    double k4 = angAcc1 + k3 * h;                                        // t0 + h

    return angAcc1 + ( h * a1 * (k1+ 2*k2 + 2*k3 + k4) );
}

// Calculation of the angular velocity 2
double DoublePendulum::func4(double gravity, double mass1, double mass2, double length1, double length2, double angle1, double angle2,
                             double angVel1, double angVel2, double u6, double u3, double u4, double u5, double u7, double u1, double a1, double h)
 {

    double M11 = (mass1+mass2) * length1;
    double M12 = (mass2 * length2) * u4;
    double M21 = length1 * u4;
    double M22 = length2;

    double F1 = ( (-1) * mass2 * length2 * angVel2 * angVel2 * u3) - (mass1 + mass2) * gravity * u1;
    double F2 = (length1 * angVel1 * angVel1 * u3 - (gravity * u7));


    double angAcc2 = 1.0 / (M11*M22 - M12*M21) * (-M21*F1 - M11*F2);

    double k1 = angAcc2;                                                 // t0
    double k2 = angAcc2 + k1 * 0.5 * h;                                  // t0 + h/2
    double k3 = angAcc2 + k2 * 0.5 * h;                                  // t0 + h/2
    double k4 = angAcc2 + k3 * h;                                        // t0 + h

    return angAcc2 + ( h * a1 * (k1+ 2*k2 + 2*k3 + k4) );
}

void DoublePendulum::iterationProcess(double angle1,double angle2, double angVel1, double angVel2, double iteration)
{
    // Define trigonometric terms in a short form
    double u1,u2,u3,u4,u5,u6,u7;
    double a1 = 0.166666;
    double cubicRootMass1 = cbrt(4);
    double cubicRootMass2 = cbrt(0.5);

    // Integrating initial angle values into angle vectors in order to use in interation
    double angle1_Value = angle1;
    double angle2_Value = angle2;
    m_angleVector1.push_back(angle1_Value);
    m_angleVector2.push_back(angle2_Value);

    // Integrating initial angular velocity values into angularVectors in order to use in interation
    double angularVeloctiy1_Value = angVel1;
    double angularVeloctiy2_Value = angVel2;
    m_angularVelocityVector1.push_back(angularVeloctiy1_Value);
    m_angularVelocityVector2.push_back(angularVeloctiy2_Value);

    double h = 0.01;

    // Iteration loop
    for(unsigned int i = 0; i<=iteration; i++)
    {
        // Iteration of angle 1 calculation
        angle1_Value = calc1(angle1_Value, a1, h);

        if(angle1_Value==0)
        {
            cout<< "The angle 1 value is approaching zero and become stable" << endl;
            break;
        }

        // Iteration of angle 2 calculation
         angle2_Value = calc2(angle2_Value, a1, h);

        if(angle2_Value==0)
        {
            cout<< "The angle 2 value is approaching zero and become stable" << endl;
            break;
        }

        u1 = sin(angle1_Value);
        u2 = sin(angle1_Value-2*angle2_Value);
        u3 = sin(angle1_Value-angle2_Value);
        u4 = cos(angle1_Value-angle2_Value);
        u5 = cos(2*angle1_Value-2*angle2_Value);
        u6 = cos(angle1_Value);
        u7 = sin(angle2_Value);

        m_angleVector1.push_back(angle1_Value);
        m_angleVector2.push_back(angle2_Value);

        // Iteration of the angular velocity 1 and 2
        m_angularVelocityVector1.push_back(func3(GRAVITY, MASS_1, MASS_2, LENGTH_1, LENGTH_2,
                    angle1_Value, angle2_Value, angularVeloctiy1_Value, angularVeloctiy2_Value, u1, u2, u3,  u4,  u5, u7, a1, h));

        m_angularVelocityVector2.push_back(func4(GRAVITY, MASS_1, MASS_2, LENGTH_1, LENGTH_2,
                    angle1_Value, angle2_Value, angularVeloctiy1_Value, angularVeloctiy2_Value, u6, u3,  u4, u5, u7, u1, a1, h));

        // Iteration of the Position of ball 1 with respect to x and y
        m_positionVectorBall1x.push_back(posCalc1x(LENGTH_1, angle1_Value));
        m_positionVectorBall1y.push_back(posCalc1y(LENGTH_1, angle1_Value));


        // Iteration of the Position of ball 2 with respect to x and y
        m_positionVectorBall2x.push_back(posCalc2x(LENGTH_2, angle2_Value, m_positionVectorBall1x.back()));
        m_positionVectorBall2y.push_back(posCalc2y(LENGTH_2, angle2_Value, m_positionVectorBall1y.back()));

        // Iteration of the Velocity of ball 1 with respect to x and y
        m_velocityVectorBall1x.push_back(velCalc1x(LENGTH_1, angle1_Value, m_angularVelocityVector1.back()));
        m_velocityVectorBall1y.push_back(velCalc1y(LENGTH_1, angle1_Value, m_angularVelocityVector1.back()));

        // Iteration of the Velocity of ball 2 WRT x and y
        m_velocityVectorBall2x.push_back(velCalc2x(LENGTH_2, angle2_Value, m_angularVelocityVector2.back(), m_velocityVectorBall1x.back()));
        m_velocityVectorBall2y.push_back(velCalc2y(LENGTH_2, angle2_Value, m_angularVelocityVector2.back(), m_velocityVectorBall1y.back()));

       // Writing VTK Files
       VtkData vtkData;

       string filename;
       filename = "file" + to_string(i) + ".vtk";

       vector < Vec3 > points;
       vector < vector <int> > cells;


       // Points Definition
       points.push_back(Vec3(0.0, 0.0, 0.0));
       points.push_back(Vec3(m_positionVectorBall1x.back(), m_positionVectorBall1y.back(), 0.0));
       points.push_back(Vec3(m_positionVectorBall2x.back(), m_positionVectorBall2y.back(), 0.0));

       // Cell Definition
       cells.push_back ( {0, 1} );
       cells.push_back ( {1, 2} );

       vector < double > scalarPointData;
       vector < Vec3 > vectorPointData;

       scalarPointData.push_back(0.0); // 0.0
       scalarPointData.push_back(cubicRootMass1); // cubic root of mass 1
       scalarPointData.push_back(cubicRootMass2);

       vectorPointData.push_back(Vec3(0.0, 0.0, 0.0)); // 0, 0, 0
       vectorPointData.push_back(Vec3(m_velocityVectorBall1x.back(), m_velocityVectorBall1y.back(), 0.0)); // velocity of first massx, y, 0
       vectorPointData.push_back(Vec3(m_velocityVectorBall2x.back(), m_velocityVectorBall2y.back(), 0.0)); // x, velocity of first massy, 0


       vtkData.setPoints (points);
       vtkData.setCells (cells);
       vtkData.setScalarPointData(scalarPointData);
       vtkData.setVectorPointData (vectorPointData);
       vtkData.writeVtkFile ( filename);

    }

}

void DoublePendulum::dumpAngle1()
{
    for(int i=0; i<m_angleVector1.size(); i++ ){
        cout << "angle 1" << "index[" << i << "]:" <<m_angleVector1[i] << endl;
    }
}

void DoublePendulum::dumpAngle2()
{
    for(int i=0; i<m_angleVector2.size(); i++ ){
        cout << "angle 2" << "index[" << i << "]:" <<m_angleVector2[i] << endl;
    }
}

void DoublePendulum::dumpAngVel1()
{
    for(int i=0; i< m_angularVelocityVector1.size(); i++ ){
        cout << "Angular Velocity 1" << "index[" << i << "]:" << m_angularVelocityVector1[i] << endl;
    }
}

void DoublePendulum::dumpAngVel2()
{
    for(int i=0; i< m_angularVelocityVector2.size(); i++ ){
        cout << "Angular Velocity 2" << "index[" << i << "]:" << m_angularVelocityVector2[i] << endl;
    }
}



void DoublePendulum::dumpAll()
{
    dumpAngle1();
    dumpAngle2();
    dumpAngVel1();
    dumpAngVel2();
}

double initialAngle1= M_PI/2; // Initial angle 1
double initialAngle2= M_PI/4; // Initial angle 2
double initialAngularVelocity1 = 0;
double initialAngularVelocity2 = 0;

const double iteration = 1000;


int main()
{

    DoublePendulum doublePend;
    doublePend.iterationProcess(initialAngle1, initialAngle2, initialAngularVelocity1, initialAngularVelocity2,  iteration);

    doublePend.dumpAngle1();
    doublePend.dumpAngle2();
    doublePend.dumpAngVel1();
    doublePend.dumpAngVel2();
}
