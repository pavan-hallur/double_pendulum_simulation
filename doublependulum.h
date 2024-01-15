#pragma once

#include <math.h>
#include <string>
#include <sstream>
#include <vector>

// pendulum parameters
#define MASS_1 4
#define MASS_2 (double)0.5

#define GRAVITY (double)9.81

#define LENGTH_1 1
#define LENGTH_2 (double)0.5

// Initial conditions
#define INITITAL_H (double)0.01


class DoublePendulum {
public:
    DoublePendulum()
    {
//        m_angularVelocityVector1.push_back(initialAngVel1);
//        m_angularVelocityVector2.push_back(initialAngVel2);
    }

    ~DoublePendulum(){}
    void iterationProcess(double angle1,double angle2, double angVel1, double angVel2, double iteration);

    void dumpAngle1();
    void dumpAngle2();
    void dumpAngVel1();
    void dumpAngVel2();

    void dumpAll();

private:

    double calc1(double angle1, double a1, double h);

    double calc2(double angle2, double a1, double h);

    // Angular Velocity for ball 1
    double func3(double gravity, double mass1, double mass2, double length1, double length2, double angle1, double angle2,
                 double angVel1, double angVel2 ,double u1, double u2, double u3, double u4, double u5, double u7, double a1, double h);

    // Angular Velocity for ball 2
    double func4(double gravity, double mass1, double mass2, double length1, double length2, double angle1, double angle2,
                 double angVel1, double angVel2, double u6, double u3, double u4, double u5, double u7, double u1, double a1, double h);

    //Position 1 with respect to x
    double posCalc1x(double length1, double angle1)
    {
        return length1 * sin(angle1);
    }

    //Position 1 with respect to y
    double posCalc1y(double length1, double angle1)
    {
        return length1 * (-1) * cos(angle1);
    }
    //Position 2 WRT x
    double posCalc2x(double length2, double angle2, double pos1x)
    {
        return pos1x + length2 * sin(angle2);
    }
    //Position 2 WRT y
    double posCalc2y(double length2, double angle2, double pos1y)
    {
        return pos1y + length2 * (-1) * cos(angle2);
    }

    // Calculation of the velocities regarding ball 1 and ball 2
    //Velocity 1 with respect to x
    double velCalc1x(double length1, double angle1, double angVel1)
    {
        return length1 * cos(angle1) * angVel1;

    }
    //Velocity 1 WRT y
    double velCalc1y(double length1, double angle1, double angVel1)
    {
        return length1 * cos(angle1) * angVel1;
    }
    //Velocity 2 WRT x
    double velCalc2x(double length2, double angle2, double angVel2, double velocity1x)
    {
        return velocity1x + length2 * cos(angle2) * angVel2;
    }
    //Velocity 2 with respect to y
    double velCalc2y(double length2, double angle2, double angVel2, double velocity1y)
    {
        return velocity1y + length2 * sin(angle2) * angVel2;
    }


private:
    std::vector <double> m_positionVectorBall1x;
    std::vector <double> m_positionVectorBall1y;
    std::vector <double> m_velocityVectorBall1x;
    std::vector <double> m_velocityVectorBall1y;
    std::vector <double> m_positionVectorBall2x;
    std::vector <double> m_positionVectorBall2y;
    std::vector <double> m_velocityVectorBall2x;
    std::vector <double> m_velocityVectorBall2y;
    std::vector <double> m_angleVector1;
    std::vector <double> m_angleVector2;
    std::vector <double> m_angularVelocityVector1;
    std::vector <double> m_angularVelocityVector2;


};




