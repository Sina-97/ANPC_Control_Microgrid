# ANPC_Control_Microgrid
The goal of this project is to use an adaptive neural predictive controller for microgrid secondary control in Matlab Simulink.
To run this code you need to change the directory of Matlab to this folder and try to use the latest version of Matlab.

In this project, the NN Predictive Controller block was used and modified to work adaptively. We try to make the frequency and voltage
of a specific point grid stable with the help of the secondary controller in case the load change.

The microgrid is composed of Inner Loop and Outer Loop. The primary controller is responsible for the voltage and current of the inner loop, and the secondary controller is responsible for the control of the outer loop. Also, the PWM control method is used for controlling the inverter.

