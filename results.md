# BLP

Solve time: `0.42` seconds.
```
Optimal objective: 40.89
Optimal motors:
    Motor("vexMotor-393", 1.637, 10.471, 14.99, 0.0945), ratio=0.3333333333333333, link1=75.0, link2=125.0, link3=200.0
    Motor("stepperMotor-GenericNEMA14", 0.098, 139.626, 12.95, 0.12), ratio=0.05263157894736842, link1=75.0, link2=125.0, link3=200.0
    Motor("stepperMotor-GenericNEMA14", 0.098, 139.626, 12.95, 0.12), ratio=0.06666666666666667, link1=75.0, link2=125.0, link3=200.0
```

# MIQP

Solve time: `0.03` seconds.
```
Optimal objective: 40.89
Optimal motors:
    Motor("vexMotor-393", 1.637, 10.471, 14.99, 0.0945), ratio=0.3333333333333333, link1=68.0, link2=195.0, link3=137.0
    Motor("stepperMotor-GenericNEMA14", 0.098, 139.626, 12.95, 0.12), ratio=0.05263157894736842, link1=68.0, link2=195.0, link3=137.0
    Motor("stepperMotor-GenericNEMA14", 0.098, 139.626, 12.95, 0.12), ratio=0.05263157894736842, link1=68.0, link2=195.0, link3=137.0
```

# GA

Solve time: `20.47` seconds.
```
Motor 1=Motor("vexMotor-393", 1.637, 10.471, 14.99, 0.0945)
Motor 2=Motor("stepperMotor-GenericNEMA14", 0.098, 139.626, 12.95, 0.12)
Motor 3=Motor("stepperMotor-GenericNEMA14", 0.098, 139.626, 12.95, 0.12)
Link 1=96.0
Link 2=134.0
Link 3=170.0
Gear ratio 1=0.3333333333333333
Gear ratio 2=0.058823529411764705
Gear ratio 3=0.09090909090909091
Fitness=-40.89
Constraint values: [-2.56604, -0.0176998, -0.244436, -0.990333, -4.92382, -6.81092, -0.0001]
Feasible.
```
