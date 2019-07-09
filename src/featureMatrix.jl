function constructMotorFeatureMatrix(motors, gearRatios)
    return hcat([
        [motor.τStall / ratio
     	 motor.ωFree * ratio
    	 motor.price
    	 motor.mass
         (motor.ωFree * ratio) / (motor.τStall / ratio)]
        for motor in motors, ratio in gearRatios]...)
end
