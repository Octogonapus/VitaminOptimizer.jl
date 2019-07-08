function constructFeatureMatrix(motors, gearRatios)
    return hcat([
        [motor.τStall / ratio
     	 motor.ωFree * ratio
    	 motor.price
    	 motor.mass]
        for motor in motors, ratio in gearRatios]...)
end
