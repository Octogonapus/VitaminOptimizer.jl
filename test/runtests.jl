tests = ["featureMatrixTest.jl"
         "optimizerTest.jl"
         "parseConstraintsTest.jl"
         "parseMotorOptionsTest.jl"]

for test in tests
  include(test)
end
