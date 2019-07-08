tests = ["parseConstraintsTest.jl",
         "parseMotorOptionsTest.jl"]

for test in tests
  include(test)
end
