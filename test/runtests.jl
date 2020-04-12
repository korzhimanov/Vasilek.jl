using Plots
using Test
using LinearAlgebra

@testset "Test everything" begin
    include("test_maxwell_sovers.jl")
    include("test_vlasov_sovers.jl")
end
