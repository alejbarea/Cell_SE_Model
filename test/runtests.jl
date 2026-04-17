using DrWatson, Test
@quickactivate "Cell_SE_Model"

# Here you include files using `srcdir`
# include(srcdir("file.jl"))
include(srcdir("cell_model.jl"))

# Run test suite
println("Starting tests")
ti = time()

@testset "Cell_SE_Model tests" begin
    simulation_parameters, domain_specs = read_parameters(datadir("sims", "parameters.toml"))

    @test simulation_parameters == SimulationParameters(1, 0.1, 10.0, 0.01, 0.01, 1.0, 0.1, 0.1, 5.0, 10.0, 1.0, 20.0)
    @test domain_specs == DomainSpecs(1000.0, 1000.0)
end

ti = time() - ti
println("\nTest took total time of:")
println(round(ti/60, digits = 3), " minutes")
