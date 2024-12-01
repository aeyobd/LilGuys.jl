@testset "Plotting Components" begin
    # Test basic plot generation
    @testset "Basic Plot Generation" begin
        # Create sample simulation data
        simulation_data = generate_sample_simulation_data()
        
        # Test trajectory plot
        @test begin
            fig = plot_trajectory(simulation_data)
            # Check that a Figure object is returned
            fig isa Figure
        end
        
        # Test energy evolution plot
        @test begin
            fig = plot_energy_evolution(simulation_data)
            fig isa Figure
        end
    end


    # Test plot customization
    @testset "Plot Customization" begin
        simulation_data = generate_sample_simulation_data()
        
        # Test custom color schemes
        @test begin
            fig = plot_trajectory(simulation_data, 
                color_scheme = :viridis,
                linewidth = 2.0
            )
            fig isa Figure
        end
        
        # Test title and label customization
        @test begin
            fig = plot_trajectory(simulation_data, 
                title = "N-Body Simulation Trajectory",
                xlabel = "X Position",
                ylabel = "Y Position"
            )
            fig isa Figure
        end
    end
    
    # Test error handling
    @testset "Error Handling" begin
        # Test with empty or invalid data
        @test_throws ArgumentError plot_trajectory([])
        @test_throws ArgumentError plot_energy_evolution(nothing)
    end
    
    # Test plot content and visualization properties
    @testset "Plot Content Validation" begin
        simulation_data = generate_sample_simulation_data()
        
        # Validate trajectory plot elements
        fig = plot_trajectory(simulation_data)
        @test begin
            # Check for specific plot elements
            axis = content(fig[1,1])
            !isempty(axis.plots)  # Ensure plots exist
        end
        
        # Check color and style consistency
        @test begin
            fig = plot_trajectory(simulation_data, color = :red)
            axis = content(fig[1,1])
            all(p.color[] == :red for p in axis.plots)
        end
    end
end
