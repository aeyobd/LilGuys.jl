using Makie
using LilGuys
import Arya

@testset "Plot xyz" begin
    # Test basic plot generation
    @testset "Basic Plot Generation" begin
        # Create sample simulation data
        pos = randn(3, 10)
        pos2 = randn(3, 5)
        
        fig = LilGuys.plot_xyz(pos)
        @test fig isa Figure

        fig = LilGuys.plot_xyz(pos, pos2)
        @test fig isa Figure
        @test size(fig.layout) == (2,2)
    end

    # Test error handling
    @testset "Error Handling" begin
        # Test with empty or invalid data
        @test_throws ArgumentError LilGuys.plot_xyz()
        @test_throws ArgumentError LilGuys.plot_xyz(randn(2, 10))
    end
    
    # Test plot content and visualization properties
    @testset "Plot Content Validation" begin
        pos = randn(3, 10)

        fig = LilGuys.plot_xyz(pos)

        # xy axis
        ax = contents(fig[1,1])[1]
        @test !isempty(plots(ax))
        p = plots(ax)[1]
        @test p[1][] ≈ Point2f.(pos[1, :], pos[2, :])
        @test p.color[] == Arya.COLORS[1]
        

        ax = contents(fig[2,2])[1]
        @test !isempty(plots(ax))
        p = plots(ax)[1]
        @test p[1][] ≈ Point2f.(pos[2, :], pos[3, :])

        ax = contents(fig[2,1])[1]
        @test !isempty(plots(ax))
        p = plots(ax)[1]
        @test p[1][] ≈ Point2f.(pos[1, :], pos[3, :])
    end
end


@testset "projecteddensity" begin
    @testset "make plot" begin
        pos = randn(3, 10)
        vel = randn(3, 10)
        snap = LilGuys.Snapshot(pos, vel, 1)
        fig, ax, p = LilGuys.projecteddensity(snap)
        @test fig isa Figure
        @test ax isa Axis
        @test p isa AbstractPlot

        # test other arguments
        fig, ax, p = LilGuys.projecteddensity(snap, bins=30)
        @test p isa AbstractPlot

        fig, ax, p = LilGuys.projecteddensity(snap, centre=false)
        @test p isa AbstractPlot

        fig, ax, p = LilGuys.projecteddensity(snap, centre=false, xdirection=3, ydirection=1)
        @test p isa AbstractPlot
    end


end
