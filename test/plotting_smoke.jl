using Test
include(joinpath(@__DIR__, "..", "src", "figures", "publication_theme.jl"))

const OUTPUT = joinpath(@__DIR__, "..", "results", "log", "plotting_smoke")
mkpath(OUTPUT)

@testset "Publication rendering" begin
  for panels in (1, 4)
    fig = publication_figure(; panels, heightwidthratio=0.9)
    @test fig.scene.viewport[].widths[1] == (panels == 1 ? 234 : 504)
    ax = Axis(fig[1, 1]; xlabel="Time (days)", ylabel="Flux (mol/s)", yscale=log10)
    @test ax.xlabelpadding[] == 6
    @test ax.ylabelpadding[] == 6
    band!(ax, 1:4, [1, 2, 3, 4], [2, 3, 4, 5]; color=(:steelblue, 0.2))
    scatterlines!(ax, 1:4, [1.5, 2.5, 3.5, 4.5]; label="Median +/- SD")
    axislegend(ax)
    if panels > 1
      bars = Axis(fig[1, 2]; xticks=(1:2, ["Low", "High"]))
      stackedbars!(bars, 1:2, [0.2 0.3; 0.4 0.1];
        colors=PUBLICATION_TAB10[1:2], labels=["A", "B"])
      axislegend(bars)
      boxes = Axis(fig[2, 1])
      boxplot!(boxes, repeat(1:2; inner=4), [1, 2, 3, 4, 2, 3, 4, 8])
      hist!(boxes, [1, 2, 2, 3]; weights=fill(0.25, 4), bins=0:4, color=(:red, 0.3))
      heat = Axis(fig[2, 2])
      hm = heatmap!(heat, [1 2; 3 4]; colormap=:viridis)
      Colorbar(fig[2, 2, Right()], hm; width=6)
    end
    for ext in ("png", "pdf")
      path = joinpath(OUTPUT, "panels_$(panels).$ext")
      savefig(path, fig)
      @test filesize(path) > 1000
    end
  end
end

@testset "Plot script loading" begin
  scripts = ["ch4_performance", "community_composition", "cross_so4_thermodynamics",
    "flux_timeseries", "efficiency_growth", "robustness_figures", "dissipation_rate",
    "flux_regime_analysis", "compare_syntrophy_seeds", "rf_analysis", "h2_srb_diagnostic",
    "trajectory_diagnostics"]
  for script in scripts
    m = Module(gensym(:PlotCheck))
    Core.eval(m, :(include(path) = Base.include(@__MODULE__, path)))
    Base.include(m, joinpath(@__DIR__, "..", "src", "figures", script * ".jl"))
    @test isdefined(m, :publication_theme)
    println("Loaded ", script)
  end
end
