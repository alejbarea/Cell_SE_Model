using Makie

Makie.activate!()

function plot_data_linear_graph(x, y, dx, dy, x_label, y_label, title)
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel=x_label, ylabel=y_label, title=title)
    scatter!(ax, x, y; yerror=dy, xerror=dx)
    return fig
end